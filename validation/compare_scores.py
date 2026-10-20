"""Is the pipeline better than sorting the Open Targets data by score?

Reads leave_one_out.py's output for one disease and compares three rankings of each held-out gene among the
candidates of its own run:
- bfwalk, rwr: the pipeline's scores;
- open_targets: the gene's highest Open Targets score for the disease over every datatype except
  genetic_association (the evidence the held-out gene was hidden from). This is "sorting the Open Targets
  data by score" without the answer.
Reports AUROC and recall@k for each (and a random ranking), a paired Wilcoxon test of the pipeline against
Open Targets on the per-gene AUROC, Spearman correlations between genetic evidence, the Open Targets score,
BFWalk and RWR over the held-out genes, and the same comparison over all candidates of the full run.
Writes validation/compare_scores.log and ends with a one-line answer per pipeline method.
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from leave_one_out import VALIDATION_DIR, auroc_from_rank, format_metrics, metrics_table, ot_scores, parse_ks

logger = logging.getLogger("compare_scores")
PIPELINE_RANKINGS = {"bfwalk": "bfwalk_rank", "rwr": "rwr_rank"}
ALPHA = 0.05


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="compare_scores.py",
        description="Compare the pipeline's leave-one-out recovery of held-out genes with sorting the Open "
        "Targets data by score, and correlate genetic evidence, Open Targets, BFWalk and RWR scores.",
    )
    parser.add_argument("--details", type=Path, default=VALIDATION_DIR / "scores_left_out_details.tsv",
                        help="leave_one_out.py's per-gene details. Default: validation/scores_left_out_details.tsv.")
    parser.add_argument("--work-dir", type=Path, default=None,
                        help="leave_one_out.py's --work-dir, for the full run. Default: validation/runs/<disease>/.")
    parser.add_argument("--k", default="10,50,100,500,1000", help="Comma-separated k values for recall@k.")
    parser.add_argument("--n-random", type=int, default=1000, help="Random rankings for the p-values. Default: 1000.")
    parser.add_argument("--seed", type=int, default=0, help="Random seed. Default: 0.")
    parser.add_argument("--log", type=Path, default=VALIDATION_DIR / "compare_scores.log",
                        help="Log file. Default: validation/compare_scores.log.")
    return parser.parse_args(argv)


def paired_comparison(ranked, method):
    """Per-gene AUROC of a pipeline ranking against the Open Targets ranking: wins, losses and a two-sided
    Wilcoxon signed-rank test on the differences."""
    pipeline = np.array([auroc_from_rank(r, n) for r, n in zip(ranked[PIPELINE_RANKINGS[method]], ranked["n_candidates"])])
    open_targets = np.array([auroc_from_rank(r, n) for r, n in zip(ranked["ot_rank"], ranked["n_candidates"])])
    difference = pipeline - open_targets
    p_value = np.nan
    if np.any(difference != 0):
        p_value = float(wilcoxon(pipeline, open_targets, zero_method="wilcox").pvalue)
    return {"pipeline_auroc": float(np.mean(pipeline)), "ot_auroc": float(np.mean(open_targets)),
            "n_better": int((difference > 0).sum()), "n_worse": int((difference < 0).sum()),
            "n_tied": int((difference == 0).sum()), "p_value": p_value}


def verdict(method, comparison):
    detail = (f"mean AUROC {comparison['pipeline_auroc']:.3f} vs {comparison['ot_auroc']:.3f}; ranks the held-out "
              f"gene higher for {comparison['n_better']}, lower for {comparison['n_worse']}, tied for "
              f"{comparison['n_tied']}; Wilcoxon p = {comparison['p_value']:.3g}")
    if not comparison["p_value"] < ALPHA:
        return f"{method}: NO significant difference from sorting Open Targets by score ({detail})."
    if comparison["pipeline_auroc"] > comparison["ot_auroc"]:
        return f"{method}: YES, better than sorting Open Targets by score ({detail})."
    return f"{method}: NO, worse than sorting Open Targets by score ({detail})."


def full_run_comparison(full_dir):
    """Over every candidate of the run without a holdout: correlations of the pipeline scores with the Open
    Targets score, and the overlap of each pipeline top-n with the Open Targets top-n (n = 100, or fewer
    candidates). Genetic evidence is left out: candidates are non-seed genes, so it is 0 for all of them."""
    candidates = pd.read_csv(full_dir / "candidate_scores.tsv", sep="\t").set_index("gene")
    ot_disease_subset = pd.read_parquet(full_dir / "ot_disease_subset.parquet")
    candidates["open_targets"] = ot_scores(ot_disease_subset).reindex(candidates.index).fillna(0.0)
    columns = ["open_targets", "bfwalk_score", "rwr_score"]
    lines = [f"Full run, no holdout ({len(candidates)} candidates, {int((candidates['open_targets'] > 0).sum())} "
             f"with Open Targets evidence; genetic evidence is 0 for every candidate, so it is not compared):"]
    lines.append("Spearman correlation:\n" + candidates[columns].corr(method="spearman").to_string(float_format=lambda v: f"{v:.2f}"))
    top_n = min(100, len(candidates))
    top_ot = set(candidates["open_targets"].sort_values(ascending=False, kind="stable").head(top_n).index)
    for score in ["bfwalk_score", "rwr_score"]:
        top = set(candidates[score].sort_values(ascending=False, kind="stable").head(top_n).index)
        lines.append(f"top-{top_n} overlap {score} vs open_targets: {len(top & top_ot)} of {top_n}")
    return "\n".join(lines)


def setup_logging(log_path):
    formatter = logging.Formatter("%(message)s")
    for handler in [logging.FileHandler(log_path, mode="w"), logging.StreamHandler()]:
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    logger.setLevel(logging.INFO)


def main(argv=None):
    args = parse_args(argv)
    setup_logging(args.log)
    details = pd.read_csv(args.details, sep="\t")
    disease = details["disease"].iloc[0]
    work_dir = args.work_dir or VALIDATION_DIR / "runs" / disease
    ranked = details[details["status"] == "ranked"]
    logger.info("Disease %s: %d held-out genes, %d ranked by the pipeline (the rest are not in their run's graph "
                "or failed and are left out of every ranking below).", disease, len(details), len(ranked))
    if len(ranked) < 2:
        logger.error("Fewer than 2 ranked held-out genes - nothing to compare.")
        return 1

    rankings = {"bfwalk": ranked["bfwalk_rank"], "rwr": ranked["rwr_rank"], "open_targets": ranked["ot_rank"]}
    table = metrics_table(rankings, ranked["n_candidates"], parse_ks(args.k), args.n_random,
                          np.random.default_rng(args.seed))
    logger.info("\n== Leave-one-out recovery of the held-out genes\n%s", format_metrics(table))

    columns = ["genetic_evidence", "ot_score", "bfwalk_score", "rwr_score"]
    logger.info("\n== Spearman correlation over the held-out genes (their hidden genetic evidence, their Open "
                "Targets score without it, and their pipeline scores when held out)\n%s",
                ranked[columns].corr(method="spearman").to_string(float_format=lambda v: f"{v:.2f}"))

    full_dir = Path(work_dir) / "full"
    if (full_dir / "candidate_scores.tsv").exists():
        logger.info("\n== %s", full_run_comparison(full_dir))
    else:
        logger.info("\n(no full run in %s - skipping the all-candidate comparison)", full_dir)

    logger.info("\n== Is the pipeline better than sorting the Open Targets data by score? (paired, per held-out gene)")
    for method in PIPELINE_RANKINGS:
        logger.info(verdict(method, paired_comparison(ranked, method)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
