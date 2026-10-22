"""Is the pipeline better than sorting the Open Targets data by score?

Reads leave_one_out.py's output (one <diseaseId>_leave_one_out/ directory per disease) and, per disease,
compares three rankings of each held-out seed among the candidates of its own leave-one-out run:
- bfwalk, rwr: the pipeline's scores;
- open_targets: the gene's highest Open Targets datasource score for the disease, over every datasource
  except the genetic_association ones (the evidence that made it a seed). This is "sorting the Open Targets data
  by score" without the answer. Scores come from --association-dir (Open Targets association_by_datasource_direct).
The catalogue's disease IDs are Open Targets 26.03's; each is looked up in --disease-file (the release of
--association-dir) by ID, else by cross-reference, since 26.06 keeps many older EFO IDs only as dbXRefs.
A disease found in neither way is left out (logged), not compared against an all-zero baseline.
Reports AUROC and recall@k for each (next to a random ranking's expectation), a paired Wilcoxon test of the pipeline against
Open Targets on the per-gene AUROC, Spearman correlations between genetic evidence, the Open Targets score,
BFWalk and RWR over the held-out seeds, and the same comparison over all candidates of the full run.

Writes into each <diseaseId>_leave_one_out/: compare_scores.log, compare_scores_details.tsv (the held-out seeds
with their Open Targets score and rank) and roc_curve_vs_open_targets.png; and into --leave-one-out-dir:
compare_scores_summary.tsv (one row per disease and pipeline method, with the verdict), compare_scores.log and
roc_curve_vs_open_targets_all_diseases.png.
"""

import argparse
import logging
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.compute as pa_compute
import pyarrow.dataset as pa_dataset
import pyarrow.parquet as pa_parquet
from scipy.stats import wilcoxon

from leave_one_out import (
    LOO_SUFFIX,
    VALIDATION_DIR,
    auroc_from_rank,
    format_metrics,
    metrics_table,
    mid_rank,
    parse_ks,
    plot_roc,
)
from stage1_ingest import DATASOURCE_TO_DATATYPE

logger = logging.getLogger("compare_scores")
PIPELINE_RANKINGS = {"bfwalk": "bfwalk_rank", "rwr": "rwr_rank"}
# The held-out gene's genetic evidence is exactly what the validation hides from the baseline.
HIDDEN_DATATYPE = "genetic_association"
ALPHA = 0.05


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="compare_scores.py",
        description="Compare the pipeline's leave-one-out recovery of held-out seeds with sorting the Open "
        "Targets data by score, per disease, and correlate genetic evidence, Open Targets, BFWalk and RWR scores.",
    )
    parser.add_argument("--leave-one-out-dir", type=Path, default=VALIDATION_DIR / "leave_one_out",
                        help="leave_one_out.py's --out-dir. Default: validation/leave_one_out/.")
    parser.add_argument("--association-dir", type=Path, required=True,
                        help="Open Targets association_by_datasource_direct parquet directory, e.g. "
                        "data/opentargets/26.06/association_by_datasource_direct.")
    parser.add_argument("--disease-file", type=Path, required=True,
                        help="Open Targets disease parquet of the same release as --association-dir, to translate "
                        "the catalogue's disease IDs, e.g. data/opentargets/26.06/disease/disease.parquet.")
    parser.add_argument("--diseases", default=None,
                        help="Only these diseases, comma-separated Open Targets IDs, e.g. EFO_0003767,MONDO_0004975.")
    parser.add_argument("--n-diseases", type=int, default=None, help="Only the first N diseases (sorted by ID).")
    parser.add_argument("--k", default="10,50,100,500,1000", help="Comma-separated k values for recall@k.")
    return parser.parse_args(argv)


def select_diseases(loo_dir, diseases, n_diseases):
    available = sorted(path.name[: -len(LOO_SUFFIX)] for path in loo_dir.glob(f"*{LOO_SUFFIX}")
                       if (path / "scores_left_out_details.tsv").exists())
    if not available:
        raise SystemExit(f"No *{LOO_SUFFIX} results in {loo_dir} - run leave_one_out.py first.")
    if diseases:
        requested = [disease.strip() for disease in diseases.split(",") if disease.strip()]
        missing = sorted(set(requested) - set(available))
        if missing:
            logger.warning("No leave-one-out results for: %s", ", ".join(missing))
        available = [disease for disease in available if disease in requested]
    if n_diseases is not None:
        if n_diseases < 1:
            raise SystemExit(f"--n-diseases must be at least 1 (got {n_diseases}).")
        available = available[:n_diseases]
    return available


def map_disease_ids(diseases, disease_file):
    """{catalogue disease ID: its ID in disease_file, or None}, by Stage 1's rule (resolve_disease_id): the same
    ID, else the first disease listing it as a cross-reference ("EFO_0003767" -> xref "EFO:0003767")."""
    disease_table = pa_parquet.read_table(disease_file, columns=["id", "dbXRefs"]).to_pandas()
    known = set(disease_table["id"])
    by_xref = {}
    for disease_id, xrefs in zip(disease_table["id"], disease_table["dbXRefs"]):
        if xrefs is not None:
            for xref in xrefs:
                by_xref.setdefault(xref, disease_id)
    mapped = {}
    for disease in diseases:
        if disease in known:
            mapped[disease] = disease
        else:
            mapped[disease] = by_xref.get(re.sub(r"^([A-Za-z]+)_", r"\1:", disease))
    return mapped


def ot_baselines(association_dir, diseases):
    """{disease: Series gene -> highest Open Targets score over every datatype except the hidden one}."""
    dataset = pa_dataset.dataset(sorted(association_dir.glob("*.parquet")), format="parquet")
    associations = dataset.to_table(
        columns=["targetId", "diseaseId", "aggregationValue", "associationScore"],
        filter=pa_compute.field("diseaseId").isin(diseases),
    ).to_pandas()
    unknown = sorted(set(associations["aggregationValue"]) - set(DATASOURCE_TO_DATATYPE))
    if unknown:
        logger.warning("Unknown Open Targets datasource(s) %s - their rows are left out.", unknown)
    associations["datatype_id"] = associations["aggregationValue"].map(DATASOURCE_TO_DATATYPE)
    visible = associations[associations["datatype_id"].notna() & (associations["datatype_id"] != HIDDEN_DATATYPE)]
    baselines = {}
    for disease, rows in visible.groupby("diseaseId"):
        baselines[disease] = rows.groupby("targetId")["associationScore"].max()
    return baselines


def add_open_targets_ranks(ranked, disease_dir, baseline):
    """The Open Targets score and rank of each held-out seed among the candidates of its own run."""
    ot_score, ot_rank = [], []
    for gene in ranked["gene"]:
        candidates = pd.read_csv(disease_dir / gene / "candidate_scores.tsv", sep="\t", index_col="gene").index
        pool = baseline.reindex(candidates).fillna(0.0)
        ot_score.append(float(pool[gene]))
        ot_rank.append(mid_rank(pool, gene))
    return ranked.assign(ot_score=ot_score, ot_rank=ot_rank)


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


def verdict(comparison):
    if not comparison["p_value"] < ALPHA:
        return "no_difference"
    if comparison["pipeline_auroc"] > comparison["ot_auroc"]:
        return "better"
    return "worse"


def verdict_text(method, comparison):
    detail = (f"mean AUROC {comparison['pipeline_auroc']:.3f} vs {comparison['ot_auroc']:.3f}; ranks the held-out "
              f"gene higher for {comparison['n_better']}, lower for {comparison['n_worse']}, tied for "
              f"{comparison['n_tied']}; Wilcoxon p = {comparison['p_value']:.3g}")
    answer = {"better": "YES, better than", "worse": "NO, worse than",
              "no_difference": "NO significant difference from"}[verdict(comparison)]
    return f"{method}: {answer} sorting Open Targets by score ({detail})."


def full_run_comparison(full_dir, baseline):
    """Over every candidate of the run without a holdout: correlations of the pipeline scores with the Open
    Targets score, and the overlap of each pipeline top-n with the Open Targets top-n (n = 100, or fewer
    candidates). Genetic evidence is left out: candidates are non-seed genes, so they have none."""
    candidates = pd.read_csv(full_dir / "candidate_scores.tsv", sep="\t").set_index("gene")
    candidates["open_targets"] = baseline.reindex(candidates.index).fillna(0.0)
    columns = ["open_targets", "bfwalk_score", "rwr_score"]
    lines = [f"Full run, no holdout ({len(candidates)} candidates, {int((candidates['open_targets'] > 0).sum())} "
             f"with Open Targets evidence other than genetic association):"]
    lines.append("Spearman correlation:\n" + candidates[columns].corr(method="spearman").to_string(float_format=lambda v: f"{v:.2f}"))
    top_n = min(100, len(candidates))
    top_ot = set(candidates["open_targets"].sort_values(ascending=False, kind="stable").head(top_n).index)
    for score in ["bfwalk_score", "rwr_score"]:
        top = set(candidates[score].sort_values(ascending=False, kind="stable").head(top_n).index)
        lines.append(f"top-{top_n} overlap {score} vs open_targets: {len(top & top_ot)} of {top_n}")
    return "\n".join(lines)


def disease_logger(log_path):
    child = logging.getLogger(f"compare_scores.{log_path.parent.name}")
    for handler in list(child.handlers):
        child.removeHandler(handler)
    handler = logging.FileHandler(log_path, mode="w")
    handler.setFormatter(logging.Formatter("%(message)s"))
    child.addHandler(handler)
    return child


def compare_disease(disease, ot_disease, disease_dir, baseline, ks):
    """Returns (summary rows, ranked seeds with Open Targets columns) or None."""
    log = disease_logger(disease_dir / "compare_scores.log")
    details = pd.read_csv(disease_dir / "scores_left_out_details.tsv", sep="\t")
    ranked = details[details["status"] == "ranked"]
    if baseline.empty:
        log.warning("%s: no Open Targets evidence other than genetic association in --association-dir - every "
                    "candidate's baseline score is 0.", disease)
    log.info("Disease %s (Open Targets %s): %d held-out seeds, %d ranked.", disease, ot_disease, len(details),
             len(ranked))
    if len(ranked) < 2:
        log.warning("Fewer than 2 ranked held-out seeds - nothing to compare.")
        return None
    ranked = add_open_targets_ranks(ranked, disease_dir, baseline)
    ranked.to_csv(disease_dir / "compare_scores_details.tsv", sep="\t", index=False, na_rep="NA")

    rankings = {"bfwalk": ranked["bfwalk_rank"], "rwr": ranked["rwr_rank"], "open_targets": ranked["ot_rank"]}
    table = metrics_table(rankings, ranked["n_candidates"], ks)
    log.info("\n== Leave-one-out recovery of the held-out seeds\n%s", format_metrics(table))
    columns = ["genetic_evidence", "ot_score", "bfwalk_score", "rwr_score"]
    log.info("\n== Spearman correlation over the held-out seeds (their genetic evidence, their Open Targets "
             "score without it, and their pipeline scores when held out)\n%s",
             ranked[columns].corr(method="spearman").to_string(float_format=lambda v: f"{v:.2f}"))
    if (disease_dir / "full" / "candidate_scores.tsv").exists():
        log.info("\n== %s", full_run_comparison(disease_dir / "full", baseline))
    plot_roc(rankings, ranked["n_candidates"], disease_dir / "roc_curve_vs_open_targets.png",
             f"{disease}: leave-one-out vs Open Targets, {len(ranked)} held-out seeds")

    log.info("\n== Is the pipeline better than sorting the Open Targets data by score? (paired, per held-out seed)")
    summary = []
    for method in PIPELINE_RANKINGS:
        comparison = paired_comparison(ranked, method)
        log.info(verdict_text(method, comparison))
        summary.append({"diseaseId": disease, "ot_diseaseId": ot_disease, "method": method, "n_held_out": len(ranked), **comparison,
                        "verdict": verdict(comparison)})
    return summary, ranked.assign(diseaseId=disease)


def setup_logging(log_path):
    formatter = logging.Formatter("%(message)s")
    for handler in [logging.FileHandler(log_path, mode="w"), logging.StreamHandler()]:
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    logger.setLevel(logging.INFO)


def main(argv=None):
    args = parse_args(argv)
    loo_dir = args.leave_one_out_dir
    setup_logging(loo_dir / "compare_scores.log")
    ks = parse_ks(args.k)
    diseases = select_diseases(loo_dir, args.diseases, args.n_diseases)
    ot_ids = map_disease_ids(diseases, args.disease_file)
    not_found = [disease for disease in diseases if ot_ids[disease] is None]
    if not_found:
        logger.warning("%d disease(s) not in %s, by ID or cross-reference - left out: %s", len(not_found),
                       args.disease_file, ", ".join(not_found))
    diseases = [disease for disease in diseases if ot_ids[disease] is not None]
    n_translated = sum(ot_ids[disease] != disease for disease in diseases)
    logger.info("Diseases: %d (%d found under another ID by cross-reference). Reading Open Targets associations "
                "from %s ...", len(diseases), n_translated, args.association_dir)
    baselines = ot_baselines(args.association_dir, sorted({ot_ids[disease] for disease in diseases}))

    summaries, all_ranked = [], []
    for disease in diseases:
        result = compare_disease(disease, ot_ids[disease], loo_dir / f"{disease}{LOO_SUFFIX}",
                                 baselines.get(ot_ids[disease], pd.Series(dtype=float)), ks)
        if result is not None:
            summaries.extend(result[0])
            all_ranked.append(result[1])
    if not summaries:
        logger.error("No disease had 2 or more ranked held-out seeds - nothing to compare.")
        return 1

    summary = pd.DataFrame(summaries)
    summary.to_csv(loo_dir / "compare_scores_summary.tsv", sep="\t", index=False)
    pooled = pd.concat(all_ranked, ignore_index=True)
    plot_roc({"bfwalk": pooled["bfwalk_rank"], "rwr": pooled["rwr_rank"], "open_targets": pooled["ot_rank"]},
             pooled["n_candidates"], loo_dir / "roc_curve_vs_open_targets_all_diseases.png",
             f"Leave-one-out vs Open Targets, {pooled['diseaseId'].nunique()} diseases, {len(pooled)} seeds")
    logger.info("\n== Is the pipeline better than sorting the Open Targets data by score? Diseases per answer "
                "(paired Wilcoxon per disease, p < %g):", ALPHA)
    counts = summary.pivot_table(index="method", columns="verdict", values="diseaseId", aggfunc="count", fill_value=0)
    logger.info(counts.reindex(columns=["better", "no_difference", "worse"], fill_value=0).to_string())
    means = summary.groupby("method")[["pipeline_auroc", "ot_auroc"]].mean()
    logger.info("\nMean AUROC over diseases:\n%s", means.to_string(float_format=lambda v: f"{v:.3f}"))
    logger.info("\nWrote %s and %s.", loo_dir / "compare_scores_summary.tsv",
                loo_dir / "roc_curve_vs_open_targets_all_diseases.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
