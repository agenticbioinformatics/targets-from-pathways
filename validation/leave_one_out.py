"""Leave-one-out validation of the pipeline for one disease.

For each held-out gene (by default every seed gene that is a node of the disease's graph), the pipeline's
Stages 2-5 are rerun with that gene treated as not associated with the disease: Stage 2's
--benchmark-holdout-file removes it from the seed genes and the GSEA signature, and Stage 4 ignores its
Open Targets evidence. The gene is then a candidate like any other, and its score and rank among that run's
candidates measure whether the pipeline recovers it.

Writes to --out-dir (default validation/):
- scores_left_out.tsv          gene, score (the --method score of each held-out gene in its own run)
- scores_left_out_details.tsv  per gene: genetic evidence, both pipeline scores and ranks, the Open Targets
                               baseline score and rank, and the candidate pool size (read by compare_scores.py)
- leave_one_out.log            AUROC and recall@k of the pipeline next to a random ranking
Runs go to --work-dir (default validation/runs/<disease>/): Stage 1 once, one full run without a holdout, and
one directory per held-out gene, deleted after scoring (kept when a stage fails, with its pipeline.log).

Metrics, over the held-out genes that are candidates in their own run:
- rank: mid-rank among the run's candidates (ties count half), 1 = best.
- AUROC: per gene, the fraction of the other candidates it outscores (ties count half); averaged over genes.
- recall@k: the fraction of held-out genes ranked within the top k.
- random: a uniformly random ranking of the same candidates; expectation, plus a permutation p-value of the
  pipeline's metric against --n-random random rankings.
"""

import argparse
import json
import logging
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

REPO_DIR = Path(__file__).resolve().parent.parent
PIPELINE_DIR = REPO_DIR / "pipeline"
VALIDATION_DIR = Path(__file__).resolve().parent

METHODS = ["bfwalk", "rwr"]
ENSEMBL_GENE_ID = re.compile(r"^ENSG\d{11}$")
# The Open Targets baseline uses every datatype except this one: the held-out gene's genetic evidence is
# exactly what the validation hides.
HIDDEN_DATATYPE = "genetic_association"
RECORD_COLUMNS = [
    "gene", "symbol", "status", "genetic_evidence", "ot_score", "n_candidates",
    "bfwalk_score", "bfwalk_rank", "rwr_score", "rwr_rank", "ot_rank",
]

logger = logging.getLogger("leave_one_out")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="leave_one_out.py",
        description="Leave-one-out validation for one disease: rerun Stages 2-5 once per held-out gene with "
        "that gene treated as not associated with the disease, and score its recovery (AUROC, recall@k) "
        "against a random ranking.",
    )
    parser.add_argument("--disease", required=True, help="Disease ID, e.g. MONDO_0005265.")
    parser.add_argument("--data-dir", type=Path, default=REPO_DIR / "data", help="Stage 1 source cache. Default: data/.")
    parser.add_argument("--from-manifest", type=Path, default=None,
                        help="Use this existing Stage 1 manifest.json for --disease instead of running Stage 1.")
    parser.add_argument("--allow-low-coverage", action="store_true", help="Passed to Stage 1.")
    parser.add_argument("--genes", type=Path, default=None,
                        help="Held-out genes, one Ensembl ID per line ('#' comments allowed). Default: every seed "
                        "gene that is a node of the disease's graph.")
    parser.add_argument("--max-genes", type=int, default=None,
                        help="Hold out at most this many genes, a random sample (--seed). Default: all.")
    parser.add_argument("--method", choices=METHODS, default="bfwalk",
                        help="Score written to scores_left_out.tsv. Both are computed. Default: bfwalk.")
    parser.add_argument("--k", default="10,50,100,500,1000", help="Comma-separated k values for recall@k.")
    parser.add_argument("--n-random", type=int, default=1000, help="Random rankings for the p-values. Default: 1000.")
    parser.add_argument("--seed", type=int, default=0, help="Random seed (gene sample, random rankings). Default: 0.")
    parser.add_argument("--jobs", type=int, default=1, help="Held-out genes processed in parallel. Default: 1.")
    parser.add_argument("--work-dir", type=Path, default=None,
                        help="Pipeline runs. Default: validation/runs/<disease>/.")
    parser.add_argument("--out-dir", type=Path, default=VALIDATION_DIR, help="Output files. Default: validation/.")
    return parser.parse_args(argv)


# ---------------------------------------------------------------------------------------------------------
# Metrics (shared with compare_scores.py)
# ---------------------------------------------------------------------------------------------------------


def mid_rank(scores, gene):
    """Rank of gene among scores (a Series indexed by gene), 1 = best; tied genes share half the places."""
    value = scores[gene]
    n_greater = int((scores > value).sum())
    n_tied_others = int((scores == value).sum()) - 1
    return 1.0 + n_greater + 0.5 * n_tied_others


def auroc_from_rank(rank, n_candidates):
    """Fraction of the other n_candidates - 1 genes that the gene outscores, ties counting half."""
    if n_candidates <= 1:
        return np.nan
    return (n_candidates - rank) / (n_candidates - 1)


def ranking_metrics(ranks, n_candidates, ks):
    """Mean AUROC and recall@k for one ranking, given each held-out gene's rank and candidate pool size."""
    ranks = np.asarray(ranks, dtype=float)
    n_candidates = np.asarray(n_candidates, dtype=float)
    metrics = {"auroc": float(np.nanmean([auroc_from_rank(r, n) for r, n in zip(ranks, n_candidates)]))}
    for k in ks:
        metrics[f"recall@{k}"] = float(np.mean(ranks <= k))
    return metrics


def random_metrics(n_candidates, ks, n_random, rng):
    """Expected metrics of a uniformly random ranking, and n_random sampled values of each metric."""
    n_candidates = np.asarray(n_candidates, dtype=int)
    expected = {"auroc": 0.5}
    for k in ks:
        expected[f"recall@{k}"] = float(np.mean(np.minimum(k, n_candidates) / n_candidates))
    # A random ranking puts each held-out gene at a uniform rank among its run's candidates.
    sampled_ranks = rng.integers(1, n_candidates + 1, size=(n_random, len(n_candidates)))
    samples = {"auroc": np.nanmean((n_candidates - sampled_ranks) / np.maximum(n_candidates - 1, 1), axis=1)}
    for k in ks:
        samples[f"recall@{k}"] = np.mean(sampled_ranks <= k, axis=1)
    return expected, samples


def permutation_p_value(observed, samples):
    """One-sided: how often a random ranking does at least as well."""
    return (1 + int(np.sum(samples >= observed - 1e-12))) / (1 + len(samples))


def metrics_table(rankings, n_candidates, ks, n_random, rng):
    """rankings: {name: per-gene ranks}. Returns a DataFrame of metric rows, with the random ranking's
    expectation and each ranking's permutation p-value against it."""
    expected, samples = random_metrics(n_candidates, ks, n_random, rng)
    rows = []
    for name, ranks in rankings.items():
        observed = ranking_metrics(ranks, n_candidates, ks)
        for metric, value in observed.items():
            rows.append({"ranking": name, "metric": metric, "value": value,
                         "random_expected": expected[metric],
                         "p_vs_random": permutation_p_value(value, samples[metric])})
    return pd.DataFrame(rows)


def format_metrics(table):
    wide = table.pivot(index="metric", columns="ranking", values="value")
    order = [m for m in ["auroc"] + [c for c in table["metric"].unique() if c.startswith("recall@")]]
    wide = wide.loc[order]
    wide["random (expected)"] = table.drop_duplicates("metric").set_index("metric").loc[order, "random_expected"]
    lines = [wide.to_string(float_format=lambda v: f"{v:.3f}")]
    p_values = table.pivot(index="metric", columns="ranking", values="p_vs_random").loc[order]
    lines.append("\npermutation p-value vs random (one-sided, chance a random ranking does at least as well):")
    lines.append(p_values.to_string(float_format=lambda v: f"{v:.4f}"))
    return "\n".join(lines)


def parse_ks(raw):
    ks = []
    for value in raw.split(","):
        if value.strip():
            ks.append(int(value))
    return ks


# ---------------------------------------------------------------------------------------------------------
# Pipeline runs
# ---------------------------------------------------------------------------------------------------------


def run_stage(script, args, log_path):
    command = [sys.executable, str(PIPELINE_DIR / script)] + [str(arg) for arg in args]
    with open(log_path, "a") as log:
        log.write(f"$ {' '.join(command)}\n")
        log.flush()
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        raise RuntimeError(f"{script} failed with exit code {result.returncode} (see {log_path})")


def stage1_files(manifest_path):
    manifest = json.loads(manifest_path.read_text())
    names = ["manifest.json"]
    for artifact in manifest["output_artifacts"]:
        names.append(Path(artifact["path"]).name)
    return names


def prepare_stage1(args, work_dir):
    """Return the Stage 1 directory for args.disease: --from-manifest's, an earlier run's, or a new one."""
    if args.from_manifest is not None:
        stage1_dir = args.from_manifest.resolve().parent
    else:
        stage1_dir = work_dir / "stage1"
        manifest_path = stage1_dir / "manifest.json"
        if manifest_path.exists():
            logger.info("Reusing Stage 1 output in %s.", stage1_dir)
        else:
            stage1_dir.mkdir(parents=True, exist_ok=True)
            logger.info("Running Stage 1 for %s (log: %s).", args.disease, stage1_dir / "stage1.log")
            stage1_args = ["--disease", args.disease, "--data-dir", args.data_dir, "--out-dir", stage1_dir]
            if args.allow_low_coverage:
                stage1_args.append("--allow-low-coverage")
            run_stage("stage1_ingest.py", stage1_args, stage1_dir / "stage1.log")
    manifest = json.loads((stage1_dir / "manifest.json").read_text())
    resolved = manifest["disease"]["efo_id"]
    requested = manifest["cli_parameters"].get("disease", resolved)
    if args.disease not in (resolved, requested):
        raise SystemExit(f"{stage1_dir / 'manifest.json'} is for {resolved}, not --disease {args.disease}.")
    logger.info("Disease: %s (%s).", manifest["disease"]["name"], resolved)
    return stage1_dir


def run_pipeline_copy(stage1_dir, run_dir, holdout_gene, bfwalk_threads):
    """Copy Stage 1's output into run_dir and run Stages 2-5 there, optionally holding out one gene."""
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True)
    for name in stage1_files(stage1_dir / "manifest.json"):
        shutil.copy(stage1_dir / name, run_dir / name)
    manifest = run_dir / "manifest.json"
    log_path = run_dir / "pipeline.log"
    stage2_args = ["--manifest", manifest]
    if holdout_gene is not None:
        (run_dir / "holdout.txt").write_text(holdout_gene + "\n")
        stage2_args += ["--benchmark-holdout-file", run_dir / "holdout.txt"]
    run_stage("stage2_gsea_discovery.py", stage2_args, log_path)
    run_stage("stage3_build_graph.py", ["--manifest", manifest], log_path)
    run_stage("stage4_genetic_evidence_weights.py", ["--manifest", manifest], log_path)
    run_stage("stage5_score_candidates.py", ["--graph-dir", run_dir, "--method", ",".join(METHODS),
                                             "--bfwalk-threads", bfwalk_threads], log_path)
    return pd.read_csv(run_dir / "candidate_scores.tsv", sep="\t").set_index("gene")


def ot_scores(ot_disease_subset):
    """Open Targets baseline: each gene's highest score for the disease over every datatype except the hidden
    one. Genes without such evidence score 0."""
    visible = ot_disease_subset[ot_disease_subset["datatype_id"] != HIDDEN_DATATYPE]
    return visible.groupby("gene_id")["score"].max()


def genetic_evidence(ot_disease_subset):
    hidden = ot_disease_subset[ot_disease_subset["datatype_id"] == HIDDEN_DATATYPE]
    return hidden.groupby("gene_id")["score"].max()


def score_held_out_gene(gene, candidates, baseline):
    """Record for one held-out gene from its run's candidate_scores (indexed by gene)."""
    record = {"gene": gene, "n_candidates": len(candidates)}
    if gene not in candidates.index:
        # Its pathways left the graph with it, so it is no candidate: the pipeline cannot rank it.
        record["status"] = "not_in_graph"
        return record
    record["status"] = "ranked"
    for method in METHODS:
        scores = candidates[f"{method}_score"]
        record[f"{method}_score"] = float(scores[gene])
        record[f"{method}_rank"] = mid_rank(scores, gene)
    pool_baseline = baseline.reindex(candidates.index).fillna(0.0)
    record["ot_rank"] = mid_rank(pool_baseline, gene)
    return record


def held_out_genes(args, full_dir, rng):
    if args.genes is not None:
        genes = []
        for line in args.genes.read_text().splitlines():
            entry = line.split("#")[0].strip()
            if entry:
                genes.append(entry)
        not_ensembl = [gene for gene in genes if not ENSEMBL_GENE_ID.match(gene)]
        if not_ensembl:
            raise SystemExit(f"--genes {args.genes}: not Ensembl gene IDs: {not_ensembl[:10]}")
        source = str(args.genes)
    else:
        seeds = pd.read_csv(full_dir / "disease_seed_genes.tsv", sep="\t")["gene_id"]
        graph_genes = set(json.loads((full_dir / "graph_gene_index.json").read_text()))
        genes = sorted(set(seeds) & graph_genes)
        n_outside = len(set(seeds) - graph_genes)
        source = f"seed genes in the graph ({n_outside} of {len(seeds)} seeds are on no disease-relevant pathway " \
                 f"and cannot be ranked by the pipeline at all)"
    genes = sorted(set(genes))
    if args.max_genes is not None and len(genes) > args.max_genes:
        genes = sorted(rng.choice(genes, args.max_genes, replace=False))
        source += f", random sample of {args.max_genes} (seed {args.seed})"
    logger.info("Held-out genes: %d - %s.", len(genes), source)
    return genes


def run_leave_one_out(args, stage1_dir, work_dir, genes, baseline, evidence, symbols):
    bfwalk_threads = max(1, (os.cpu_count() or 1) // max(1, args.jobs))

    def one(gene):
        run_dir = work_dir / "loo" / gene
        try:
            candidates = run_pipeline_copy(stage1_dir, run_dir, gene, bfwalk_threads)
        except RuntimeError as error:
            return {"gene": gene, "status": f"failed: {error}"}
        shutil.rmtree(run_dir)
        return score_held_out_gene(gene, candidates, baseline)

    records = {}
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        futures = [pool.submit(one, gene) for gene in genes]
        for n_done, future in enumerate(as_completed(futures), start=1):
            record = future.result()
            record["symbol"] = symbols.get(record["gene"], record["gene"])
            record["genetic_evidence"] = float(evidence.get(record["gene"], 0.0))
            record["ot_score"] = float(baseline.get(record["gene"], 0.0))
            records[record["gene"]] = record
            rank_text = ""
            if record["status"] == "ranked":
                rank_text = f"{args.method} rank {record[f'{args.method}_rank']:.0f}/{record['n_candidates']}"
            logger.info("[%d/%d] %s (%s): %s %s", n_done, len(genes), record["gene"], record["symbol"],
                        record["status"], rank_text)
    return pd.DataFrame([records[gene] for gene in genes]).reindex(columns=RECORD_COLUMNS)


def setup_logging(log_path):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s", "%Y-%m-%d %H:%M:%S")
    for handler in [logging.FileHandler(log_path, mode="w"), logging.StreamHandler()]:
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    logger.setLevel(logging.INFO)


def main(argv=None):
    args = parse_args(argv)
    work_dir = (args.work_dir or VALIDATION_DIR / "runs" / args.disease).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(args.out_dir / "leave_one_out.log")
    ks = parse_ks(args.k)
    rng = np.random.default_rng(args.seed)
    logger.info("leave_one_out.py %s", " ".join(sys.argv[1:] if argv is None else argv))

    stage1_dir = prepare_stage1(args, work_dir)
    full_dir = work_dir / "full"
    if (full_dir / "candidate_scores.tsv").exists():
        logger.info("Reusing the full run (no holdout) in %s.", full_dir)
    else:
        logger.info("Full run (no holdout) in %s.", full_dir)
        run_pipeline_copy(stage1_dir, full_dir, None, os.cpu_count() or 1)

    ot_disease_subset = pd.read_parquet(full_dir / "ot_disease_subset.parquet")
    baseline = ot_scores(ot_disease_subset)
    evidence = genetic_evidence(ot_disease_subset)
    genes_df = pd.read_parquet(full_dir / "genes.parquet", columns=["gene_id", "symbol"])
    symbols = dict(zip(genes_df["gene_id"], genes_df["symbol"]))

    genes = held_out_genes(args, full_dir, rng)
    records = run_leave_one_out(args, stage1_dir, work_dir, genes, baseline, evidence, symbols)

    scores = records[["gene", f"{args.method}_score"]].rename(columns={f"{args.method}_score": "score"})
    scores.to_csv(args.out_dir / "scores_left_out.tsv", sep="\t", index=False, na_rep="NA")
    records.insert(0, "disease", args.disease)
    records.to_csv(args.out_dir / "scores_left_out_details.tsv", sep="\t", index=False, na_rep="NA")
    logger.info("Wrote %s and %s.", args.out_dir / "scores_left_out.tsv", args.out_dir / "scores_left_out_details.tsv")

    ranked = records[records["status"] == "ranked"]
    n_failed = int(records["status"].str.startswith("failed").sum())
    logger.info("Held-out genes ranked: %d of %d (%d not in their run's graph, %d failed).", len(ranked),
                len(records), int((records["status"] == "not_in_graph").sum()), n_failed)
    if ranked.empty:
        logger.error("No held-out gene could be ranked - no metrics.")
        return 1
    table = metrics_table({method: ranked[f"{method}_rank"] for method in METHODS},
                          ranked["n_candidates"], ks, args.n_random, rng)
    logger.info("Leave-one-out metrics (%d held-out genes; median candidate pool %d genes; %s scores the output):\n%s",
                len(ranked), int(ranked["n_candidates"].median()), args.method, format_metrics(table))
    return 0


if __name__ == "__main__":
    sys.exit(main())
