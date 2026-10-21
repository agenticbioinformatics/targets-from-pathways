"""Leave-one-out validation over the disease graph catalogue (validation/build_graph_catalogue.py).

For each disease, its graph and seeds come from <catalogue>/<diseaseId>_disease_graph/. Every seed that is a
graph node is held out in turn: it is removed from the seeds, every graph gene that is not one of the remaining
seeds is a candidate (the held-out gene among them), and the candidates are scored with Stage 5's BFWalk and
random walk with restart (pipeline/stage5_score_candidates.py). The held-out gene's rank among the candidates
measures whether the pipeline recovers it. The graph itself is not rebuilt: its pathways were selected by a GSEA
that included the held-out gene, so this measures recovery given the disease's pathways.

Writes into --out-dir/<diseaseId>_leave_one_out/:
- full/candidate_scores.tsv       all seeds, nothing held out (read by compare_scores.py)
- <gene>/candidate_scores.tsv     every candidate's bfwalk_score and rwr_score with <gene> held out
- <gene>/result.json, seed.log    the held-out gene's scores and ranks, and what was run
- scores_left_out_details.tsv     one row per held-out seed (read by compare_scores.py)
- metrics.tsv, leave_one_out.log  AUROC and recall@k of BFWalk and RWR next to a random ranking
- roc_curve.png                   ROC curves of BFWalk, RWR and a random ranking
and into --out-dir: leave_one_out_summary.tsv (one row per disease and method), leave_one_out.log, and
roc_curve_all_diseases.png (every held-out seed of every disease pooled).

Metrics, over the held-out seeds:
- rank: mid-rank among the candidates (ties count half), 1 = best.
- AUROC: per gene, the fraction of the other candidates it outscores (ties count half); averaged over genes.
  The ROC curve plots, for each fraction x of the candidates taken from the top, the fraction of held-out genes
  ranked within it; its area is this AUROC.
- recall@k: the fraction of held-out genes ranked within the top k.
- random: a uniformly random ranking of the same candidates; expectation, plus a permutation p-value of each
  method's metric against --n-random random rankings.
"""

import argparse
import json
import logging
import os
import signal
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

VALIDATION_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(VALIDATION_DIR.parent / "pipeline"))
import stage5_score_candidates as sc  # noqa: E402
from stage0_schemas import SeedGenesSchema, load_graph  # noqa: E402

METHODS = ["bfwalk", "rwr"]
GRAPH_SUFFIX = "_disease_graph"
LOO_SUFFIX = "_leave_one_out"
RECORD_COLUMNS = ["gene", "status", "genetic_evidence", "n_seeds", "n_candidates",
                  "bfwalk_score", "bfwalk_rank", "rwr_score", "rwr_rank"]

logger = logging.getLogger("leave_one_out")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="leave_one_out.py",
        description="Leave-one-out validation over the disease graph catalogue: hold out each seed of a disease "
        "in turn, score the candidates with BFWalk and RWR, and measure the held-out gene's recovery (AUROC, "
        "recall@k) against a random ranking.",
    )
    parser.add_argument("--catalogue", type=Path, required=True,
                        help="build_graph_catalogue.py's output, e.g. "
                        "validation/Disease_graph_catalogue/disease_graph_catalogue.")
    parser.add_argument("--out-dir", type=Path, default=VALIDATION_DIR / "leave_one_out",
                        help="Output directory. Default: validation/leave_one_out/.")
    parser.add_argument("--diseases", default=None,
                        help="Only these diseases, comma-separated Open Targets IDs, e.g. EFO_0003767,MONDO_0004975.")
    parser.add_argument("--n-diseases", type=int, default=None,
                        help="Only the first N diseases of the catalogue (sorted by ID).")
    parser.add_argument("--max-seeds", type=int, default=None,
                        help="Hold out at most this many seeds per disease, a random sample (--seed). Default: all.")
    parser.add_argument("--k", default="10,50,100,500,1000", help="Comma-separated k values for recall@k.")
    parser.add_argument("--n-random", type=int, default=1000, help="Random rankings for the p-values. Default: 1000.")
    parser.add_argument("--seed", type=int, default=0, help="Random seed (seed sample, random rankings). Default: 0.")
    parser.add_argument("--restart-prob", type=float, default=0.5, help="RWR restart probability. Default: 0.5.")
    parser.add_argument("--bfwalk-alpha", type=float, default=0.5, help="BFWalk attenuation. Default: 0.5.")
    parser.add_argument("--jobs", type=int, default=1, help="Held-out seeds scored in parallel. Default: 1.")
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


def random_ranks(n_candidates, size, rng):
    """size random rankings: each held-out gene at a uniform rank among its own candidates."""
    n_candidates = np.asarray(n_candidates, dtype=int)
    return rng.integers(1, n_candidates + 1, size=(size, len(n_candidates)))


def random_metrics(n_candidates, ks, n_random, rng):
    """Expected metrics of a uniformly random ranking, and n_random sampled values of each metric."""
    n_candidates = np.asarray(n_candidates, dtype=int)
    expected = {"auroc": 0.5}
    for k in ks:
        expected[f"recall@{k}"] = float(np.mean(np.minimum(k, n_candidates) / n_candidates))
    sampled_ranks = random_ranks(n_candidates, n_random, rng)
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


def plot_roc(rankings, n_candidates, path, title, rng):
    """ROC curve of each ranking: x = fraction of the other candidates ranked above the held-out gene, y = the
    fraction of held-out genes recovered within that top fraction. Area under it = the mean AUROC."""
    n_candidates = np.asarray(n_candidates, dtype=float)
    rankings = dict(rankings)
    rankings["random"] = random_ranks(n_candidates, 1, rng)[0]
    figure, axis = plt.subplots(figsize=(5, 5))
    for name, ranks in rankings.items():
        false_positive_fraction = (np.asarray(ranks, dtype=float) - 1) / np.maximum(n_candidates - 1, 1)
        x = np.concatenate([[0.0], np.sort(false_positive_fraction), [1.0]])
        y = np.concatenate([[0.0], np.arange(1, len(false_positive_fraction) + 1) / len(false_positive_fraction), [1.0]])
        auroc = 1.0 - float(np.mean(false_positive_fraction))
        style = {"linestyle": ":", "color": "grey"} if name == "random" else {}
        axis.step(x, y, where="post", label=f"{name} (AUROC {auroc:.3f})", **style)
    axis.plot([0, 1], [0, 1], linestyle="--", color="lightgrey", linewidth=1)
    axis.set_xlabel("fraction of candidates ranked above the held-out gene")
    axis.set_ylabel("fraction of held-out genes recovered")
    axis.set_title(title, fontsize=10)
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.legend(loc="lower right", fontsize=9)
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    plt.close(figure)


# ---------------------------------------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------------------------------------


def select_diseases(catalogue, diseases, n_diseases):
    available = sorted(path.name[: -len(GRAPH_SUFFIX)] for path in catalogue.glob(f"*{GRAPH_SUFFIX}") if path.is_dir())
    if not available:
        raise SystemExit(f"No *{GRAPH_SUFFIX} directories in {catalogue} - run build_graph_catalogue.py first.")
    if diseases:
        requested = [disease.strip() for disease in diseases.split(",") if disease.strip()]
        missing = sorted(set(requested) - set(available))
        if missing:
            logger.warning("Not in the catalogue (no significant Reactome pathway, or not built): %s", ", ".join(missing))
        available = [disease for disease in available if disease in requested]
    if n_diseases is not None:
        if n_diseases < 1:
            raise SystemExit(f"--n-diseases must be at least 1 (got {n_diseases}).")
        available = available[:n_diseases]
    return available


def score_candidates(seeds, weight, gene_index, graph, restart_prob, bfwalk_alpha, bfwalk_threads):
    """Stage 5's BFWalk and RWR scores of every graph gene that is not a seed."""
    scores = pd.DataFrame({"bfwalk_score": pd.Series(sc.bfwalk_scores(weight, gene_index, seeds, bfwalk_alpha,
                                                                         bfwalk_threads))})
    scores["rwr_score"] = pd.Series(sc.rwr_scores(graph, seeds, restart_prob))
    scores.index.name = "gene"
    return scores.sort_values("bfwalk_score", ascending=False, kind="stable")


# Each worker process loads the disease's graph once (set by init_worker) instead of receiving it per seed.
_WORKER = {}


def init_worker(disease_dir, restart_prob, bfwalk_alpha, bfwalk_threads):
    # Ctrl+C reaches every process of the terminal's process group; only the parent handles it (see
    # run_seeds), so workers do not each print a traceback.
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    weight, _, gene_index, _ = load_graph(disease_dir)
    _WORKER.update(weight=weight, gene_index=gene_index, graph=sc.build_digraph(weight, gene_index),
                   restart_prob=restart_prob, bfwalk_alpha=bfwalk_alpha, bfwalk_threads=bfwalk_threads)


def hold_out(gene, seeds, seed_dir):
    """Score the candidates with gene removed from the seeds; write seed_dir and return gene's record."""
    started = time.time()
    remaining = sorted(set(seeds) - {gene})
    seed_dir.mkdir(parents=True, exist_ok=True)
    record = {"gene": gene, "n_seeds": len(remaining)}
    if not remaining:
        record.update(status="no_seed_left")
    else:
        scores = score_candidates(remaining, _WORKER["weight"], _WORKER["gene_index"], _WORKER["graph"],
                                  _WORKER["restart_prob"], _WORKER["bfwalk_alpha"], _WORKER["bfwalk_threads"])
        scores.to_csv(seed_dir / "candidate_scores.tsv", sep="\t")
        record.update(status="ranked", n_candidates=len(scores))
        for method in METHODS:
            record[f"{method}_score"] = float(scores.loc[gene, f"{method}_score"])
            record[f"{method}_rank"] = mid_rank(scores[f"{method}_score"], gene)
    (seed_dir / "result.json").write_text(json.dumps(record, indent=2))
    lines = [f"held-out gene: {gene}", f"remaining seeds: {len(remaining)}", f"status: {record['status']}"]
    if record["status"] == "ranked":
        for method in METHODS:
            lines.append(f"{method}: score {record[f'{method}_score']:.6g}, rank {record[f'{method}_rank']:g} "
                         f"of {record['n_candidates']} candidates")
    lines.append(f"restart_prob {_WORKER['restart_prob']}, bfwalk_alpha {_WORKER['bfwalk_alpha']}, "
                 f"{time.time() - started:.1f} s")
    (seed_dir / "seed.log").write_text("\n".join(lines) + "\n")
    return record


def disease_logger(log_path):
    """A logger writing to the disease's own leave_one_out.log as well as the run's log and the console."""
    child = logging.getLogger(f"leave_one_out.{log_path.parent.name}")
    for handler in list(child.handlers):
        child.removeHandler(handler)
    handler = logging.FileHandler(log_path, mode="w")
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s", "%Y-%m-%d %H:%M:%S"))
    child.addHandler(handler)
    return child


def run_disease(disease, args, ks, rng):
    """Leave-one-out for one disease. Returns (summary rows, details of the ranked seeds) or None."""
    disease_dir = args.catalogue / f"{disease}{GRAPH_SUFFIX}"
    out_dir = args.out_dir / f"{disease}{LOO_SUFFIX}"
    out_dir.mkdir(parents=True, exist_ok=True)
    log = disease_logger(out_dir / "leave_one_out.log")

    _, _, gene_index, _ = load_graph(disease_dir)
    seed_table = SeedGenesSchema.validate(pd.read_csv(disease_dir / "disease_seed_genes.tsv", sep="\t"))
    evidence = dict(zip(seed_table["gene_id"], seed_table["score"]))
    seeds = sorted(set(seed_table["gene_id"]) & set(gene_index))
    held_out = seeds
    if args.max_seeds is not None and len(seeds) > args.max_seeds:
        held_out = sorted(rng.choice(seeds, args.max_seeds, replace=False))
    log.info("%s: %d genes in the graph, %d seeds (%d of the disease's %d scored genes are not in the graph and "
             "cannot be held out); holding out %d.", disease, len(gene_index), len(seeds),
             len(seed_table) - len(seeds), len(seed_table), len(held_out))
    if len(seeds) < 2:
        log.warning("Fewer than 2 seeds in the graph: holding one out leaves nothing to propagate from - skipped.")
        return None

    bfwalk_threads = max(1, (os.cpu_count() or 1) // max(1, args.jobs))
    init_args = (disease_dir, args.restart_prob, args.bfwalk_alpha, bfwalk_threads)
    started = time.time()
    records = run_seeds(seeds, held_out, out_dir, args.jobs, init_args)
    details = pd.DataFrame(records)
    details["genetic_evidence"] = details["gene"].map(evidence)
    details = details.reindex(columns=RECORD_COLUMNS)
    details.to_csv(out_dir / "scores_left_out_details.tsv", sep="\t", index=False, na_rep="NA")

    ranked = details[details["status"] == "ranked"]
    rankings = {method: ranked[f"{method}_rank"] for method in METHODS}
    table = metrics_table(rankings, ranked["n_candidates"], ks, args.n_random, rng)
    table.to_csv(out_dir / "metrics.tsv", sep="\t", index=False)
    plot_roc(rankings, ranked["n_candidates"], out_dir / "roc_curve.png",
             f"{disease}: leave-one-out, {len(ranked)} held-out seeds", rng)
    log.info("Leave-one-out metrics (%d held-out seeds, median %d candidates, %.0f s):\n%s", len(ranked),
             int(ranked["n_candidates"].median()), time.time() - started, format_metrics(table))

    summary = []
    for method in METHODS:
        row = {"diseaseId": disease, "method": method, "n_held_out": len(ranked),
               "median_candidates": float(ranked["n_candidates"].median())}
        for metric_row in table[table["ranking"] == method].itertuples():
            row[metric_row.metric] = metric_row.value
            row[f"{metric_row.metric}_p_vs_random"] = metric_row.p_vs_random
        summary.append(row)
    return summary, ranked.assign(diseaseId=disease)


def run_seeds(seeds, held_out, out_dir, jobs, init_args):
    """The full run and one run per held-out seed, in worker processes. On Ctrl+C the queued runs are cancelled
    and the workers are terminated at once, instead of the pool waiting for every queued seed to finish."""
    pool = ProcessPoolExecutor(max_workers=max(1, jobs), initializer=init_worker, initargs=init_args)
    try:
        full = pool.submit(score_candidates_full, seeds, out_dir / "full")
        futures = [pool.submit(hold_out, gene, seeds, out_dir / gene) for gene in held_out]
        records = [future.result() for future in futures]
        full.result()
    except BaseException:
        # pool._processes is private, but it is the only handle on the workers; terminating them is what
        # stops a BFWalk run that is already in progress.
        workers = list((pool._processes or {}).values())
        pool.shutdown(wait=False, cancel_futures=True)
        for worker in workers:
            worker.terminate()
        raise
    pool.shutdown()
    return records


def score_candidates_full(seeds, full_dir):
    """All seeds, nothing held out: the run compare_scores.py compares with Open Targets."""
    full_dir.mkdir(parents=True, exist_ok=True)
    score_candidates(seeds, _WORKER["weight"], _WORKER["gene_index"], _WORKER["graph"], _WORKER["restart_prob"],
                     _WORKER["bfwalk_alpha"], _WORKER["bfwalk_threads"]).to_csv(full_dir / "candidate_scores.tsv",
                                                                               sep="\t")


def setup_logging(log_path):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s", "%Y-%m-%d %H:%M:%S")
    for handler in [logging.FileHandler(log_path, mode="w"), logging.StreamHandler()]:
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    logger.setLevel(logging.INFO)


def main(argv=None):
    args = parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(args.out_dir / "leave_one_out.log")
    ks = parse_ks(args.k)
    rng = np.random.default_rng(args.seed)
    logger.info("leave_one_out.py %s", " ".join(sys.argv[1:] if argv is None else argv))

    diseases = select_diseases(args.catalogue, args.diseases, args.n_diseases)
    logger.info("Diseases: %d.", len(diseases))
    summaries, all_ranked = [], []
    try:
        for n, disease in enumerate(diseases, start=1):
            logger.info("[%d/%d] %s", n, len(diseases), disease)
            result = run_disease(disease, args, ks, rng)
            if result is not None:
                summaries.extend(result[0])
                all_ranked.append(result[1])
    except KeyboardInterrupt:
        logger.error("Interrupted (Ctrl+C): all jobs stopped. Diseases finished before it keep their results; "
                     "no summary was written.")
        return 130

    if not summaries:
        logger.error("No disease could be validated - no metrics.")
        return 1
    summary = pd.DataFrame(summaries)
    summary.to_csv(args.out_dir / "leave_one_out_summary.tsv", sep="\t", index=False)
    pooled = pd.concat(all_ranked, ignore_index=True)
    plot_roc({method: pooled[f"{method}_rank"] for method in METHODS}, pooled["n_candidates"],
             args.out_dir / "roc_curve_all_diseases.png",
             f"Leave-one-out, {pooled['diseaseId'].nunique()} diseases, {len(pooled)} held-out seeds", rng)
    per_method = summary.groupby("method")[["auroc"] + [c for c in summary.columns if c.startswith("recall@")
                                                         and not c.endswith("p_vs_random")]].mean()
    logger.info("Mean over %d diseases (each disease weighted equally):\n%s", summary["diseaseId"].nunique(),
                per_method.to_string(float_format=lambda v: f"{v:.3f}"))
    logger.info("Wrote %s and %s.", args.out_dir / "leave_one_out_summary.tsv",
                args.out_dir / "roc_curve_all_diseases.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
