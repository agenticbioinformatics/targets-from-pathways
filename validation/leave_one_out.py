"""Leave-one-out validation over the disease graph catalogue (validation/build_graph_catalogue.py).

For each disease, its graph and seeds come from <catalogue>/<diseaseId>_disease_graph/. Every seed that is a
graph node is held out in turn: it is removed from the seeds, every graph gene that is not one of the remaining
seeds is a candidate (the held-out gene among them), and the candidates are scored with Stage 5's BFWalk and
random walk with restart (pipeline/stage5_score_candidates.py), all remaining seeds weighted equally. The held-out gene's rank among the candidates
measures whether the pipeline recovers it. The graph itself is not rebuilt: its pathways were selected by a GSEA
that included the held-out gene, so this measures recovery given the disease's pathways.

Writes into --out-dir/<diseaseId>_leave_one_out/:
- full/candidate_scores.tsv       all seeds, nothing held out (read by compare_scores.py)
- <gene>/candidate_scores.tsv     every candidate's bfwalk_score and rwr_score with <gene> held out
- <gene>/result.json, seed.log    the held-out gene's scores and ranks, and what was run
- scores_left_out_details.tsv     one row per held-out seed (read by compare_scores.py)
- cache_key.json                  what the scores were computed from (see below)
- metrics.tsv, leave_one_out.log  AUROC and recall@k of BFWalk and RWR next to a random ranking's expectation
- roc_curve.png                   ROC curves of BFWalk and RWR, with the random ranking (AUROC 0.500) as diagonal
and into --out-dir: leave_one_out_summary.tsv (one row per disease and method), leave_one_out.log,
roc_curve_all_diseases.png (every held-out seed of every disease pooled), roc_curves_per_disease.png (one curve
per disease and method, each method in shades of its own colour) and roc_curves_per_disease_<method>.png (one
curve per disease for that method); the per-disease plots also show the pooled curve of all diseases.

Caching: a disease is scored again only if something its scores depend on changed: the sha256 of its graph and
seed files (CACHE_INPUT_FILES), --restart-prob, --bfwalk-alpha, the BFWalk version or the held-out seeds (set by
--max-seeds and --seed). Otherwise its scores are read back from scores_left_out_details.tsv; the metrics and
plots are always recomputed, so --k takes effect. cache_key.json is written last, so an interrupted disease is
scored again. --recompute ignores the cache.

Metrics, over the held-out seeds:
- rank: mid-rank among the candidates (ties count half), 1 = best.
- AUROC: per gene, the fraction of the other candidates it outscores (ties count half); averaged over genes.
  The ROC curve plots, for each fraction x of the candidates taken from the top, the fraction of held-out genes
  ranked within it; its area is this AUROC.
- recall@k: the fraction of held-out genes ranked within the top k.
- random: the expectation of a uniformly random ranking of the same candidates (AUROC 0.500; recall@k = k /
  number of candidates, averaged); its ROC curve is the diagonal. No random ranking is sampled.
"""

import argparse
import importlib.metadata
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
from stage0_schemas import SeedGenesSchema, load_graph, sha256_of_file  # noqa: E402

METHODS = ["bfwalk", "rwr"]
GRAPH_SUFFIX = "_disease_graph"
LOO_SUFFIX = "_leave_one_out"
CACHE_FILE = "cache_key.json"
# The files of a catalogue disease that its scores are computed from (graph_metadata.json is only a summary).
CACHE_INPUT_FILES = ["graph_weight.npz", "graph_sign.npz", "graph_gene_index.json", "disease_seed_genes.tsv"]
RECORD_COLUMNS = ["gene", "status", "genetic_evidence", "n_seeds", "n_candidates",
                  "bfwalk_score", "bfwalk_rank", "rwr_score", "rwr_rank"]
# Per-disease ROC plots: each method's diseases are shades of its colour map.
METHOD_COLOUR_MAPS = {"bfwalk": "Blues", "rwr": "Oranges"}
MAX_LEGEND_DISEASES = 20
# A uniformly random ranking's expected AUROC; its ROC curve is the diagonal.
RANDOM_AUROC = 0.5

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
    parser.add_argument("--seed", type=int, default=0, help="Random seed for the --max-seeds sample. Default: 0.")
    parser.add_argument("--restart-prob", type=float, default=0.5, help="RWR restart probability. Default: 0.5.")
    parser.add_argument("--bfwalk-alpha", type=float, default=0.5, help="BFWalk attenuation. Default: 0.5.")
    parser.add_argument("--jobs", type=int, default=1, help="Held-out seeds scored in parallel. Default: 1.")
    parser.add_argument("--recompute", action="store_true",
                        help="Score every disease again, even if its cached results match its inputs.")
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


def random_expected(n_candidates, ks):
    """Expected metrics of a uniformly random ranking: AUROC 0.5, and recall@k = k / candidates, averaged."""
    n_candidates = np.asarray(n_candidates, dtype=int)
    expected = {"auroc": RANDOM_AUROC}
    for k in ks:
        expected[f"recall@{k}"] = float(np.mean(np.minimum(k, n_candidates) / n_candidates))
    return expected


def metrics_table(rankings, n_candidates, ks):
    """rankings: {name: per-gene ranks}. Returns a DataFrame of metric rows, with the random ranking's
    expectation next to each."""
    expected = random_expected(n_candidates, ks)
    rows = []
    for name, ranks in rankings.items():
        for metric, value in ranking_metrics(ranks, n_candidates, ks).items():
            rows.append({"ranking": name, "metric": metric, "value": value, "random_expected": expected[metric]})
    return pd.DataFrame(rows)


def format_metrics(table):
    wide = table.pivot(index="metric", columns="ranking", values="value")
    order = [m for m in ["auroc"] + [c for c in table["metric"].unique() if c.startswith("recall@")]]
    wide = wide.loc[order]
    wide["random (expected)"] = table.drop_duplicates("metric").set_index("metric").loc[order, "random_expected"]
    return wide.to_string(float_format=lambda v: f"{v:.3f}")


def parse_ks(raw):
    ks = []
    for value in raw.split(","):
        if value.strip():
            ks.append(int(value))
    return ks


def roc_curve(ranks, n_candidates):
    """(x, y, AUROC): x = fraction of the other candidates ranked above the held-out gene, y = the fraction of
    held-out genes recovered within that top fraction. Area under the curve = the mean AUROC."""
    n_candidates = np.asarray(n_candidates, dtype=float)
    false_positive_fraction = (np.asarray(ranks, dtype=float) - 1) / np.maximum(n_candidates - 1, 1)
    x = np.concatenate([[0.0], np.sort(false_positive_fraction), [1.0]])
    y = np.concatenate([[0.0], np.arange(1, len(false_positive_fraction) + 1) / len(false_positive_fraction), [1.0]])
    return x, y, 1.0 - float(np.mean(false_positive_fraction))


def style_roc_axis(axis, title):
    """Axes, labels, and the random ranking: the diagonal, AUROC 0.500."""
    axis.plot([0, 1], [0, 1], linestyle=":", color="grey", linewidth=1.2, label=f"random (AUROC {RANDOM_AUROC:.3f})")
    axis.set_xlabel("fraction of candidates ranked above the held-out gene")
    axis.set_ylabel("fraction of held-out genes recovered")
    axis.set_title(title, fontsize=10)
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)


def plot_roc(rankings, n_candidates, path, title):
    """ROC curve of each ranking (see roc_curve), next to the random ranking's diagonal."""
    figure, axis = plt.subplots(figsize=(5, 5))
    for name, ranks in rankings.items():
        x, y, auroc = roc_curve(ranks, n_candidates)
        axis.step(x, y, where="post", label=f"{name} (AUROC {auroc:.3f})")
    style_roc_axis(axis, title)
    axis.legend(loc="lower right", fontsize=9)
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    plt.close(figure)


def plot_roc_per_disease(pooled, methods, path, title):
    """One ROC curve per disease and method, plus a thick curve of all diseases pooled per method.

    pooled: the ranked held-out seeds of every disease (diseaseId, n_candidates, <method>_rank). With several
    methods, each method has its own colour and its diseases are shades of it; a disease has the same shade in
    every method, so its curves can be matched. With one method, the diseases have distinct colours and the
    pooled curve is black. Diseases are listed in the legend only up to MAX_LEGEND_DISEASES."""
    diseases = sorted(pooled["diseaseId"].unique())
    list_diseases = len(diseases) <= MAX_LEGEND_DISEASES
    figure, axis = plt.subplots(figsize=(8.5, 5.5))
    for method in methods:
        if len(methods) == 1:
            prefix = ""
            colours = [plt.cm.tab20(i % 20) for i in range(len(diseases))]
            pooled_colour = "black"
        else:
            prefix = f"{method}: "
            colour_map = plt.get_cmap(METHOD_COLOUR_MAPS[method])
            colours = colour_map(np.linspace(0.3, 0.75, len(diseases)) if len(diseases) > 1 else [0.55])
            pooled_colour = colour_map(1.0)
        x, y, auroc = roc_curve(pooled[f"{method}_rank"], pooled["n_candidates"])
        axis.step(x, y, where="post", color=pooled_colour, linewidth=2.8, zorder=3,
                  label=f"{prefix}all diseases (AUROC {auroc:.3f})")
        for i, disease in enumerate(diseases):
            rows = pooled[pooled["diseaseId"] == disease]
            x, y, auroc = roc_curve(rows[f"{method}_rank"], rows["n_candidates"])
            if list_diseases:
                label = f"{prefix}{disease} (AUROC {auroc:.3f})"
            elif i == 0:
                label = f"{prefix}{len(diseases)} diseases, one curve each"
            else:
                label = None
            axis.step(x, y, where="post", color=colours[i], linewidth=1, alpha=0.85 if list_diseases else 0.4,
                      label=label)
    style_roc_axis(axis, title)
    axis.legend(loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=7.5, frameon=False)
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


def package_version(name):
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


def cache_key(disease_dir, held_out, args):
    """Everything a disease's leave-one-out scores depend on."""
    input_sha256 = {}
    for name in CACHE_INPUT_FILES:
        input_sha256[name] = sha256_of_file(disease_dir / name)
    return {"input_sha256": input_sha256, "restart_prob": args.restart_prob, "bfwalk_alpha": args.bfwalk_alpha,
            "bfwalk_version": package_version("bfwalk"), "held_out": [str(gene) for gene in held_out]}


def cached_details(out_dir, key):
    """The disease's scores_left_out_details.tsv if out_dir holds a complete run for exactly this key, else None."""
    try:
        stored = json.loads((out_dir / CACHE_FILE).read_text())
    except (FileNotFoundError, json.JSONDecodeError):
        return None
    details_path = out_dir / "scores_left_out_details.tsv"
    if stored != key or not details_path.exists() or not (out_dir / "full" / "candidate_scores.tsv").exists():
        return None
    details = pd.read_csv(details_path, sep="\t")
    for gene in details.loc[details["status"] == "ranked", "gene"]:
        if not (out_dir / gene / "candidate_scores.tsv").exists():
            return None
    return details


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

    started = time.time()
    key = cache_key(disease_dir, held_out, args)
    details = None if args.recompute else cached_details(out_dir, key)
    if details is not None:
        log.info("Cached: %s already holds the scores for exactly these input files and parameters - scoring "
                 "skipped (--recompute to score again).", out_dir)
    else:
        # Removed first, so a run interrupted from here on is not taken for a complete one.
        (out_dir / CACHE_FILE).unlink(missing_ok=True)
        bfwalk_threads = max(1, (os.cpu_count() or 1) // max(1, args.jobs))
        init_args = (disease_dir, args.restart_prob, args.bfwalk_alpha, bfwalk_threads)
        records = run_seeds(seeds, held_out, out_dir, args.jobs, init_args)
        details = pd.DataFrame(records)
        details["genetic_evidence"] = details["gene"].map(evidence)
        details = details.reindex(columns=RECORD_COLUMNS)
        details.to_csv(out_dir / "scores_left_out_details.tsv", sep="\t", index=False, na_rep="NA")
        (out_dir / CACHE_FILE).write_text(json.dumps(key, indent=2))

    ranked = details[details["status"] == "ranked"]
    rankings = {method: ranked[f"{method}_rank"] for method in METHODS}
    table = metrics_table(rankings, ranked["n_candidates"], ks)
    table.to_csv(out_dir / "metrics.tsv", sep="\t", index=False)
    plot_roc(rankings, ranked["n_candidates"], out_dir / "roc_curve.png",
             f"{disease}: leave-one-out, {len(ranked)} held-out seeds")
    log.info("Leave-one-out metrics (%d held-out seeds, median %d candidates, %.0f s):\n%s", len(ranked),
             int(ranked["n_candidates"].median()), time.time() - started, format_metrics(table))

    summary = []
    for method in METHODS:
        row = {"diseaseId": disease, "method": method, "n_held_out": len(ranked),
               "median_candidates": float(ranked["n_candidates"].median())}
        for metric_row in table[table["ranking"] == method].itertuples():
            row[metric_row.metric] = metric_row.value
            row[f"{metric_row.metric}_random_expected"] = metric_row.random_expected
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
             f"Leave-one-out, {pooled['diseaseId'].nunique()} diseases, {len(pooled)} held-out seeds")
    per_disease_title = f"Leave-one-out per disease, {pooled['diseaseId'].nunique()} diseases"
    plot_roc_per_disease(pooled, METHODS, args.out_dir / "roc_curves_per_disease.png", per_disease_title)
    for method in METHODS:
        plot_roc_per_disease(pooled, [method], args.out_dir / f"roc_curves_per_disease_{method}.png",
                             f"{per_disease_title}: {method}")
    per_method = summary.groupby("method")[["auroc"] + [c for c in summary.columns if c.startswith("recall@")
                                                         and not c.endswith("random_expected")]].mean()
    logger.info("Mean over %d diseases (each disease weighted equally):\n%s", summary["diseaseId"].nunique(),
                per_method.to_string(float_format=lambda v: f"{v:.3f}"))
    logger.info("Wrote %s, %s and roc_curves_per_disease*.png.", args.out_dir / "leave_one_out_summary.tsv",
                args.out_dir / "roc_curve_all_diseases.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
