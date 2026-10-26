"""Leave-one-out validation over the disease graph catalogue (validation/build_graph_catalogue.py).

For each disease (<catalogue>/<diseaseId>_disease_graph/), every seed (in the graph)is held out iteratively:
BFWalk and RWR (pipeline/stage5_score_candidates.py) score other genes in the graph.
Five rank aggregation methods ask whether BFWalk and RWR are
complementary and which aggregation method is best. Each held-out gene's rank under each ranking gives AUROC, AUPRC,
precision@k and recall@k. Combined_score_tests check each aggregation against BFWalk, RWR
and random.

Caching: a disease is scored again only if something its scores depend on changed: the sha256 of its graph and
seed files (CACHE_INPUT_FILES), --restart-prob, --bfwalk-alpha, the BFWalk version or the held-out seeds (set by
--max-seeds and --seed). Otherwise its scores are read back from scores_left_out_details.tsv using cache_key.json
(--recompute ignores the cache and recomputes).
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
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import digamma
from scipy.stats import beta, wilcoxon

VALIDATION_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(VALIDATION_DIR.parent / "pipeline"))
import stage5_score_candidates as sc
from stage0_schemas import SeedGenesSchema, load_graph, sha256_of_file

METHODS = ["bfwalk", "rwr"]
# Rank aggregations of BFWalk and RWR (see combined_scores): whether the two methods are complementary.
COMBINATIONS = ["rrf", "borda", "rra", "union", "intersection"]
RANKINGS = METHODS + COMBINATIONS
RANKING_LABELS = {"bfwalk": "BFWalk", "rwr": "RWR", "rrf": "RRF", "borda": "Borda", "rra": "RRA",
                  "union": "union top-k", "intersection": "intersection top-k", "open_targets": "Open Targets"}
# One colour per ranking in every plot.
RANKING_COLOURS = {"bfwalk": "tab:blue", "rwr": "tab:orange", "rrf": "tab:green", "borda": "tab:red",
                   "rra": "tab:purple", "union": "tab:brown", "intersection": "tab:pink", "open_targets": "black"}
# Reciprocal rank fusion's constant (Cormack et al. 2009): damps the weight of the very top ranks.
RRF_K = 60
GRAPH_SUFFIX = "_disease_graph"
LOO_SUFFIX = "_leave_one_out"
CACHE_FILE = "cache_key.json"
# The files of a catalogue disease that its scores are computed from (graph_metadata.json is only a summary).
CACHE_INPUT_FILES = ["graph_weight.npz", "graph_sign.npz", "graph_gene_index.json", "disease_seed_genes.tsv"]
RECORD_COLUMNS = ["gene", "status", "genetic_evidence", "n_seeds", "n_candidates",
                  "bfwalk_score", "bfwalk_rank", "rwr_score", "rwr_rank"] + [f"{combination}_rank" for combination in COMBINATIONS]
# A uniformly random ranking's expected AUROC; its ROC curve is the diagonal.
RANDOM_AUROC = 0.5
# Random rankings sampled for a random ranking's expected AUPRC (see random_auprc).
RANDOM_DRAWS = 200

logger = logging.getLogger("leave_one_out")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="leave_one_out.py",
        description="Leave-one-out validation over the disease graph catalogue: hold out each seed of a disease "
        "in turn, score the candidates with BFWalk and RWR, and measure the held-out gene's recovery (AUROC, "
        "AUPRC, recall@k) against a random ranking.",
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
    parser.add_argument("--k", default="10,50,100",
                        help="Comma-separated k values for recall@k and precision@k.")
    parser.add_argument("--follow-up-k", type=int, default=20,
                        help="The depth whose genes would realistically be followed up: precision@k at it is "
                        "marked on the precision-recall curves. Default: 20.")
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


def roc_curve(ranks, n_candidates):
    """(x, y, AUROC): x = fraction of the other candidates ranked above the held-out gene, y = the fraction of
    held-out genes recovered within that top fraction. Area under the curve = the mean AUROC."""
    n_candidates = np.asarray(n_candidates, dtype=float)
    false_positive_fraction = (np.asarray(ranks, dtype=float) - 1) / np.maximum(n_candidates - 1, 1)
    x = np.concatenate([[0.0], np.sort(false_positive_fraction), [1.0]])
    y = np.concatenate([[0.0], np.arange(1, len(false_positive_fraction) + 1) / len(false_positive_fraction), [1.0]])
    return x, y, 1.0 - float(np.mean(false_positive_fraction))


def top_k_size(depths, n_candidates, enrichment=False):
    """At each depth k: the candidates in all the runs' top-k lists, sum of min(k, n). With enrichment, the
    held-out genes a random ranking expects there instead, sum of min(k, n) / n."""
    depths = np.atleast_1d(np.asarray(depths, dtype=float))
    n_sorted = np.sort(np.asarray(n_candidates, dtype=float))
    n_below = np.searchsorted(n_sorted, depths)  # runs shorter than k contribute all their candidates
    if enrichment:
        inverse_tail = np.concatenate([np.cumsum((1.0 / n_sorted)[::-1])[::-1], [0.0]])
        return n_below + depths * inverse_tail[n_below]
    head = np.concatenate([[0.0], np.cumsum(n_sorted)])
    return head[n_below] + depths * (len(n_sorted) - n_below)


def pr_curve(ranks, n_candidates, enrichment=False):
    """(recall, precision, AUPRC) of the held-out genes pooled, every run cut at the same depth k: recall = the
    fraction of held-out genes within the top k, precision = those genes / all candidates in the top-k lists. One
    point per distinct rank; AUPRC = average precision, sum of recall step x precision (no interpolation).
    With enrichment, precision / prevalence (see top_k_size): 1 = random, comparable between diseases; within
    one disease (the same n in every run) it is precision x n, and AUPRC becomes AUPRC / prevalence."""
    ranks = np.sort(np.asarray(ranks, dtype=float))
    depths = np.unique(ranks)
    recovered = np.searchsorted(ranks, depths, side="right")
    recall = recovered / len(ranks)
    precision = recovered / top_k_size(depths, n_candidates, enrichment)
    return recall, precision, float(np.sum(np.diff(recall, prepend=0.0) * precision))


def precision_at_k(ranks, n_candidates, k, enrichment=False):
    """The held-out genes in the top-k lists / the candidates in them (/ prevalence with enrichment)."""
    return float(np.sum(np.asarray(ranks, dtype=float) <= k) / top_k_size(k, n_candidates, enrichment)[0])


def prevalence(n_candidates):
    """The fraction of candidates that are held-out genes: a random ranking's precision at any depth."""
    return float(len(n_candidates) / np.sum(np.asarray(n_candidates, dtype=float)))


def random_ranks(n_candidates, rng):
    """A random ranking's sample: each held-out gene a uniformly random rank among its candidates."""
    n_candidates = np.asarray(n_candidates, dtype=float)
    return np.floor(rng.random(len(n_candidates)) * n_candidates) + 1


def random_auprc(n_candidates, rng):
    """A random ranking's expected AUPRC, the mean over RANDOM_DRAWS samples. It is close to the prevalence only
    with many held-out genes: with few, the first points of the curve are noisy and lift it (about 1.2 x the
    prevalence with 50 held-out genes, 2.5 x with 3)."""
    return float(np.mean([pr_curve(random_ranks(n_candidates, rng), n_candidates)[2] for _ in range(RANDOM_DRAWS)]))


def metric_names(ks, follow_up_k):
    names = ["auroc", "auprc", "auprc_fold_enrichment"]
    if follow_up_k is not None:
        names += [f"precision@{follow_up_k}", f"precision@{follow_up_k}_fold_enrichment"]
    return names + [f"recall@{k}" for k in ks]


def metrics_table(rankings, n_candidates, ks, follow_up_k=None, seed=0):
    """One row per ranking ({name: the held-out genes' ranks}) and a last row, random, with a uniformly random
    ranking's expectation. Columns (see metric_names): AUROC; AUPRC and AUPRC / prevalence (pr_curve);
    precision@follow_up_k, as is and / prevalence; recall@k."""
    n_candidates = np.asarray(n_candidates, dtype=float)
    base = prevalence(n_candidates)
    rows = {}
    for name, ranks in rankings.items():
        ranks = np.asarray(ranks, dtype=float)
        auprc = pr_curve(ranks, n_candidates)[2]
        row = {"auroc": float(np.nanmean([auroc_from_rank(r, n) for r, n in zip(ranks, n_candidates)])),
               "auprc": auprc, "auprc_fold_enrichment": auprc / base}
        if follow_up_k is not None:
            row[f"precision@{follow_up_k}"] = precision_at_k(ranks, n_candidates, follow_up_k)
            row[f"precision@{follow_up_k}_fold_enrichment"] = precision_at_k(ranks, n_candidates, follow_up_k, True)
        for k in ks:
            row[f"recall@{k}"] = float(np.mean(ranks <= k))
        rows[name] = row
    auprc = random_auprc(n_candidates, np.random.default_rng(seed))
    row = {"auroc": RANDOM_AUROC, "auprc": auprc, "auprc_fold_enrichment": auprc / base}
    if follow_up_k is not None:
        row[f"precision@{follow_up_k}"] = float(top_k_size(follow_up_k, n_candidates, True)[0]
                                                / top_k_size(follow_up_k, n_candidates)[0])
        row[f"precision@{follow_up_k}_fold_enrichment"] = 1.0
    for k in ks:
        row[f"recall@{k}"] = float(np.mean(np.minimum(k, n_candidates) / n_candidates))
    rows["random"] = row
    table = pd.DataFrame.from_dict(rows, orient="index")[metric_names(ks, follow_up_k)]
    table.index.name = "ranking"
    return table


def format_metrics(table):
    return table.to_string(float_format=lambda v: f"{v:.3g}")


def parse_ks(raw):
    ks = []
    for value in raw.split(","):
        if value.strip():
            ks.append(int(value))
    return ks


def draw_curve(axis, name, x, y, label):
    """A ranking's curve in its own colour; the two methods thicker than the combinations."""
    (line,) = axis.step(x, y, where="post", color=RANKING_COLOURS.get(name), label=label,
                        linewidth=2.0 if name not in COMBINATIONS else 1.2)
    return line


def draw_roc(axis, rankings, n_candidates, title):
    for name, ranks in rankings.items():
        x, y, auroc = roc_curve(ranks, n_candidates)
        draw_curve(axis, name, x, y, f"{RANKING_LABELS.get(name, name)}  AUROC {auroc:.3f}")
    axis.plot([0, 1], [0, 1], linestyle=":", color="grey", linewidth=1.2, label=f"random  AUROC {RANDOM_AUROC:.3f}")
    axis.set_xlabel("fraction of candidates ranked above the held-out gene")
    axis.set_ylabel("fraction of held-out genes recovered")
    axis.set_title(title, fontsize=10)
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)


def draw_pr(axis, rankings, n_candidates, title, follow_up_k, enrichment=False):
    """Precision-recall curve of each ranking (see pr_curve) with a dot at the follow-up depth, and the random
    ranking's expected precision, flat at the prevalence. With enrichment, precision / prevalence ("x" = fold).
    Precision on a log scale: it spans orders of magnitude, and the random line stays visible."""
    baseline = 1.0 if enrichment else prevalence(n_candidates)
    value = (lambda v: f"x{v:.3g}") if enrichment else (lambda v: f"{v:.3g}")
    top = bottom = baseline
    for name, ranks in rankings.items():
        recall, precision, auprc = pr_curve(ranks, n_candidates, enrichment)
        at_k = precision_at_k(ranks, n_candidates, follow_up_k, enrichment)
        # Precision at recall r holds from the previous point up to r (the average-precision area).
        line = draw_curve(axis, name, np.concatenate([[0.0], recall]), np.concatenate([precision, precision[-1:]]),
                          f"{RANKING_LABELS.get(name, name)}  AUPRC {value(auprc)}, P@{follow_up_k} {value(at_k)}")
        if at_k > 0:  # none of the held-out genes in the top k: no dot on a log scale
            axis.plot(np.mean(np.asarray(ranks, dtype=float) <= follow_up_k), at_k, "o", color=line.get_color(),
                      markersize=5, markeredgecolor="black", markeredgewidth=0.6, zorder=3)
        top, bottom = max(top, precision.max()), min(bottom, precision.min())
    axis.axhline(baseline, linestyle=":", color="grey", linewidth=1.2, label=f"random  precision {value(baseline)}")
    axis.plot([], [], "o", color="white", markeredgecolor="black", markeredgewidth=0.6,
              label=f"dot: top {follow_up_k} (follow-up depth)")
    axis.set_xlabel("recall: fraction of held-out genes within the top k")
    axis.set_ylabel("precision / prevalence (fold over random)" if enrichment
                    else "precision: held-out genes / candidates in the top k")
    axis.set_title(title, fontsize=10)
    axis.set_xlim(0, 1)
    axis.set_yscale("log")
    axis.set_ylim(bottom / 1.5, top * 1.5)


def plot_roc(rankings, n_candidates, path, title):
    """ROC curve of each ranking (see roc_curve), next to the random ranking's diagonal."""
    figure, axis = plt.subplots(figsize=(5.5, 5.5))
    draw_roc(axis, rankings, n_candidates, title)
    axis.legend(loc="lower right", fontsize=8.5)
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    plt.close(figure)


def plot_roc_pr(rankings, n_candidates, path, title, follow_up_k, enrichment=False):
    """Two panels: the ROC curves with their AUROC, and the precision-recall curves with their AUPRC and
    precision@follow_up_k (both / prevalence with enrichment, to pool diseases)."""
    figure, (roc_axis, pr_axis) = plt.subplots(1, 2, figsize=(14, 7))
    draw_roc(roc_axis, rankings, n_candidates, "ROC")
    draw_pr(pr_axis, rankings, n_candidates, "Precision-recall" + (", fold over random" if enrichment else ""),
            follow_up_k, enrichment)
    for axis in (roc_axis, pr_axis):
        axis.legend(loc="upper center", bbox_to_anchor=(0.5, -0.11), ncol=2, fontsize=8, frameon=False)
    figure.suptitle(title, fontsize=11)
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    plt.close(figure)


def average_precision_vs_random(ranks, n_candidates):
    """Each held-out gene's average precision (its own run's AUPRC, 1 / rank) / a random ranking's expected
    average precision (the mean of 1 / rank over ranks 1..n, harmonic(n) / n): 1 = random."""
    n_candidates = np.asarray(n_candidates, dtype=float)
    harmonic = digamma(n_candidates + 1) + np.euler_gamma
    return n_candidates / (np.asarray(ranks, dtype=float) * harmonic)


# The per-gene measures the combinations are tested on.
PER_GENE = {
    "auroc": lambda ranks, n_candidates: np.array([auroc_from_rank(r, n) for r, n in zip(ranks, n_candidates)],
                                                  dtype=float),
    "average_precision_vs_random": average_precision_vs_random,
}


def combined_score_tests(ranked, rng):
    """For each per-gene measure (PER_GENE) and combination: its mean over the held-out genes, the means under
    BFWalk, RWR and a random ranking (one sample drawn with rng), and one-sided Wilcoxon signed-rank p-values,
    paired by gene, that the combination's is greater. Genes scored equally by both are dropped
    (zero_method="wilcox"); NA when all are. ranked: the ranked held-out seeds (n_candidates, <ranking>_rank)."""
    n_candidates = ranked["n_candidates"]
    other_ranks = {method: ranked[f"{method}_rank"] for method in METHODS}
    other_ranks["random"] = random_ranks(n_candidates, rng)
    rows = []
    for measure, per_gene in PER_GENE.items():
        others = {name: per_gene(ranks, n_candidates) for name, ranks in other_ranks.items()}
        for combination in COMBINATIONS:
            combined = per_gene(ranked[f"{combination}_rank"], n_candidates)
            row = {"measure": measure, "combination": combination, "n_held_out": len(combined),
                   "mean": float(np.mean(combined))}
            for name, other in others.items():
                row[f"mean_{name}"] = float(np.mean(other))
            for name, other in others.items():
                p_value = np.nan
                if np.any(combined != other):
                    p_value = float(wilcoxon(combined, other, zero_method="wilcox", alternative="greater").pvalue)
                row[f"p_vs_{name}"] = p_value
            rows.append(row)
    return rows


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


def robust_rank_aggregation(ranks):
    """Kolde et al. 2012's rho score of each candidate (lower is better). Its normalised ranks (rank / number of
    candidates) sorted, r(1) <= ... <= r(m); rho = min over j of P(Beta(j, m - j + 1) <= r(j)), the probability
    that the j-th best of m uniformly random normalised ranks is at least as good as r(j). Kolde's p-value
    multiplies rho by m (Bonferroni), which does not change the order. ranks: candidates x methods."""
    normalised = np.sort(ranks.to_numpy() / len(ranks), axis=1)
    order = np.arange(1, normalised.shape[1] + 1)
    return pd.Series(beta.cdf(normalised, order, normalised.shape[1] - order + 1).min(axis=1), index=ranks.index)


def combined_scores(scores):
    """Each candidate's rank aggregations of BFWalk and RWR, from its two ranks (mid-ranks among the candidates,
    1 = best); higher is better, like the method scores. Returns {combination: Series}:
    - rrf: reciprocal rank fusion, the sum of 1 / (RRF_K + rank);
    - borda: Borda count, the sum of the candidates ranked below it (n - rank);
    - rra: robust rank aggregation, rho negated;
    - union: the better rank, negated: a gene joins the union of the two top-k lists at k = its better rank;
    - intersection: the worse rank, negated: a gene joins the intersection at k = its worse rank."""
    ranks = pd.concat([scores[f"{method}_score"].rank(ascending=False, method="average") for method in METHODS],
                      axis=1)
    return {"rrf": (1.0 / (RRF_K + ranks)).sum(axis=1), "borda": (len(ranks) - ranks).sum(axis=1),
            "rra": -robust_rank_aggregation(ranks), "union": -ranks.min(axis=1),
            "intersection": -ranks.max(axis=1)}


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
        for combination, combined in combined_scores(scores).items():
            record[f"{combination}_rank"] = mid_rank(combined, gene)
    (seed_dir / "result.json").write_text(json.dumps(record, indent=2))
    lines = [f"held-out gene: {gene}", f"remaining seeds: {len(remaining)}", f"status: {record['status']}"]
    if record["status"] == "ranked":
        for method in METHODS:
            lines.append(f"{method}: score {record[f'{method}_score']:.6g}, rank {record[f'{method}_rank']:g} "
                         f"of {record['n_candidates']} candidates")
        for combination in COMBINATIONS:
            lines.append(f"{combination} ({RANKING_LABELS[combination]}): rank {record[f'{combination}_rank']:g}")
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


def add_combined_ranks(details, out_dir):
    """The combinations' ranks, computed from each held-out seed's candidate_scores.tsv; rewrites
    scores_left_out_details.tsv. Run on every cache hit: the combinations are cheap and not part of the cache
    key, so results cached by an earlier definition of a combination are brought up to date."""
    ranks = {combination: [] for combination in COMBINATIONS}
    for gene, status in zip(details["gene"], details["status"]):
        combined = {}
        if status == "ranked":
            # round_trip: the default parser can change a score's last digit, which breaks exact ties between
            # candidates and so changes their ranks.
            combined = combined_scores(pd.read_csv(out_dir / gene / "candidate_scores.tsv", sep="\t",
                                                   index_col="gene", float_precision="round_trip"))
        for combination in COMBINATIONS:
            ranks[combination].append(mid_rank(combined[combination], gene) if combined else np.nan)
    details = details.assign(**{f"{combination}_rank": ranks[combination] for combination in COMBINATIONS})
    details = details.reindex(columns=RECORD_COLUMNS)
    details.to_csv(out_dir / "scores_left_out_details.tsv", sep="\t", index=False, na_rep="NA")
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
    log.info("%s: %d genes in the graph, %d seeds (%d / %d disease-associated genes are not in the graph); "
             "holding out %d.", disease, len(gene_index), len(seeds),
             len(seed_table) - len(seeds), len(seed_table), len(held_out))
    if len(seeds) < 2:
        log.warning("Fewer than 2 seeds in the graph: holding one out leaves nothing to propagate from - skipped.")
        return None

    started = time.time()
    key = cache_key(disease_dir, held_out, args)
    details = None if args.recompute else cached_details(out_dir, key)
    if details is not None:
        log.info("Cached: scores already exist (--recompute to score again).")
        details = add_combined_ranks(details, out_dir)
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
    rankings = {ranking: ranked[f"{ranking}_rank"] for ranking in RANKINGS}
    table = metrics_table(rankings, ranked["n_candidates"], ks, args.follow_up_k, args.seed)
    table.to_csv(out_dir / "metrics.tsv", sep="\t")
    plot_roc_pr(rankings, ranked["n_candidates"], out_dir / "roc_pr_curve.png",
                f"{disease}: leave-one-out, {len(ranked)} held-out seeds, median "
                f"{int(ranked['n_candidates'].median())} candidates", args.follow_up_k)
    log.info("Leave-one-out metrics (%d held-out seeds, median %d candidates, %.0f s):\n%s", len(ranked),
             int(ranked["n_candidates"].median()), time.time() - started, format_metrics(table))

    summary = table.reset_index()
    summary.insert(0, "diseaseId", disease)
    summary.insert(2, "n_held_out", len(ranked))
    summary.insert(3, "n_candidates", float(ranked["n_candidates"].median()))
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
    ks = sorted(set(parse_ks(args.k)) | {args.follow_up_k})
    rng = np.random.default_rng(args.seed)
    logger.info("leave_one_out.py %s", " ".join(sys.argv[1:] if argv is None else argv))

    diseases = select_diseases(args.catalogue, args.diseases, args.n_diseases)
    summaries, all_ranked = [], []
    try:
        for n, disease in enumerate(diseases, start=1):
            logger.info("[%d/%d] %s", n, len(diseases), disease)
            result = run_disease(disease, args, ks, rng)
            if result is not None:
                summaries.append(result[0])
                all_ranked.append(result[1])
    except KeyboardInterrupt:
        logger.error("Interrupted (Ctrl+C): all jobs stopped. Diseases finished before it keep their results; "
                     "no summary was written.")
        return 130

    if not summaries:
        logger.error("No disease could be validated - no metrics.")
        return 1
    summary = pd.concat(summaries, ignore_index=True)
    summary.to_csv(args.out_dir / "leave_one_out_summary.tsv", sep="\t", index=False)
    pooled = pd.concat(all_ranked, ignore_index=True)
    title = f"Leave-one-out ({pooled['diseaseId'].nunique()} diseases, {len(pooled)} held-out genes)"
    plot_roc({method: pooled[f"{method}_rank"] for method in METHODS}, pooled["n_candidates"],
             args.out_dir / "roc_curve_all_diseases.png", title)
    plot_roc_pr({ranking: pooled[f"{ranking}_rank"] for ranking in RANKINGS}, pooled["n_candidates"],
                args.out_dir / "roc_pr_curve_combined_all_diseases.png", title, args.follow_up_k, enrichment=True)
    tests = write_combined_score_tests(pooled, args.out_dir / "combined_score_tests.tsv", args.seed)
    # AUPRC and precision@k depend on each disease's prevalence: across diseases only their fold enrichments.
    metrics = (["auroc", "auprc_fold_enrichment", f"precision@{args.follow_up_k}_fold_enrichment"]
               + [f"recall@{k}" for k in ks])
    log_summary(summary, tests, metrics)
    logger.info("Wrote %s, roc_curve_all_diseases.png, roc_pr_curve_combined_all_diseases.png and "
                "combined_score_tests.tsv.", args.out_dir / "leave_one_out_summary.tsv")
    return 0


def write_combined_score_tests(pooled, path, seed):
    """combined_score_tests over all diseases pooled, then per disease; the random ranking's sample is drawn
    with its own generator, so the --max-seeds sample does not change."""
    rng = np.random.default_rng(seed)
    rows = [{"diseaseId": "all", **row} for row in combined_score_tests(pooled, rng)]
    for disease, ranked in pooled.groupby("diseaseId", sort=True):
        rows.extend({"diseaseId": disease, **row} for row in combined_score_tests(ranked, rng))
    tests = pd.DataFrame(rows)
    tests.to_csv(path, sep="\t", index=False, na_rep="NA")
    return tests


def method_wins(summary, metrics):
    """BFWalk vs RWR, disease by disease. Returns (per metric: the number of diseases where BFWalk scores
    higher, where RWR does, and where they tie; per disease: the number of metrics where each scores higher)."""
    by_disease = summary.set_index(["ranking", "diseaseId"])[metrics]
    bfwalk, rwr = by_disease.loc["bfwalk"], by_disease.loc["rwr"]
    per_metric = pd.DataFrame({"bfwalk_higher": (bfwalk > rwr).sum(), "rwr_higher": (rwr > bfwalk).sum(),
                               "tied": (bfwalk == rwr).sum()})
    per_metric.index.name = "metric"
    per_disease = pd.DataFrame({"bfwalk_higher": (bfwalk > rwr).sum(axis=1), "rwr_higher": (rwr > bfwalk).sum(axis=1)})
    return per_metric, per_disease


def log_summary(summary, tests, metrics):
    """Is any rank aggregation better than BFWalk, RWR, or random: each ranking's metrics averaged over the
    diseases, then the pooled p-values of combined_score_tests; last, BFWalk vs RWR disease by disease."""
    n_diseases = summary["diseaseId"].nunique()
    per_ranking = summary.groupby("ranking", sort=False)[metrics].mean()
    pooled = tests[tests["diseaseId"] == "all"]
    sections = [f"Mean over {n_diseases} diseases:\n{format_metrics(per_ranking)}"]
    for other in ["random"] + METHODS:
        p_values = pooled.pivot(index="combination", columns="measure", values=f"p_vs_{other}")
        sections.append(f"rank combination method vs {other} (p-value, one-sided Wilcoxon signed-rank test)\n"
                        + format_metrics(p_values.reindex(index=COMBINATIONS, columns=list(PER_GENE))))
    logger.info("Is any rank aggregation better than BFWalk, RWR, or random?\n\n%s", "\n\n".join(sections))

    per_metric, per_disease = method_wins(summary, metrics)
    logger.info("BFWalk vs RWR: the number of diseases (of %d) where each scores higher, per metric:\n%s",
                n_diseases, per_metric.to_string())


if __name__ == "__main__":
    sys.exit(main())
