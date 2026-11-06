"""Leave-one-out validation over the disease graph catalogue (validation/build_graph_catalogue.py).

For each disease (<catalogue>/<diseaseId>_disease_graph/), every seed (in the graph)is held out iteratively:
BFWalk and RWR (pipeline/stage5_score_candidates.py) score other genes in the graph.
Five rank aggregation methods ask whether BFWalk and RWR are
complementary and which aggregation method is best. Each held-out gene's rank under each ranking gives AUROC, AUPRC,
precision@k and recall@k. Combined_score_tests check each aggregation against BFWalk, RWR
and random.

A degree-only ranking (each candidate's number of graph neighbours, no seeds used) is the baseline any
propagation has to beat: on these graphs seeds are hubs, so recovering them partly means ranking hubs high.

--refit-gsea removes the circularity of holding a seed out of the propagation only: the seed also helped GSEA
choose the pathways the graph is built from. Each held-out seed's GSEA is rerun without it, the pathways are
selected again and the graph rebuilt from them (gsea_refit.py); a seed left outside the refitted graph is ranked
with score 0, as unreachable.

--holdout known_drug evaluates on a second, independent evidence line instead: the disease's Open Targets
known_drug targets that are graph genes but not seeds, each ranked by the full-seed scores against the
candidates that are not known_drug targets. No seed is removed, so nothing leaks into GSEA or the graph.

Caching (seed hold-outs): a disease is scored again only if something its scores depend on changed: the sha256 of
its graph and seed files (CACHE_INPUT_FILES), --restart-prob, --bfwalk-alpha, the BFWalk version, the held-out
seeds (set by --max-seeds and --seed) or the refit settings. Otherwise its scores are read back from
scores_left_out_details.tsv using cache_key.json (--recompute ignores the cache and recomputes).
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
import scipy.sparse as sp
from scipy.special import digamma
from scipy.stats import beta, wilcoxon

VALIDATION_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(VALIDATION_DIR.parent / "pipeline"))
import stage3_build_graph as bg
import stage5_score_candidates as sc
from stage0_schemas import SeedGenesSchema, load_graph, sha256_of_file

sys.path.insert(0, str(VALIDATION_DIR))

METHODS = ["bfwalk", "rwr"]
# Rank aggregations of BFWalk and RWR (see combined_scores): whether the two methods are complementary.
COMBINATIONS = ["rrf", "borda", "rra", "union", "intersection"]
# Rankings that ignore the seeds: what a propagation method has to beat.
BASELINES = ["degree"]
RANKINGS = METHODS + COMBINATIONS + BASELINES
RANKING_LABELS = {"bfwalk": "BFWalk", "rwr": "RWR", "rrf": "RRF", "borda": "Borda", "rra": "RRA",
                  "union": "union top-k", "intersection": "intersection top-k", "degree": "degree only",
                  "open_targets": "Open Targets"}
# One colour per ranking in every plot.
RANKING_COLOURS = {"bfwalk": "tab:blue", "rwr": "tab:orange", "rrf": "tab:green", "borda": "tab:red",
                   "rra": "tab:purple", "union": "tab:brown", "intersection": "tab:pink", "degree": "tab:gray",
                   "open_targets": "black"}
# Reciprocal rank fusion's constant (Cormack et al. 2009): damps the weight of the very top ranks.
RRF_K = 60
GRAPH_SUFFIX = "_disease_graph"
LOO_SUFFIX = "_leave_one_out"
CACHE_FILE = "cache_key.json"
# The files of a catalogue disease that its scores are computed from (graph_metadata.json is only a summary).
CACHE_INPUT_FILES = ["graph_weight.npz", "graph_sign.npz", "graph_gene_index.json", "disease_seed_genes.tsv"]
RECORD_COLUMNS = (["gene", "status", "in_graph", "genetic_evidence", "n_seeds", "n_candidates",
                   "bfwalk_score", "bfwalk_rank", "rwr_score", "rwr_rank", "degree_score", "degree_rank"]
                  + [f"{combination}_rank" for combination in COMBINATIONS]
                  + ["n_pathways", "n_pathways_lost", "n_pathways_gained", "pathway_jaccard"])
# Bumped when the records change, so caches written by an earlier version are recomputed.
RECORD_VERSION = 2
HOLDOUTS = ["seeds", "known_drug"]
DEGREE_QUARTILES = ["lowest 25%", "25-50%", "50-75%", "top 25%"]
KNOWN_DRUG_SUFFIX = "_known_drug"
# BFWalk's peak memory grows with the square of the graph's nodes and not with its threads (measured: 1.4 GB at 2,411
# nodes, 5.4 GB at 5,001, the same with 1 and 2 threads). Every --jobs worker holds its own copy, so parallel workers
# can exhaust a laptop's memory and freeze it: at most this fraction of the machine's memory is planned for BFWalk.
BFWALK_BYTES_PER_NODE_SQUARED = 230
MEMORY_FRACTION = 0.6
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
    parser.add_argument("--min-pooled-genes", type=int, default=200,
                        help="Leave diseases whose graph has fewer genes out of the all-diseases outputs: their numbers "
                        "are too unstable to compare. Their own directories are still written. Default: 200.")
    parser.add_argument("--min-pooled-seeds", type=int, default=5,
                        help="As --min-pooled-genes, for the seeds in the graph. Default: 5.")
    parser.add_argument("--holdout", choices=HOLDOUTS, default="seeds",
                        help="seeds (default): hold out each seed in turn. known_drug: rank the disease's known_drug "
                        "targets that are graph genes but not seeds, by the full-seed scores (needs --associations).")
    parser.add_argument("--associations", type=Path, default=None,
                        help="--holdout known_drug: Open Targets association_by_datatype_indirect directory, e.g. "
                        "data/association_by_datatype_indirect.")
    parser.add_argument("--refit-gsea", action="store_true",
                        help="Rerun each held-out seed's GSEA without it and rebuild the graph from the pathways "
                        "selected then (see gsea_refit.py). Needs --gsea-catalogue-dir, --target-dir and --fi-file.")
    parser.add_argument("--gsea-catalogue-dir", type=Path, default=None,
                        help="--refit-gsea: the DisPathCatalogue v2 directory with merged_gene_sets.gmt and "
                        "genetic_gene_scores_input.parquet.")
    parser.add_argument("--target-dir", type=Path, default=None,
                        help="--refit-gsea: Open Targets target parquet directory, e.g. data/opentargets/26.06/target.")
    parser.add_argument("--fi-file", type=Path, default=None,
                        help="--refit-gsea: Reactome functional-interaction file, unzipped.")
    parser.add_argument("--fdr", type=float, default=0.05,
                        help="--refit-gsea: FDR cutoff, as build_graph_catalogue.py --fdr. Default: 0.05.")
    args = parser.parse_args(argv)
    if args.holdout == "known_drug" and args.associations is None:
        parser.error("--holdout known_drug needs --associations.")
    if args.holdout == "known_drug" and args.refit_gsea:
        parser.error("--refit-gsea applies to seed hold-outs only: known_drug targets are not in the GSEA ranking.")
    if args.refit_gsea and None in (args.gsea_catalogue_dir, args.target_dir, args.fi_file):
        parser.error("--refit-gsea needs --gsea-catalogue-dir, --target-dir and --fi-file.")
    return args


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
              label=f"after top {follow_up_k} genes\n(better ranking is higher and further right)")
    axis.set_xlabel("recall: fraction of held-out genes within the top k")
    axis.set_ylabel("precision / prevalence (x-fold over random)" if enrichment
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
    draw_pr(pr_axis, rankings, n_candidates, "PRC" + (" (x-fold enrichment over random)" if enrichment else ""),
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
    BFWalk, RWR, degree only and a random ranking (one sample drawn with rng), and one-sided Wilcoxon signed-rank p-values,
    paired by gene, that the combination's is greater. Genes scored equally by both are dropped
    (zero_method="wilcox"); NA when all are. ranked: the ranked held-out seeds (n_candidates, <ranking>_rank)."""
    n_candidates = ranked["n_candidates"]
    other_ranks = {method: ranked[f"{method}_rank"] for method in METHODS + BASELINES}
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


def node_degree(weight, gene_index):
    """Each gene's number of distinct graph neighbours, either direction, self-loops not counted."""
    adjacency = ((weight != 0) + (weight.T != 0)).tolil()
    adjacency.setdiag(0)
    return pd.Series(np.diff(adjacency.tocsr().indptr), index=list(gene_index), dtype=float, name="degree_score")


def score_candidates(seeds, weight, gene_index, graph, restart_prob, bfwalk_alpha, bfwalk_threads, degree=None):
    """Stage 5's BFWalk and RWR scores of every graph gene that is not a seed, and its degree (the baseline)."""
    scores = pd.DataFrame({"bfwalk_score": pd.Series(sc.bfwalk_scores(weight, gene_index, seeds, bfwalk_alpha,
                                                                         bfwalk_threads))})
    scores["rwr_score"] = pd.Series(sc.rwr_scores(graph, seeds, restart_prob))
    degree = node_degree(weight, gene_index) if degree is None else degree
    scores["degree_score"] = degree.reindex(scores.index)
    scores.index.name = "gene"
    return scores.sort_values("bfwalk_score", ascending=False, kind="stable")


def add_unreachable(scores, gene):
    """scores with gene added at score 0 for every ranking: a gene outside the graph cannot be reached."""
    if gene in scores.index:
        return scores
    unreachable = pd.DataFrame(0.0, index=pd.Index([gene], name=scores.index.name), columns=scores.columns)
    return pd.concat([scores, unreachable])


# Each worker process loads the disease's graph once (set by init_worker) instead of receiving it per seed.
_WORKER = {}


def init_worker(disease_dir, restart_prob, bfwalk_alpha, bfwalk_threads, refit=None):
    """refit: None, or the disease's refit context (see refit_context); disease_dir is then the graph of the
    full refitted selection."""
    # Ctrl+C reaches every process of the terminal's process group; only the parent handles it (see
    # run_seeds), so workers do not each print a traceback.
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    weight, _, gene_index, _ = load_graph(disease_dir)
    _WORKER.update(weight=weight, gene_index=gene_index, graph=sc.build_digraph(weight, gene_index),
                   degree=node_degree(weight, gene_index), restart_prob=restart_prob, bfwalk_alpha=bfwalk_alpha,
                   bfwalk_threads=bfwalk_threads, refit=refit)


def refitted_graph(gene):
    """(weight, gene_index, digraph, degree, selection summary) for gene's held-out run with --refit-gsea:
    the full graph when its GSEA selects the same pathways without it, else the graph of the new selection
    (weight None when no pathway is left)."""
    import build_graph_catalogue as gc
    import gsea_refit

    refit = _WORKER["refit"]
    selection = refit["gsea"].selection_without(refit["symbols"].get(gene, set()))
    summary = gsea_refit.selection_summary(refit["gsea"].full_selection, selection)
    if selection == refit["gsea"].full_selection:
        return _WORKER["weight"], _WORKER["gene_index"], _WORKER["graph"], _WORKER["degree"], summary
    if not selection:
        return None, [], None, None, summary
    graph = gc.build_disease_graph(selection, refit["gene_sets"], refit["interactions_by_gene_a"])
    weight, _, gene_index = bg.graph_to_sparse(graph)
    return weight, gene_index, sc.build_digraph(weight, gene_index), node_degree(weight, gene_index), summary


def hold_out(gene, seeds, seed_dir):
    """Score the candidates with gene removed from the seeds (those in the graph used for gene); write seed_dir
    and return gene's record."""
    started = time.time()
    weight, gene_index, graph, degree = _WORKER["weight"], _WORKER["gene_index"], _WORKER["graph"], _WORKER["degree"]
    record = {"gene": gene}
    if _WORKER["refit"] is not None:
        weight, gene_index, graph, degree, summary = refitted_graph(gene)
        record.update(summary)
    remaining = sorted((set(seeds) - {gene}) & set(gene_index))
    seed_dir.mkdir(parents=True, exist_ok=True)
    record.update(n_seeds=len(remaining), in_graph=gene in set(gene_index))
    if weight is None:
        record.update(status="no_pathway_left")
    elif not remaining:
        record.update(status="no_seed_left")
    else:
        scores = score_candidates(remaining, weight, gene_index, graph, _WORKER["restart_prob"],
                                  _WORKER["bfwalk_alpha"], _WORKER["bfwalk_threads"], degree)
        scores = add_unreachable(scores, gene)
        scores.to_csv(seed_dir / "candidate_scores.tsv", sep="\t")
        record.update(status="ranked", **rank_record(scores, gene))
    (seed_dir / "result.json").write_text(json.dumps(record, indent=2))
    lines = [f"held-out gene: {gene}", f"remaining seeds: {len(remaining)}", f"status: {record['status']}",
             f"in the graph: {record['in_graph']}"]
    if "n_pathways" in record:
        lines.append(f"refitted GSEA: {record['n_pathways']} pathways ({record['n_pathways_lost']} lost, "
                     f"{record['n_pathways_gained']} gained, Jaccard {record['pathway_jaccard']:.3f})")
    if record["status"] == "ranked":
        for method in METHODS + BASELINES:
            lines.append(f"{method}: score {record[f'{method}_score']:.6g}, rank {record[f'{method}_rank']:g} "
                         f"of {record['n_candidates']} candidates")
        for combination in COMBINATIONS:
            lines.append(f"{combination} ({RANKING_LABELS[combination]}): rank {record[f'{combination}_rank']:g}")
    lines.append(f"restart_prob {_WORKER['restart_prob']}, bfwalk_alpha {_WORKER['bfwalk_alpha']}, "
                 f"{time.time() - started:.1f} s")
    (seed_dir / "seed.log").write_text("\n".join(lines) + "\n")
    return record


def rank_record(scores, gene):
    """gene's score and rank under every ranking, among scores' candidates."""
    record = {"n_candidates": len(scores)}
    for method in METHODS + BASELINES:
        record[f"{method}_score"] = float(scores.loc[gene, f"{method}_score"])
        record[f"{method}_rank"] = mid_rank(scores[f"{method}_score"], gene)
    for combination, combined in combined_scores(scores).items():
        record[f"{combination}_rank"] = mid_rank(combined, gene)
    return record


def jobs_that_fit(n_genes, jobs, memory_bytes):
    """How many of the requested workers fit in memory_bytes on a graph of n_genes; 0 if not even one does."""
    per_run = BFWALK_BYTES_PER_NODE_SQUARED * n_genes ** 2
    return min(max(1, jobs), int(MEMORY_FRACTION * memory_bytes // per_run))


def fitting_jobs(n_genes, jobs, log):
    """jobs_that_fit on this machine, logging when it lowers --jobs or leaves nothing to run."""
    memory_bytes = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    fitted = jobs_that_fit(n_genes, jobs, memory_bytes)
    needed_gb = BFWALK_BYTES_PER_NODE_SQUARED * n_genes ** 2 / 1e9
    if fitted == 0:
        log.warning("Skipped: BFWalk on %d genes needs about %.1f GB, more than %.0f%% of this machine's %.0f GB.",
                    n_genes, needed_gb, 100 * MEMORY_FRACTION, memory_bytes / 1e9)
    elif fitted < jobs:
        log.warning("--jobs %d lowered to %d: each BFWalk run on %d genes needs about %.1f GB, and this machine has "
                    "%.0f GB.", jobs, fitted, n_genes, needed_gb, memory_bytes / 1e9)
    return fitted


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
    key = {"record_version": RECORD_VERSION, "input_sha256": input_sha256, "restart_prob": args.restart_prob,
           "bfwalk_alpha": args.bfwalk_alpha, "bfwalk_version": package_version("bfwalk"),
           "held_out": [str(gene) for gene in held_out], "refit_gsea": args.refit_gsea}
    if args.refit_gsea:
        key["refit"] = {"fdr": args.fdr, "gmt_sha256": sha256_of_file(args.gsea_catalogue_dir / "merged_gene_sets.gmt"),
                        "scores_sha256": sha256_of_file(args.gsea_catalogue_dir / "genetic_gene_scores_input.parquet"),
                        "fi_sha256": sha256_of_file(args.fi_file)}
    return key


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


def load_refit_inputs(args):
    """What every disease's --refit-gsea shares: the GSEA library, the Ensembl gene sets and interactions the
    graphs are built from, and the symbol index. Loaded once."""
    import build_graph_catalogue as gc
    import gsea_refit
    from stage1_ingest import (CoverageTracker, build_genes_table, build_interactions_table, build_symbol_index,
                               read_target_table)

    target_paths = sorted(args.target_dir.glob("*.parquet"))
    index = build_symbol_index(build_genes_table(read_target_table(target_paths)))
    gene_sets, _ = gc.read_reactome_gene_sets(args.gsea_catalogue_dir / "merged_gene_sets.gmt", index)
    interactions = build_interactions_table(args.fi_file, index, "reactome", CoverageTracker(source="reactome_fi"),
                                            curated_only=True, min_score=None)
    return {"library": gsea_refit.CatalogueGsea(args.gsea_catalogue_dir / "merged_gene_sets.gmt", target_paths),
            "index": index, "gene_sets": gene_sets, "interactions_by_gene_a": interactions.set_index("gene_a")}


def refit_context(disease, disease_dir, out_dir, args, log):
    """The disease's GSEA refitted with every seed (its null sampled once here), and the graph of the pathways
    it selects, written to out_dir/refit_full_graph. Returns (graph directory, context for the workers)."""
    import build_graph_catalogue as gc
    import gsea_refit

    shared = args.refit_inputs
    started = time.time()
    scores = gsea_refit.read_disease_scores(args.gsea_catalogue_dir / "genetic_gene_scores_input.parquet", disease)
    gsea = gsea_refit.DiseaseRefit(shared["library"], scores, args.fdr)
    catalogue = gsea_refit.catalogue_selection(disease_dir, shared["gene_sets"])
    log.info("Refit GSEA with every seed: %d pathways selected (%.0f s); the catalogue graph has %d (Jaccard %.3f). "
             "Held-out genes are compared with this selection's graph.", len(gsea.full_selection),
             time.time() - started, len(catalogue), gsea_refit.jaccard(gsea.full_selection, catalogue))
    graph_dir = out_dir / "refit_full_graph"
    graph = gc.build_disease_graph(gsea.full_selection, shared["gene_sets"], shared["interactions_by_gene_a"])
    graph_dir.mkdir(parents=True, exist_ok=True)
    weight, sign, gene_index = bg.graph_to_sparse(graph)
    sp.save_npz(graph_dir / "graph_weight.npz", weight)
    sp.save_npz(graph_dir / "graph_sign.npz", sign)
    (graph_dir / "graph_gene_index.json").write_text(json.dumps(gene_index, indent=2))
    (graph_dir / "graph_metadata.json").write_text(json.dumps(
        {"n_nodes": graph.number_of_nodes(), "n_edges": graph.number_of_edges(),
         "pathways": sorted(gsea.full_selection)}, indent=2))
    context = {"gsea": gsea, "symbols": gsea_refit.symbols_by_gene(scores, shared["index"]),
               "gene_sets": shared["gene_sets"], "interactions_by_gene_a": shared["interactions_by_gene_a"]}
    return graph_dir, context


def known_drug_targets(associations, disease):
    """{target: known_drug score > 0} of one disease (association_by_datatype_indirect)."""
    import pyarrow.compute as pc
    import pyarrow.dataset as ds

    table = ds.dataset(str(associations), format="parquet").to_table(
        filter=(pc.field("diseaseId") == disease) & (pc.field("datatypeId") == "known_drug") & (pc.field("score") > 0),
        columns=["targetId", "score"]).to_pandas()
    return dict(zip(table["targetId"], table["score"]))


def run_known_drug(disease, args, ks):
    """--holdout known_drug for one disease: the full-seed scores, and each known_drug target that is a graph
    gene but not a seed ranked among the candidates that are not known_drug targets (so every target is ranked
    against the same negatives, one positive at a time, as a held-out seed is). Returns as run_disease."""
    disease_dir = args.catalogue / f"{disease}{GRAPH_SUFFIX}"
    out_dir = args.out_dir / f"{disease}{KNOWN_DRUG_SUFFIX}"
    out_dir.mkdir(parents=True, exist_ok=True)
    log = disease_logger(out_dir / "leave_one_out.log")
    weight, _, gene_index, _ = load_graph(disease_dir)
    seed_table = SeedGenesSchema.validate(pd.read_csv(disease_dir / "disease_seed_genes.tsv", sep="\t"))
    seeds = sorted(set(seed_table["gene_id"]) & set(gene_index))
    drug = known_drug_targets(args.associations, disease)
    targets = sorted(set(drug) & (set(gene_index) - set(seeds)))
    log.info("%s: %d genes in the graph, %d seeds; %d known_drug targets, %d of them graph genes that are not seeds.",
             disease, len(gene_index), len(seeds), len(drug), len(targets))
    if not seeds or not targets:
        log.warning("No seed or no known_drug target to rank - skipped.")
        return None
    if not fitting_jobs(len(gene_index), 1, log):
        return None
    started = time.time()
    scores = score_candidates(seeds, weight, gene_index, sc.build_digraph(weight, gene_index), args.restart_prob,
                              args.bfwalk_alpha, os.cpu_count() or 1)
    scores.to_csv(out_dir / "candidate_scores.tsv", sep="\t")
    negatives = scores.drop(index=targets)
    records = []
    for target in targets:
        records.append({"gene": target, "status": "ranked", "in_graph": True, "genetic_evidence": np.nan,
                        "n_seeds": len(seeds), **rank_record(pd.concat([negatives, scores.loc[[target]]]), target)})
    details = pd.DataFrame(records).reindex(columns=RECORD_COLUMNS)
    details["known_drug_score"] = details["gene"].map(drug)
    details.to_csv(out_dir / "scores_known_drug_details.tsv", sep="\t", index=False, na_rep="NA")
    return summarise_disease(disease, details, out_dir, ks, args, log, started, "known_drug targets",
                             len(gene_index), len(seeds))


def summarise_disease(disease, details, out_dir, ks, args, log, started, what, n_genes, n_seeds):
    """Metrics, curves and the disease's summary rows from its records. n_genes, n_seeds: in the graph."""
    ranked = details[details["status"] == "ranked"]
    rankings = {ranking: ranked[f"{ranking}_rank"] for ranking in RANKINGS}
    table = metrics_table(rankings, ranked["n_candidates"], ks, args.follow_up_k, args.seed)
    table.to_csv(out_dir / "metrics.tsv", sep="\t")
    plot_roc_pr(rankings, ranked["n_candidates"], out_dir / "roc_pr_curve.png",
                f"{disease}: {len(ranked)} held-out {what}, median "
                f"{int(ranked['n_candidates'].median())} candidates", args.follow_up_k)
    not_in_graph = int((~ranked["in_graph"].astype(bool)).sum())
    log.info("Metrics (%d held-out %s, %d outside the graph, median %d candidates, %.0f s):\n%s", len(ranked), what,
             not_in_graph, int(ranked["n_candidates"].median()), time.time() - started, format_metrics(table))

    summary = table.reset_index()
    summary.insert(0, "diseaseId", disease)
    summary.insert(2, "n_held_out", len(ranked))
    summary.insert(3, "n_candidates", float(ranked["n_candidates"].median()))
    summary.insert(4, "n_not_in_graph", not_in_graph)
    summary.insert(5, "n_genes", n_genes)
    summary.insert(6, "n_seeds_in_graph", n_seeds)
    return summary, ranked.assign(diseaseId=disease)


def run_disease(disease, args, ks, rng):
    """Leave-one-out for one disease. Returns (summary rows, details of the ranked seeds) or None."""
    if args.holdout == "known_drug":
        return run_known_drug(disease, args, ks)
    disease_dir = args.catalogue / f"{disease}{GRAPH_SUFFIX}"
    out_dir = args.out_dir / f"{disease}{LOO_SUFFIX}"
    out_dir.mkdir(parents=True, exist_ok=True)
    log = disease_logger(out_dir / "leave_one_out.log")

    graph_dir, refit = disease_dir, None
    if args.refit_gsea:
        graph_dir, refit = refit_context(disease, disease_dir, out_dir, args, log)
    _, _, gene_index, _ = load_graph(graph_dir)
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
        jobs = fitting_jobs(len(gene_index), args.jobs, log)
        if not jobs:
            return None
        # Removed first, so a run interrupted from here on is not taken for a complete one.
        (out_dir / CACHE_FILE).unlink(missing_ok=True)
        bfwalk_threads = max(1, (os.cpu_count() or 1) // jobs)
        init_args = (graph_dir, args.restart_prob, args.bfwalk_alpha, bfwalk_threads, refit)
        records = run_seeds(seeds, held_out, out_dir, jobs, init_args, sorted(set(seed_table["gene_id"])))
        details = pd.DataFrame(records)
        details["genetic_evidence"] = details["gene"].map(evidence)
        details = details.reindex(columns=RECORD_COLUMNS)
        details.to_csv(out_dir / "scores_left_out_details.tsv", sep="\t", index=False, na_rep="NA")
        (out_dir / CACHE_FILE).write_text(json.dumps(key, indent=2))

    return summarise_disease(disease, details, out_dir, ks, args, log, started,
                             "seeds (GSEA refitted)" if args.refit_gsea else "seeds", len(gene_index), len(seeds))


def run_seeds(seeds, held_out, out_dir, jobs, init_args, all_seeds):
    """The full run on seeds (those in the graph) and one run per held-out seed, in worker processes. all_seeds:
    every seed of the disease, because a refitted graph can hold seeds the full graph does not. On Ctrl+C the queued runs are cancelled
    and the workers are terminated at once, instead of the pool waiting for every queued seed to finish."""
    pool = ProcessPoolExecutor(max_workers=max(1, jobs), initializer=init_worker, initargs=init_args)
    try:
        full = pool.submit(score_candidates_full, seeds, out_dir / "full")
        futures = [pool.submit(hold_out, gene, all_seeds, out_dir / gene) for gene in held_out]
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
    args.refit_inputs = load_refit_inputs(args) if args.refit_gsea else None
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
    summary = mark_pooled(pd.concat(summaries, ignore_index=True), args.min_pooled_genes, args.min_pooled_seeds)
    summary.to_csv(args.out_dir / "leave_one_out_summary.tsv", sep="\t", index=False)
    pooled_diseases = summary.loc[summary["in_pooled"], "diseaseId"].unique()
    pooled = pd.concat(all_ranked, ignore_index=True)
    pooled = pooled[pooled["diseaseId"].isin(pooled_diseases)]
    quartiles = degree_quartile_table(pooled)
    quartiles.to_csv(args.out_dir / "auroc_by_degree_quartile.tsv", sep="\t")
    logger.info("AUROC by the held-out gene's degree quartile:\n%s", format_metrics(quartiles))
    what = {"seeds": "Leave-one-out", "known_drug": "known_drug targets"}[args.holdout]
    title = f"{what} ({pooled['diseaseId'].nunique()} diseases, {len(pooled)} held-out genes)"
    plot_roc({method: pooled[f"{method}_rank"] for method in METHODS + BASELINES}, pooled["n_candidates"],
             args.out_dir / "roc_curve_all_diseases.png", title)
    plot_roc_pr({ranking: pooled[f"{ranking}_rank"] for ranking in RANKINGS}, pooled["n_candidates"],
                args.out_dir / "roc_pr_curve_combined_all_diseases.png", title, args.follow_up_k, enrichment=True)
    tests = write_combined_score_tests(pooled, args.out_dir / "combined_score_tests.tsv", args.seed)
    # AUPRC and precision@k depend on each disease's prevalence: across diseases only their fold enrichments.
    metrics = (["auroc", "auprc_fold_enrichment", f"precision@{args.follow_up_k}_fold_enrichment"]
               + [f"recall@{k}" for k in ks])
    log_summary(summary[summary["in_pooled"]], tests, metrics)
    logger.info("Wrote %s, auroc_by_degree_quartile.tsv, roc_curve_all_diseases.png, "
                "roc_pr_curve_combined_all_diseases.png and combined_score_tests.tsv.",
                args.out_dir / "leave_one_out_summary.tsv")
    return 0


def mark_pooled(summary, min_genes, min_seeds):
    """summary with in_pooled: whether a disease's graph is large enough for the all-diseases outputs. If none is,
    all are pooled, so that a run on small graphs still gets its outputs."""
    keep = (summary["n_genes"] >= min_genes) & (summary["n_seeds_in_graph"] >= min_seeds)
    left_out = sorted(summary.loc[~keep, "diseaseId"].unique())
    if not keep.any():
        logger.warning("No graph has at least %d genes and %d seeds: pooling all %d diseases, though their numbers "
                       "are unstable.", min_genes, min_seeds, len(left_out))
        return summary.assign(in_pooled=True)
    if left_out:
        logger.info("Left out of the all-diseases outputs, with fewer than %d genes or %d seeds in the graph: %s.",
                    min_genes, min_seeds, ", ".join(left_out))
    return summary.assign(in_pooled=keep)


def degree_quartile_table(pooled):
    """Mean AUROC per ranking by how connected the held-out gene is: quartiles of its degree percentile among its
    disease's candidates (auroc_from_rank of its degree rank). A method that only follows hubs tracks the degree
    baseline from quartile to quartile. Quartiles pool genes, so a disease counts by its number of held-out genes."""
    percentile = pd.Series([auroc_from_rank(rank, n) for rank, n in zip(pooled["degree_rank"], pooled["n_candidates"])],
                           index=pooled.index)
    table = pd.DataFrame({"quartile": pd.cut(percentile, [-np.inf, 0.25, 0.5, 0.75, np.inf], labels=DEGREE_QUARTILES)})
    for ranking in METHODS + BASELINES:
        table[f"auroc_{ranking}"] = [auroc_from_rank(rank, n) for rank, n in
                                     zip(pooled[f"{ranking}_rank"], pooled["n_candidates"])]
    grouped = table.groupby("quartile", observed=False)
    result = grouped.mean()
    result.insert(0, "n_held_out", grouped.size())
    return result


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
    for other in ["random"] + METHODS + BASELINES:
        p_values = pooled.pivot(index="combination", columns="measure", values=f"p_vs_{other}")
        sections.append(f"rank combination method vs {other} (p-value, one-sided Wilcoxon signed-rank test)\n"
                        + format_metrics(p_values.reindex(index=COMBINATIONS, columns=list(PER_GENE))))
    logger.info("Is any rank aggregation better than BFWalk, RWR, or random?\n\n%s", "\n\n".join(sections))

    per_metric, per_disease = method_wins(summary, metrics)
    logger.info("BFWalk vs RWR: the number of diseases (of %d) where each scores higher, per metric:\n%s",
                n_diseases, per_metric.to_string())


if __name__ == "__main__":
    sys.exit(main())
