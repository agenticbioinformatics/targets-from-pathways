"""Stage 5 - candidate scoring by network propagation from the disease's seed genes.

Seeds are the disease's genetic-association genes (Stage 2's disease_seed_genes.tsv) that are nodes of the
Stage 3 or Stage 4 graph. Every other graph gene is a candidate, scored by its proximity to the seeds, and
written to candidate_scores.tsv. Methods (--method):
- topology (always computed): (1 + number of seed neighbours) / (1 + weighted shortest-path distance from
  the nearest seed, each edge costing 1 / weight; a zero-weight edge costs as much as the weakest
  positive-weight edge). Unreachable genes score 0.
- rwr (opt-in): random walk with restart, restarting uniformly on the seeds (personalised PageRank).
- bfwalk (opt-in): BFWalk backtrack-free network propagation from the seeds
  (https://github.com/jedrzejkubica/BFWalk). Needs the bfwalk package (for now from TestPyPI, see README).
The output is ranked by the last requested method in the order topology, rwr, bfwalk.
"""

import argparse
import logging
from pathlib import Path

import networkx as nx
import pandas as pd

from stage0_schemas import ENSEMBL_GENE_ID_PATTERN, CandidateScoresSchema, SeedGenesSchema, fail, load_graph

logger = logging.getLogger("score_candidates")

METHODS = ["topology", "rwr", "bfwalk"]


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="stage5_score_candidates.py",
        description="Stage 5 - candidate scoring. Propagates from the disease's seed genes over the graph and "
        "ranks every non-seed graph gene: topology (always) and, opt-in via --method, random-walk-with-restart "
        "and/or BFWalk. Reads the graph files in --graph-dir (Stage 3 or Stage 4) and the seed genes "
        "(--seed-genes); writes candidate_scores.tsv.",
    )
    parser.add_argument("--graph-dir", required=True, type=Path, help="Directory with the Stage 3/4 graph files.")
    parser.add_argument("--seed-genes", type=Path, default=None,
                        help="Seed genes: a TSV with a gene_id column (Stage 2's disease_seed_genes.tsv). "
                        "Default: <graph-dir>/disease_seed_genes.tsv.")
    parser.add_argument("--method", default="topology",
                        help="Comma-separated subset of {topology, rwr, bfwalk}. Default: topology (always computed).")
    parser.add_argument("--restart-prob", type=float, default=0.5,
                        help="RWR restart probability (0 < p < 1); PageRank alpha = 1 - p. Default: 0.5.")
    parser.add_argument("--bfwalk-alpha", type=float, default=0.5,
                        help="BFWalk attenuation coefficient (0 < alpha < 1). Default: 0.5.")
    parser.add_argument("--bfwalk-threads", type=int, default=0,
                        help="BFWalk threads; 0 = all available cores. Default: 0.")
    parser.add_argument("--out", type=Path, default=None, help="Output TSV. Default: <graph-dir>/candidate_scores.tsv.")
    return parser.parse_args(argv)


def parse_methods(raw_methods):
    methods = ["topology"]
    for method in raw_methods.split(","):
        method = method.strip().lower()
        if not method:
            continue
        if method not in METHODS:
            fail(f"--method has unknown value {method!r}; valid: {METHODS}.")
        if method not in methods:
            methods.append(method)
    return methods


def read_seeds(seeds_path, gene_index):
    """Return the seed genes that are graph nodes, sorted. Seeds outside the graph (in no size-capped
    pathway) cannot propagate anything, so they are dropped with a log line."""
    if not seeds_path.exists():
        fail(f"--seed-genes {seeds_path} does not exist (run stage2_gsea_discovery.py first?).")
    all_seeds = set(SeedGenesSchema.validate(pd.read_csv(seeds_path, sep="\t"))["gene_id"])
    seeds = sorted(all_seeds & set(gene_index))
    logger.info("Seeds: %d/%d seed genes are graph nodes (%s).", len(seeds), len(all_seeds), seeds_path)
    if not seeds:
        fail(f"None of the {len(all_seeds)} seed genes in {seeds_path} is a node of the graph.")
    return seeds


def build_digraph(weight, gene_index):
    """Each edge gets "weight" (used by RWR, 0 stays 0 so evidence-free edges carry no flow) and
    "cost" (1 / weight, used for shortest paths). Stage 4 stores weight 0 on edges between genes without
    evidence. Such an edge still exists, so it is kept, but it costs as much as the weakest positive-weight
    edge: it must never be cheaper than an edge that has evidence."""
    graph = nx.DiGraph()
    graph.add_nodes_from(gene_index)
    edges = weight.tocoo()
    positive = edges.data > 0
    zero_weight_cost = 1.0
    if positive.any():
        zero_weight_cost = 1.0 / float(edges.data[positive].min())
    for i, j, edge_weight in zip(edges.row, edges.col, edges.data):
        edge_weight = max(float(edge_weight), 0.0)
        if edge_weight > 0:
            cost = 1.0 / edge_weight
        else:
            cost = zero_weight_cost
        graph.add_edge(gene_index[i], gene_index[j], weight=edge_weight, cost=cost)
    return graph


def topology_scores(graph, seeds):
    distance_from_seeds = nx.multi_source_dijkstra_path_length(graph, set(seeds), weight="cost")
    # Co-membership links all genes of a pathway, so seed neighbours approximate pathways shared with seeds.
    undirected = graph.to_undirected(as_view=True)
    seed_set = set(seeds)

    scores = {}
    for gene in graph.nodes:
        if gene in seed_set:
            continue
        if gene not in distance_from_seeds:
            scores[gene] = 0.0
            continue
        n_seed_neighbours = len(seed_set & set(undirected.neighbors(gene)))
        scores[gene] = (1.0 + n_seed_neighbours) / (1.0 + distance_from_seeds[gene])
    return scores


def rwr_scores(graph, seeds, restart_prob):
    if not 0.0 < restart_prob < 1.0:
        fail(f"--restart-prob must be strictly between 0 and 1 (got {restart_prob}).")
    restart = {}
    for gene in seeds:
        restart[gene] = 1.0 / len(seeds)
    # networkx PageRank is a deterministic power iteration, so no seed is needed.
    pagerank = nx.pagerank(graph, alpha=1.0 - restart_prob, personalization=restart, weight="weight")
    seed_set = set(seeds)
    scores = {}
    for gene in graph.nodes:
        if gene not in seed_set:
            scores[gene] = pagerank[gene]
    return scores


def bfwalk_scores(weight, gene_index, seeds, alpha, threads):
    """BFWalk from the seeds, on the directed weighted graph. BFWalk needs edge weights in ]0, 1], so
    zero-weight edges are dropped (as in RWR, they carry no flow) and the rest are divided by the largest."""
    if not 0.0 < alpha < 1.0:
        fail(f"--bfwalk-alpha must be strictly between 0 and 1 (got {alpha}).")
    try:
        from BFWalk.bfwalk import calculate_scores
    except ImportError:
        fail("--method bfwalk needs the bfwalk package: pip install --index-url https://test.pypi.org/simple/ "
             "--extra-index-url https://pypi.org/simple/ bfwalk")

    edges = weight.tocoo()
    positive = edges.data > 0
    max_weight = edges.data[positive].max() if positive.any() else 1.0
    network = []
    for i, j, edge_weight in zip(edges.row[positive], edges.col[positive], edges.data[positive]):
        network.append((int(i), int(j), float(edge_weight) / max_weight))

    node2idx = {}
    for i, gene in enumerate(gene_index):
        node2idx[gene] = i
    seeds_vector = [0.0] * len(gene_index)
    for gene in seeds:
        seeds_vector[node2idx[gene]] = 1.0
    logger.info("BFWalk: %d edges with weight > 0 (scaled by max weight %.3g), alpha=%.2f.",
                len(network), max_weight, alpha)
    raw_scores = calculate_scores(network, node2idx, seeds_vector, alpha=alpha, threads=threads)

    seed_set = set(seeds)
    scores = {}
    for gene, score in zip(gene_index, raw_scores):
        if gene not in seed_set:
            scores[gene] = float(score)
    return scores


def run(args):
    methods = parse_methods(args.method)
    weight, _, gene_index, _ = load_graph(args.graph_dir)
    for gene in gene_index:
        if not ENSEMBL_GENE_ID_PATTERN.match(gene):
            fail(f"Graph gene index has a non-Ensembl ID: {gene}")

    seeds = read_seeds(args.seed_genes or args.graph_dir / "disease_seed_genes.tsv", gene_index)
    seed_set = set(seeds)
    candidates = []
    for gene in gene_index:
        if gene not in seed_set:
            candidates.append(gene)
    logger.info("Scoring %d candidates from %d seeds (methods: %s).", len(candidates), len(seeds), ",".join(methods))
    if not candidates:
        logger.warning("Every graph gene is a seed - writing an empty ranking.")

    graph = build_digraph(weight, gene_index)
    result = pd.DataFrame({"gene": pd.Series(candidates, dtype="object")})
    result["topology_score"] = result["gene"].map(topology_scores(graph, seeds)).astype(float)
    rank_column = "topology_score"
    if "rwr" in methods:
        result["rwr_score"] = result["gene"].map(rwr_scores(graph, seeds, args.restart_prob)).astype(float)
        rank_column = "rwr_score"
    if "bfwalk" in methods:
        scores = bfwalk_scores(weight, gene_index, seeds, args.bfwalk_alpha, args.bfwalk_threads)
        result["bfwalk_score"] = result["gene"].map(scores).astype(float)
        rank_column = "bfwalk_score"

    result = result.sort_values(rank_column, ascending=False, kind="stable").reset_index(drop=True)
    result = CandidateScoresSchema.validate(result)
    out_path = args.out or args.graph_dir / "candidate_scores.tsv"
    result.to_csv(out_path, sep="\t", index=False)
    logger.info("Stage 5 complete (%d candidates ranked by %s): %s", len(result), rank_column, out_path)
    return out_path


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    run(parse_args(argv))


if __name__ == "__main__":
    main()
