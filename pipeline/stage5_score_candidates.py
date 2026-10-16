"""Stage 5 - candidate scoring.

Scores every non-target gene of the Stage 3 or Stage 4 graph by its proximity to the target and writes
candidate_scores.tsv. Methods (--method):
- topology (always computed): (1 + common neighbours with the target) / (1 + weighted shortest-path distance
  from the target, each edge costing 1 / weight). Unreachable genes score 0.
- rwr (opt-in): random walk with restart on the target (personalised PageRank).
"""

import argparse
import logging
from pathlib import Path

import networkx as nx
import pandas as pd

from stage0_schemas import ENSEMBL_GENE_ID_PATTERN, CandidateScoresSchema, fail, load_graph

logger = logging.getLogger("score_candidates")

METHODS = ["topology", "rwr"]
# Stage 4 stores weight 0 on edges between genes without evidence. Such an edge still exists, so it costs
# one plain hop instead of disconnecting the graph.
ZERO_WEIGHT_HOP_COST = 1.0


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="stage5_score_candidates.py",
        description="Stage 5 - candidate scoring. Ranks every non-target graph gene by proximity to the "
        "target: topology (always) and, opt-in via --method, random-walk-with-restart. Reads the graph files "
        "in --graph-dir (Stage 3 or Stage 4); writes candidate_scores.tsv.",
    )
    parser.add_argument("--graph-dir", required=True, type=Path, help="Directory with the Stage 3/4 graph files.")
    parser.add_argument("--method", default="topology",
                        help="Comma-separated subset of {topology, rwr}. Default: topology (always computed).")
    parser.add_argument("--target", default=None,
                        help="Target Ensembl gene ID. Default: resolved_target.gene_id from graph_metadata.json.")
    parser.add_argument("--restart-prob", type=float, default=0.5,
                        help="RWR restart probability (0 < p < 1); PageRank alpha = 1 - p. Default: 0.5.")
    parser.add_argument(
        "--node-weight-mix", type=float, default=0.5,
        help="RWR on a Stage 4 graph only: extra restart mass (in units of the target's) spread over genes by "
        "genetic-evidence score. 0 = restart on the target alone. Default: 0.5.",
    )
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


def build_digraph(weight, gene_index):
    """Each edge gets "weight" (used by RWR, 0 stays 0 so evidence-free edges carry no flow) and
    "cost" (1 / weight, used for shortest paths)."""
    graph = nx.DiGraph()
    graph.add_nodes_from(gene_index)
    edges = weight.tocoo()
    for i, j, edge_weight in zip(edges.row, edges.col, edges.data):
        edge_weight = max(float(edge_weight), 0.0)
        if edge_weight > 0:
            cost = 1.0 / edge_weight
        else:
            cost = ZERO_WEIGHT_HOP_COST
        graph.add_edge(gene_index[i], gene_index[j], weight=edge_weight, cost=cost)
    return graph


def topology_scores(graph, target):
    distance_from_target = nx.single_source_dijkstra_path_length(graph, target, weight="cost")
    # Co-membership links all genes of a pathway, so shared neighbours approximate shared pathways.
    undirected = graph.to_undirected(as_view=True)
    target_neighbours = set(undirected.neighbors(target))

    scores = {}
    for gene in graph.nodes:
        if gene == target:
            continue
        if gene not in distance_from_target:
            scores[gene] = 0.0
            continue
        n_common = len(target_neighbours & set(undirected.neighbors(gene)))
        scores[gene] = (1.0 + n_common) / (1.0 + distance_from_target[gene])
    return scores


def restart_distribution(metadata, gene_index, target, node_weight_mix):
    """The target always gets 1 unit of restart mass. On a Stage 4 graph, node_weight_mix more units are
    spread over genes in proportion to their evidence score, so the target keeps at least 1 / (1 + mix)."""
    restart = {}
    for gene in gene_index:
        restart[gene] = 0.0
    restart[target] = 1.0

    evidence_scores = metadata.get("genetic_evidence_score")
    if not evidence_scores:
        return restart
    total_evidence = 0.0
    for gene in gene_index:
        total_evidence += max(0.0, float(evidence_scores.get(gene, 0.0)))
    if total_evidence <= 0:
        return restart

    mix = max(node_weight_mix, 0.0)
    for gene in gene_index:
        restart[gene] += mix * max(0.0, float(evidence_scores.get(gene, 0.0))) / total_evidence
    total_mass = sum(restart.values())
    for gene in gene_index:
        restart[gene] = restart[gene] / total_mass
    logger.info("RWR restart: target keeps %.2f of the mass (--node-weight-mix=%.2f).", restart[target], mix)
    return restart


def rwr_scores(graph, gene_index, target, metadata, restart_prob, node_weight_mix):
    if not 0.0 < restart_prob < 1.0:
        fail(f"--restart-prob must be strictly between 0 and 1 (got {restart_prob}).")
    restart = restart_distribution(metadata, gene_index, target, node_weight_mix)
    # networkx PageRank is a deterministic power iteration, so no seed is needed.
    pagerank = nx.pagerank(graph, alpha=1.0 - restart_prob, personalization=restart, weight="weight")
    scores = {}
    for gene in graph.nodes:
        if gene != target:
            scores[gene] = pagerank[gene]
    return scores


def run(args):
    methods = parse_methods(args.method)
    weight, _, gene_index, metadata = load_graph(args.graph_dir)
    for gene in gene_index:
        if not ENSEMBL_GENE_ID_PATTERN.match(gene):
            fail(f"Graph gene index has a non-Ensembl ID: {gene}")

    target = args.target or metadata.get("resolved_target", {}).get("gene_id")
    if not target:
        fail("No --target given and graph_metadata.json has no resolved_target.gene_id.")
    if target not in gene_index:
        fail(f"Target {target} is not a node in this graph ({len(gene_index)} nodes).")
    logger.info("Scoring %d candidates against %s (methods: %s).", len(gene_index) - 1, target, ",".join(methods))

    graph = build_digraph(weight, gene_index)
    candidates = []
    for gene in gene_index:
        if gene != target:
            candidates.append(gene)
    if not candidates:
        logger.warning("Graph has no nodes other than the target - writing an empty ranking.")

    result = pd.DataFrame({"gene": pd.Series(candidates, dtype="object")})
    result["topology_score"] = result["gene"].map(topology_scores(graph, target)).astype(float)
    rank_column = "topology_score"
    if "rwr" in methods:
        scores = rwr_scores(graph, gene_index, target, metadata, args.restart_prob, args.node_weight_mix)
        result["rwr_score"] = result["gene"].map(scores).astype(float)
        rank_column = "rwr_score"

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
