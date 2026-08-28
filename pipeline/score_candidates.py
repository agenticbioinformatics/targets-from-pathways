"""Stage 5 — Candidate scoring.

Loads a Stage 3 *or* Stage 4 pathway-based gene-gene graph (the two share
one serialization format — ``graph_weight.npz``, ``graph_sign.npz``,
``graph_gene_index.json``, ``graph_metadata.json``, read via ``--graph-dir``)
plus the target node, and scores every other node for proximity to the
target. Never reads a raw source file or a Stage 1/2 artifact — everything
it needs is in the graph directory.

``--method`` is a comma-separated list drawn from ``{topology, rwr}``,
**defaulting to ``topology`` only**:

- **topology** (always computed — the output always carries a
  ``topology_score`` column): direction- and weight-aware proximity from
  two graph features combined into one score,
  ``(1 + common_neighbours) / (1 + shortest_path_distance)`` —

  - *shortest-path distance* from the target to the node along edge
    directions, each edge costing ``1 / weight`` (``weight`` from
    ``graph_weight.npz``: comembership + interaction-confidence on a Stage 3
    graph, genetic-evidence-derived on a Stage 4 graph). A present edge with
    a stored weight of ``0`` — which Stage 4 produces for an edge between
    two genes with no genetic evidence — costs ``1.0`` (one plain hop), so
    an unweighted or all-zero graph degrades to hop count rather than
    becoming disconnected. A node the target cannot reach scores ``0.0``.
  - *common-neighbour count* with the target in the undirected projection —
    a graph-derived stand-in for shared-pathway membership, since Stage 3's
    co-membership edges make every pair of genes in a pathway adjacent, so
    two genes in the same pathway(s) share that pathway's other members as
    neighbours.

  Deterministic, cheap, no iteration — hence always-on.

- **rwr** (opt-in): personalised random-walk-with-restart, i.e. personalised
  PageRank (``networkx.pagerank``, *not* hand-rolled), restarting on the
  target node, walking ``graph_weight.npz`` as edge weights. ``alpha`` (the
  PageRank damping) is ``1 - --restart-prob``. The restart distribution
  always puts a full unit of mass on the target; when the graph carries
  node weights (a Stage 4 graph's ``genetic_evidence_score``), a further
  ``--node-weight-mix`` units are spread over genes in proportion to their
  genetic-evidence score before renormalising, so genetic evidence seeds
  the walk without ever displacing the target as the primary restart point.
  On a Stage 3 graph (no node weights) the restart is the target alone. The
  result is a probability distribution over nodes (non-negative, sums to ~1).

  ``networkx.pagerank`` is a deterministic power iteration with no random
  step, so the ``seed`` recorded in ``graph_metadata.json`` has nothing to
  seed here; it is logged for provenance and otherwise unused.

Backtrack-free walk (non-backtracking RWR) was considered as a third method
and removed from scope — see README.md plan update 5.

**Output** — ``candidate_scores.tsv`` (or ``--out``), validated against
``schemas.CandidateScoresSchema``: one row per non-target graph node,
columns ``gene`` (Ensembl gene ID, matching the graph index and the Open
Targets convention, so Stage 6 joins on it directly) and ``topology_score``
always, plus ``rwr_score`` only when ``rwr`` was in ``--method``. Rows are
ranked by ``rwr_score`` when present, else ``topology_score``, descending.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import networkx as nx
import pandas as pd
import scipy.sparse as sp

from stage0_schemas import CandidateScoresSchema, ENSEMBL_GENE_ID_PATTERN

logger = logging.getLogger("score_candidates")

METHODS = ("topology", "rwr")
# A present edge whose stored weight is 0 (Stage 4: an edge between two genes
# with no genetic evidence) still exists structurally; give it a plain
# one-hop traversal cost so the graph never silently disconnects.
ZERO_WEIGHT_HOP_COST = 1.0


# ==========================================================================
# CLI
# ==========================================================================


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="score_candidates.py",
        description="Stage 5: score candidate alternative targets by proximity to the target "
        "node over a Stage 3 or Stage 4 pathway-based gene-gene graph.",
    )
    p.add_argument(
        "--graph-dir",
        required=True,
        type=Path,
        help="Directory with graph_weight.npz / graph_sign.npz / graph_gene_index.json / "
        "graph_metadata.json, as written by stage3_build_graph.py or stage4_genetic_evidence_weights.py.",
    )
    p.add_argument(
        "--method",
        default="topology",
        help="Comma-separated subset of {topology, rwr}. Default: topology. rwr is opt-in "
        "(e.g. --method topology,rwr or --method rwr); topology_score is always emitted regardless.",
    )
    p.add_argument(
        "--target",
        default=None,
        help="Target gene's Ensembl ID. Defaults to graph_metadata.json's resolved_target.gene_id.",
    )
    p.add_argument(
        "--restart-prob",
        type=float,
        default=0.5,
        help="RWR restart probability (0 < p < 1); PageRank alpha = 1 - restart_prob. Default: 0.5.",
    )
    p.add_argument(
        "--node-weight-mix",
        type=float,
        default=0.5,
        help="RWR only, and only when the graph carries node weights (Stage 4): extra restart "
        "mass (>= 0, in units of the target's) spread over genes by genetic-evidence score. "
        "0 = restart on the target alone; 0.5 (default) = target keeps ~2/3 of the restart mass. "
        "Ignored on a Stage 3 graph.",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output TSV path. Defaults to <graph-dir>/candidate_scores.tsv.",
    )
    return p


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def _fail(msg: str) -> None:
    logger.error(msg)
    sys.exit(1)


def parse_methods(raw: str) -> list[str]:
    """Ordered, de-duplicated method list. ``topology`` is always included
    (the output always needs a ``topology_score``), even if only ``rwr`` was
    asked for."""
    requested = [m.strip().lower() for m in raw.split(",") if m.strip()]
    unknown = [m for m in requested if m not in METHODS]
    if unknown:
        _fail(f"--method has unknown value(s) {unknown}; valid: {list(METHODS)}.")
    methods = list(dict.fromkeys(requested)) or ["topology"]
    if "topology" not in methods:
        logger.info("topology not requested explicitly, but topology_score is always emitted — adding it.")
        methods = ["topology", *methods]
    return methods


# ==========================================================================
# Graph loading (Stage 3 / Stage 4 shared format)
# ==========================================================================


def load_graph(graph_dir: Path) -> tuple[sp.csr_matrix, sp.csr_matrix, list[str], dict]:
    weight_path = graph_dir / "graph_weight.npz"
    sign_path = graph_dir / "graph_sign.npz"
    index_path = graph_dir / "graph_gene_index.json"
    metadata_path = graph_dir / "graph_metadata.json"
    for path in (weight_path, sign_path, index_path, metadata_path):
        if not path.exists():
            _fail(f"Graph artifact {path} does not exist (run stage3_build_graph.py first?).")

    weight = sp.load_npz(weight_path).tocsr()
    sign = sp.load_npz(sign_path).tocsr()
    gene_index = json.loads(index_path.read_text())
    metadata = json.loads(metadata_path.read_text())

    n = len(gene_index)
    if weight.shape != (n, n) or sign.shape != (n, n):
        _fail(
            f"Graph artifacts in {graph_dir} are misaligned: weight.shape={weight.shape}, "
            f"sign.shape={sign.shape}, gene index has {n} entries."
        )
    return weight, sign, gene_index, metadata


def resolve_target(args_target: str | None, metadata: dict, gene_index: list[str]) -> str:
    target = args_target or metadata.get("resolved_target", {}).get("gene_id")
    if not target:
        _fail("No --target given and graph_metadata.json has no resolved_target.gene_id.")
    if target not in set(gene_index):
        _fail(f"Target {target} is not a node in this graph ({len(gene_index)} nodes).")
    return target


def build_digraph(weight: sp.csr_matrix, gene_index: list[str]) -> nx.DiGraph:
    """Directed graph over *every* gene in the index (isolated nodes kept),
    each edge carrying:

    - ``weight``: the raw stored value, used as-is by RWR. A stored ``0``
      (Stage 4: no genetic evidence on either endpoint) stays ``0`` — that
      edge simply carries no diffusion flow, which is the intended meaning;
      forcing it to a positive value would let no-evidence edges dominate a
      Stage 4 walk.
    - ``cost``: ``1 / weight`` for shortest paths, with a stored ``0``
      costing one plain hop (``ZERO_WEIGHT_HOP_COST``) so a present-but-
      unweighted edge keeps the graph connected for the topology method.
    """
    graph = nx.DiGraph()
    graph.add_nodes_from(gene_index)
    coo = weight.tocoo()
    for i, j, w in zip(coo.row, coo.col, coo.data):
        w = float(w)
        graph.add_edge(
            gene_index[i], gene_index[j],
            weight=max(w, 0.0),
            cost=(1.0 / w) if w > 0 else ZERO_WEIGHT_HOP_COST,
        )
    return graph


# ==========================================================================
# Scoring methods
# ==========================================================================


def topology_scores(graph: nx.DiGraph, target: str) -> dict[str, float]:
    """``(1 + common_neighbours) / (1 + shortest_path_distance)`` per node;
    0.0 for a node unreachable from the target. See module docstring."""
    dist = nx.single_source_dijkstra_path_length(graph, target, weight="cost")
    undirected = graph.to_undirected(as_view=True)
    target_neighbours = set(undirected.neighbors(target))

    scores: dict[str, float] = {}
    for node in graph.nodes:
        if node == target:
            continue
        if node not in dist:
            scores[node] = 0.0
            continue
        common = len(target_neighbours.intersection(undirected.neighbors(node)))
        scores[node] = (1.0 + common) / (1.0 + dist[node])
    return scores


def _restart_distribution(
    metadata: dict, gene_index: list[str], target: str, node_weight_mix: float
) -> dict[str, float]:
    """The RWR restart distribution: always a full unit of mass on the
    target, plus — when the graph carries node weights (Stage 4) —
    ``--node-weight-mix`` extra units spread over genes in proportion to
    their ``genetic_evidence_score``, then renormalised. The target
    therefore always keeps at least ``1 / (1 + mix)`` of the restart mass,
    so "restart on the target" holds for any ``mix``; higher ``mix`` just
    lets genetic evidence seed the walk more strongly. Returns a
    normalised dict over every gene in ``gene_index``."""
    restart = {g: (1.0 if g == target else 0.0) for g in gene_index}
    node_weights = metadata.get("genetic_evidence_score")
    if not node_weights:
        return restart

    total = float(sum(max(0.0, float(node_weights.get(g, 0.0))) for g in gene_index))
    if total <= 0:
        logger.info("Graph carries node weights but all are 0 — RWR restarts on the target alone.")
        return restart

    mix = max(node_weight_mix, 0.0)
    for g in gene_index:
        restart[g] += mix * max(0.0, float(node_weights.get(g, 0.0))) / total
    norm = sum(restart.values())
    restart = {g: v / norm for g, v in restart.items()}
    logger.info(
        "RWR restart distribution: target keeps %.2f of the mass, %.2f spread over genes by "
        "genetic-evidence score (--node-weight-mix=%.2f).",
        restart[target], 1.0 - restart[target], mix,
    )
    return restart


def rwr_scores(
    graph: nx.DiGraph,
    gene_index: list[str],
    target: str,
    metadata: dict,
    restart_prob: float,
    node_weight_mix: float,
) -> dict[str, float]:
    if not 0.0 < restart_prob < 1.0:
        _fail(f"--restart-prob must be strictly between 0 and 1 (got {restart_prob}).")
    personalization = _restart_distribution(metadata, gene_index, target, node_weight_mix)
    pagerank = nx.pagerank(
        graph, alpha=1.0 - restart_prob, personalization=personalization, weight="weight"
    )
    return {node: pagerank[node] for node in graph.nodes if node != target}


# ==========================================================================
# Orchestration
# ==========================================================================


def run(args: argparse.Namespace) -> Path:
    methods = parse_methods(args.method)
    # graph_sign.npz is loaded only for the shape/alignment check — neither
    # topology nor RWR is sign-aware in this prototype (edge signs stay
    # available in the artifact for a later stage that wants them).
    weight, _sign, gene_index, metadata = load_graph(args.graph_dir)
    bad = [g for g in gene_index if not ENSEMBL_GENE_ID_PATTERN.match(g)]
    if bad:
        _fail(f"Graph gene index has non-Ensembl ID(s): {bad[:10]}")

    target = resolve_target(args.target, metadata, gene_index)
    seed = metadata.get("seed")
    logger.info(
        "Scoring %d candidate genes against target %s. Methods: %s. Graph seed: %s "
        "(no stochastic step in topology/RWR — recorded, not consumed).",
        len(gene_index) - 1, target, ",".join(methods), seed,
    )

    graph = build_digraph(weight, gene_index)

    candidates = [g for g in gene_index if g != target]
    if not candidates:
        logger.warning("Graph has no nodes other than the target — writing an empty ranking.")
    result = pd.DataFrame({"gene": pd.Series(candidates, dtype="object")})

    topo = topology_scores(graph, target)
    result["topology_score"] = result["gene"].map(topo).astype(float)

    rank_col = "topology_score"
    if "rwr" in methods:
        rwr = rwr_scores(graph, gene_index, target, metadata, args.restart_prob, args.node_weight_mix)
        result["rwr_score"] = result["gene"].map(rwr).astype(float)
        rank_col = "rwr_score"

    result = result.sort_values(rank_col, ascending=False, kind="stable").reset_index(drop=True)
    result = CandidateScoresSchema.validate(result)

    out_path = args.out or (args.graph_dir / "candidate_scores.tsv")
    result.to_csv(out_path, sep="\t", index=False)

    top = result.head(5).to_string(index=False)
    logger.info("Ranked by %s (descending). Top %d:\n%s", rank_col, min(5, len(result)), top)
    logger.info("Stage 5 complete: %s", out_path)
    return out_path


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    run(parse_args(argv))


if __name__ == "__main__":
    main()
