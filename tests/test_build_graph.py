"""pytest tests for stage3_build_graph.py (Stage 3: pathway-based gene-gene graph
construction).

Hand-built DataFrames, calling stage3_build_graph.py's functions directly (no
manifest.json/parquet round-trip — that's Stage 1/2 plumbing, already
covered by tests/test_ingest_*.py, tests/test_gsea_discovery.py, and
stage3_build_graph.py's own ``_resolve_artifact_path``). Gene IDs are fake but
Ensembl-pattern-valid (``graph_to_sparse``'s index validation requires it).

Pathways:
- ``R-SMALL`` = {A, B} (2 genes) -> co-membership weight 1/(2-1) = 1.0.
- ``R-LARGE`` = {A, B, C, D} (4 genes) -> weight 1/(4-1) = 1/3 per pair.
  (A, B) sits in *both*, so its total weight is the sum, 1.0 + 1/3 = 4/3 —
  the pair that actually exercises "summed over the pathways that induced
  them," not just the per-pathway formula in isolation.
- ``R-ISOLATE`` = {F} alone (1 gene) — no other gene shares a pathway with
  F, so it must come out of graph construction as an isolated node, not be
  silently dropped for having no edges.

Union membership: R-SMALL/R-LARGE are Stage 2 "disease-relevant" pathways;
R-ISOLATE is the target's own pathway (target = F) and is deliberately
*not* disease-relevant, exercising the union's other half.

Interactions (each deliberately touching only one direction of a
co-membership pair, so the untouched reverse direction proves
co-membership's sign=0 default survives independently):
- A -> B: sign=+1, confidence=0.9 — overrides (A, B) only, not (B, A).
- C -> D: sign=-1, confidence=0.7 — overrides (C, D) only, not (D, C).
- B -> C: sign=0, confidence=0.5 — sign is numerically unchanged from
  co-membership's default, but confidence/evidence_type must still update;
  catches a bug class where a sign==0 interaction row gets skipped as if
  it carried no information.
- Z -> A: Z is not a member of any pathway at all, so this row must be
  dropped as out of scope — which also makes the "exactly these edges"
  assertion below a real test of that scoping, not just an incidental pass.

Expected result: exactly 12 directed co-membership-derived edges (the 6
unordered pairs from R-SMALL ∪ R-LARGE, x2 directions each), F present
with degree 0, Z absent from the graph entirely.
"""

from __future__ import annotations

import pandas as pd
import pytest

import stage3_build_graph as bg

A, B, C, D, F, Z = (f"ENSG{n:011d}" for n in (1, 2, 3, 4, 6, 99))
TARGET_GENE_ID = F

GENE_SETS = pd.DataFrame(
    [
        {"set_id": "R-SMALL", "gene_id": A, "source_db": "reactome"},
        {"set_id": "R-SMALL", "gene_id": B, "source_db": "reactome"},
        {"set_id": "R-LARGE", "gene_id": A, "source_db": "reactome"},
        {"set_id": "R-LARGE", "gene_id": B, "source_db": "reactome"},
        {"set_id": "R-LARGE", "gene_id": C, "source_db": "reactome"},
        {"set_id": "R-LARGE", "gene_id": D, "source_db": "reactome"},
        {"set_id": "R-ISOLATE", "gene_id": F, "source_db": "reactome"},
    ]
)

DISEASE_PATHWAYS = pd.DataFrame({"set_id": ["R-SMALL", "R-LARGE"]})

INTERACTIONS = pd.DataFrame(
    [
        {"gene_a": A, "gene_b": B, "directed": True, "sign": 1, "source_db": "reactome",
         "evidence_type": "curated", "confidence": 0.9},
        {"gene_a": C, "gene_b": D, "directed": True, "sign": -1, "source_db": "reactome",
         "evidence_type": "curated", "confidence": 0.7},
        {"gene_a": B, "gene_b": C, "directed": False, "sign": 0, "source_db": "reactome",
         "evidence_type": "predicted", "confidence": 0.5},
        # Out of scope: Z is not in any pathway, so this row must be dropped.
        {"gene_a": Z, "gene_b": A, "directed": True, "sign": 1, "source_db": "reactome",
         "evidence_type": "curated", "confidence": 0.99},
    ]
)

WEIGHT_SMALL = 1.0 / (2 - 1)  # 1.0
WEIGHT_LARGE = 1.0 / (4 - 1)  # 0.333...

# (gene_a, gene_b) -> expected comembership_weight, sign, confidence, evidence_type, source_db
EXPECTED_EDGES = {
    (A, B): (WEIGHT_SMALL + WEIGHT_LARGE, 1, 0.9, "curated", "reactome"),
    (B, A): (WEIGHT_SMALL + WEIGHT_LARGE, 0, None, None, "reactome"),
    (A, C): (WEIGHT_LARGE, 0, None, None, "reactome"),
    (C, A): (WEIGHT_LARGE, 0, None, None, "reactome"),
    (A, D): (WEIGHT_LARGE, 0, None, None, "reactome"),
    (D, A): (WEIGHT_LARGE, 0, None, None, "reactome"),
    (B, C): (WEIGHT_LARGE, 0, 0.5, "predicted", "reactome"),
    (C, B): (WEIGHT_LARGE, 0, None, None, "reactome"),
    (B, D): (WEIGHT_LARGE, 0, None, None, "reactome"),
    (D, B): (WEIGHT_LARGE, 0, None, None, "reactome"),
    (C, D): (WEIGHT_LARGE, -1, 0.7, "curated", "reactome"),
    (D, C): (WEIGHT_LARGE, 0, None, None, "reactome"),
}


@pytest.fixture(scope="module")
def graph_and_stats():
    return bg.build_pathway_based_gene_gene_graph(GENE_SETS, DISEASE_PATHWAYS, INTERACTIONS, TARGET_GENE_ID)


def test_exact_directed_signed_edge_list_with_source_db(graph_and_stats):
    graph, _stats = graph_and_stats

    actual_edges = {(u, v): (
        attrs["comembership_weight"], attrs["sign"], attrs["confidence"], attrs["evidence_type"], attrs["source_db"]
    ) for u, v, attrs in graph.edges(data=True)}

    assert set(actual_edges) == set(EXPECTED_EDGES), "edge set must match exactly — no extra, none missing"
    for pair, (exp_weight, exp_sign, exp_conf, exp_evtype, exp_source_db) in EXPECTED_EDGES.items():
        weight, sign, conf, evtype, source_db = actual_edges[pair]
        assert weight == pytest.approx(exp_weight), f"{pair}: comembership_weight"
        assert sign == exp_sign, f"{pair}: sign"
        assert conf == exp_conf, f"{pair}: confidence"
        assert evtype == exp_evtype, f"{pair}: evidence_type"
        assert source_db == exp_source_db, f"{pair}: source_db"

    # Z must never have entered the graph at all — the out-of-scope row was dropped.
    assert Z not in graph.nodes


def test_comembership_weights_match_hand_computed_formula(graph_and_stats):
    graph, _stats = graph_and_stats
    # A pair present in only the large pathway: exactly 1/(4-1).
    assert graph[A][C]["comembership_weight"] == pytest.approx(1.0 / 3.0)
    # A pair present in both pathways: the *sum* of both pathways' contributions.
    assert graph[A][B]["comembership_weight"] == pytest.approx(1.0 / (2 - 1) + 1.0 / (4 - 1))


def test_isolated_gene_present_not_dropped(graph_and_stats):
    graph, _stats = graph_and_stats
    assert F in graph.nodes
    assert graph.degree(F) == 0
    assert graph.number_of_nodes() == 5  # A, B, C, D, F — not Z


def test_sparse_serialization_matches_graph(graph_and_stats):
    """Light cross-check that graph_to_sparse's combined weight/sign matrices
    agree with the networkx-level attributes exercised above."""
    graph, _stats = graph_and_stats
    weight, sign, gene_index = bg.graph_to_sparse(graph)
    i, j = gene_index.index(A), gene_index.index(B)
    assert weight[i, j] == pytest.approx(graph[A][B]["comembership_weight"] + graph[A][B]["confidence"])
    assert sign[i, j] == graph[A][B]["sign"]
    assert F in gene_index  # the isolated node must survive serialization too
