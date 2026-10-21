"""pytest tests for stage3_build_graph.py (Stage 3: pathway-based gene-gene graph
construction).

Hand-built DataFrames, calling stage3_build_graph.py's functions directly (no
manifest.json/parquet round-trip — that's Stage 1/2 plumbing, already
covered by tests/test_ingest_*.py, tests/test_gsea_discovery.py, and
stage0_schemas.py's ``find_artifact``). Gene IDs are fake but
Ensembl-pattern-valid (``graph_to_sparse``'s index validation requires it).

Pathways:
- ``R-SMALL`` = {A, B} and ``R-LARGE`` = {A, B, C, D, E} are the Stage 2
  (GSEA) disease-relevant pathways, so A-E are the graph's nodes. E has no
  interaction, so it must come out as an isolated node, not be silently
  dropped for having no edges.
- ``R-SEED`` = {F} contains a seed gene but is not disease-relevant, and
  ``R-OTHER`` = {D, Y} is neither: F and Y must stay out of the graph (even
  though D, a gene of R-LARGE, is in it).

Interactions (the only edges; sharing a pathway adds none):
- A -> B: sign=+1, confidence=0.9 — directed, so (B, A) is not an edge.
- C -> D: sign=-1, confidence=0.7.
- B - C: undirected, sign=0, confidence=0.5 — an edge in both directions.
- Z -> A: Z is in no pathway, so this row is dropped as out of scope.
- F -> A: F is only on the non-disease seed pathway, so this row is dropped too.

Expected result: exactly 4 directed edges (A->B, C->D, B->C, C->B), each with
weight 1 in the serialized matrix (the graph is unweighted).
"""

from __future__ import annotations

import pandas as pd
import pytest

import stage3_build_graph as bg

A, B, C, D, E, F, Z = (f"ENSG{n:011d}" for n in (1, 2, 3, 4, 5, 6, 99))
Y = f"ENSG{98:011d}"

GENE_SETS = pd.DataFrame(
    [
        {"set_id": "R-SMALL", "gene_id": A, "source_db": "reactome"},
        {"set_id": "R-SMALL", "gene_id": B, "source_db": "reactome"},
        {"set_id": "R-LARGE", "gene_id": A, "source_db": "reactome"},
        {"set_id": "R-LARGE", "gene_id": B, "source_db": "reactome"},
        {"set_id": "R-LARGE", "gene_id": C, "source_db": "reactome"},
        {"set_id": "R-LARGE", "gene_id": D, "source_db": "reactome"},
        {"set_id": "R-LARGE", "gene_id": E, "source_db": "reactome"},
        {"set_id": "R-SEED", "gene_id": F, "source_db": "reactome"},
        {"set_id": "R-OTHER", "gene_id": D, "source_db": "reactome"},
        {"set_id": "R-OTHER", "gene_id": Y, "source_db": "reactome"},
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
        # Out of scope: Z is not in any pathway.
        {"gene_a": Z, "gene_b": A, "directed": True, "sign": 1, "source_db": "reactome",
         "evidence_type": "curated", "confidence": 0.99},
        # Out of scope: F is only on a pathway GSEA did not detect.
        {"gene_a": F, "gene_b": A, "directed": True, "sign": 1, "source_db": "reactome",
         "evidence_type": "curated", "confidence": 0.99},
    ]
)

# (gene_a, gene_b) -> expected sign, confidence, evidence_type, source_db
EXPECTED_EDGES = {
    (A, B): (1, 0.9, "curated", "reactome"),
    (C, D): (-1, 0.7, "curated", "reactome"),
    (B, C): (0, 0.5, "predicted", "reactome"),
    (C, B): (0, 0.5, "predicted", "reactome"),
}


@pytest.fixture(scope="module")
def graph_and_stats():
    return bg.build_graph(GENE_SETS, DISEASE_PATHWAYS, INTERACTIONS)


def test_edges_are_exactly_the_in_scope_interactions(graph_and_stats):
    graph, stats = graph_and_stats

    actual_edges = {(u, v): (
        attrs["sign"], attrs["confidence"], attrs["evidence_type"], attrs["source_db"]
    ) for u, v, attrs in graph.edges(data=True)}

    # No co-membership edges: e.g. A and C share R-LARGE but do not interact.
    assert actual_edges == EXPECTED_EDGES
    assert stats["n_interaction_edges_added"] == 4
    assert stats["n_interaction_edges_out_of_scope"] == 2


def test_nodes_are_the_genes_of_the_gsea_pathways_only(graph_and_stats):
    graph, stats = graph_and_stats
    assert set(graph.nodes) == {A, B, C, D, E}  # not F (seed pathway only), Y (other pathway) or Z (no pathway)
    assert graph.degree(E) == 0  # the isolated gene survives
    assert stats["n_disease_pathways"] == 2
    assert stats["n_nodes"] == 5


def test_sparse_serialization_is_unweighted(graph_and_stats):
    graph, _stats = graph_and_stats
    weight, sign, gene_index = bg.graph_to_sparse(graph)
    assert weight.nnz == 4
    assert set(weight.data) == {1.0}  # interaction confidence does not become a weight
    i, j = gene_index.index(A), gene_index.index(B)
    assert weight[i, j] == 1.0 and weight[j, i] == 0.0  # direction kept
    assert sign[i, j] == 1
    assert sign[gene_index.index(C), gene_index.index(D)] == -1
    assert E in gene_index  # the isolated node must survive serialization too


def test_undirected_interaction_connects_both_directions():
    """An undirected interaction between genes of two different pathways must come out symmetric, or
    propagation would cross it only one way, chosen by the FI file's row order."""
    gene_sets = pd.DataFrame([
        {"set_id": "R-P1", "gene_id": A, "source_db": "reactome"},
        {"set_id": "R-P1", "gene_id": C, "source_db": "reactome"},
        {"set_id": "R-P2", "gene_id": B, "source_db": "reactome"},
        {"set_id": "R-P2", "gene_id": D, "source_db": "reactome"},
    ])
    interactions = pd.DataFrame([
        {"gene_a": A, "gene_b": B, "directed": False, "sign": 0, "source_db": "reactome",
         "evidence_type": "curated", "confidence": 1.0},
    ])
    graph, stats = bg.build_graph(gene_sets, pd.DataFrame({"set_id": ["R-P1", "R-P2"]}), interactions)
    weight, _sign, gene_index = bg.graph_to_sparse(graph)
    i, j = gene_index.index(A), gene_index.index(B)
    assert weight[i, j] == pytest.approx(1.0)
    assert weight[j, i] == pytest.approx(1.0)
    assert stats["n_interaction_edges_added"] == 2
