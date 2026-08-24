"""pytest tests for genetic_evidence_weights.py (Stage 4: genetic-evidence
node/edge weighting).

Hand-built inputs, calling the module's two pure functions directly — no
manifest.json/graph-directory round-trip (that's Stage 1/3 plumbing,
already covered by tests/test_ingest_*.py, tests/test_build_graph.py, and
this module's own ``_resolve_artifact_path``/``load_stage3_graph``).

Small synthetic graph (a chain, as a raw ``scipy.sparse`` weight matrix +
gene index — exactly the shape ``load_stage3_graph`` hands to
``compute_edge_weights``; the actual numeric values in this matrix are
irrelevant to Stage 4, only *where* it has entries, since edge weights are
derived purely from endpoint node scores):

    A <-> B <-> C <-> D

Open Targets association fixture (a mix of datatypes across genes A-D):
- A: ``genetic_association``=0.9, ``known_drug``=0.5, and a decoy
  ``affected_pathway``=0.99 — the max of the *allowed* rows must win
  (0.9), never the higher disallowed one.
- B: ``known_drug``=0.7 only.
- C: ``affected_pathway``=0.99 only — the case the task calls out
  explicitly: a gene with *only* excluded-datatype evidence must still get
  genetic_evidence_score 0.0 by default, not 0.99 and not null.
- D: no rows at all — the plain "no evidence" case, kept distinct from C's
  "has evidence, wrong datatype" case.

Expected node scores with the default --datatypes (genetic_association,
known_drug): A=0.9, B=0.7, C=0.0, D=0.0.

Expected edge weights (avg mode) at the graph's existing edges: (A,B)=0.8,
(B,C)=0.35, (C,D)=0.0 — the last one a real, present-but-zero edge, not a
missing one, since both C and D score 0.
"""

from __future__ import annotations

import pandas as pd
import pytest
import scipy.sparse as sp

import genetic_evidence_weights as gew

A, B, C, D = (f"ENSG{n:011d}" for n in (1, 2, 3, 4))
GENE_INDEX = [A, B, C, D]

OT_DISEASE_SUBSET = pd.DataFrame(
    [
        {"gene_id": A, "disease_id": "EFO_TEST", "datatype_id": "genetic_association",
         "datasource_id": "eva", "score": 0.9},
        {"gene_id": A, "disease_id": "EFO_TEST", "datatype_id": "known_drug",
         "datasource_id": "clinical_precedence", "score": 0.5},
        # Decoy: A's disallowed-datatype evidence is higher than its allowed max.
        {"gene_id": A, "disease_id": "EFO_TEST", "datatype_id": "affected_pathway",
         "datasource_id": "reactome", "score": 0.99},
        {"gene_id": B, "disease_id": "EFO_TEST", "datatype_id": "known_drug",
         "datasource_id": "clinical_precedence", "score": 0.7},
        # C has *only* excluded-datatype evidence.
        {"gene_id": C, "disease_id": "EFO_TEST", "datatype_id": "affected_pathway",
         "datasource_id": "reactome", "score": 0.99},
        # D has no rows at all.
    ]
)


def _chain_weight_matrix() -> sp.csr_matrix:
    # A<->B<->C<->D, arbitrary structural weights (values don't matter here).
    rows = [0, 1, 1, 2, 2, 3]
    cols = [1, 0, 2, 1, 3, 2]
    data = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    return sp.csr_matrix((data, (rows, cols)), shape=(4, 4))


def test_node_weights_reflect_only_allowed_datatypes():
    scores = gew.compute_genetic_evidence_scores(OT_DISEASE_SUBSET, gew.DEFAULT_DATATYPES, GENE_INDEX)
    assert scores[A] == pytest.approx(0.9)  # max of allowed rows, not the 0.99 decoy
    assert scores[B] == pytest.approx(0.7)
    assert scores[C] == pytest.approx(0.0)
    assert scores[D] == pytest.approx(0.0)

    # The filter is genuinely datatype-driven, not hardcoded: flip the
    # allowlist to affected_pathway only and the picture inverts.
    scores_pathway_only = gew.compute_genetic_evidence_scores(OT_DISEASE_SUBSET, ["affected_pathway"], GENE_INDEX)
    assert scores_pathway_only[A] == pytest.approx(0.99)
    assert scores_pathway_only[B] == pytest.approx(0.0)
    assert scores_pathway_only[C] == pytest.approx(0.99)
    assert scores_pathway_only[D] == pytest.approx(0.0)


def test_gene_with_only_affected_pathway_evidence_gets_zero_weight_by_default():
    scores = gew.compute_genetic_evidence_scores(OT_DISEASE_SUBSET, gew.DEFAULT_DATATYPES, GENE_INDEX)
    assert scores[C] == 0.0


def test_edge_weights_computed_from_endpoint_node_weights():
    scores = gew.compute_genetic_evidence_scores(OT_DISEASE_SUBSET, gew.DEFAULT_DATATYPES, GENE_INDEX)
    node_scores = scores.to_numpy()
    weight_matrix = _chain_weight_matrix()

    avg_weights = gew.compute_edge_weights(weight_matrix, node_scores, "avg").tocsr()
    i, j, k, l = (GENE_INDEX.index(g) for g in (A, B, C, D))
    assert avg_weights[i, j] == pytest.approx((0.9 + 0.7) / 2)  # (A, B)
    assert avg_weights[j, k] == pytest.approx((0.7 + 0.0) / 2)  # (B, C)
    # Both endpoints score 0, but the edge must still be a real, present
    # zero entry, not silently absent.
    assert (k, l) in set(zip(*weight_matrix.tocoo().nonzero()))
    assert avg_weights[k, l] == pytest.approx(0.0)

    product_weights = gew.compute_edge_weights(weight_matrix, node_scores, "product").tocsr()
    assert product_weights[i, j] == pytest.approx(0.9 * 0.7)
