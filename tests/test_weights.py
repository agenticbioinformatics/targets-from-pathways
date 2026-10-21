"""pytest tests for stage4_weights.py (Stage 4: optional user
node/edge weights, and genetic-evidence gene scores).

Hand-built inputs, calling ``compute_genetic_evidence_scores``,
``apply_user_weights`` and ``read_weights_file`` directly — no
manifest.json/graph-directory round-trip (that's Stage 1/3 plumbing,
already covered by tests/test_ingest_*.py, tests/test_build_graph.py, and
stage0_schemas.py's ``find_artifact``/``load_graph``). Genes A-D stand for
the graph's gene index; the user-weight tests use a small chain graph
A<->B, B->C, C->D (see the second half of the file).

Open Targets association fixture (a mix of datatypes across genes A-D):
- A: ``genetic_association``=0.9, ``known_drug``=0.5, and a decoy
  ``affected_pathway``=0.99 — the max of the *allowed* rows must win
  (0.9), never the higher disallowed one.
- B: ``genetic_association``=0.7 and a decoy ``known_drug``=0.95 — only
  genetic_association counts, so B scores 0.7.
- C: ``affected_pathway``=0.99 only — the case the task calls out
  explicitly: a gene with *only* excluded-datatype evidence must still get
  genetic_evidence_score 0.0 by default, not 0.99 and not null.
- D: no rows at all — the plain "no evidence" case, kept distinct from C's
  "has evidence, wrong datatype" case.

Expected node scores (genetic_association only): A=0.9, B=0.7, C=0.0, D=0.0.
"""

from __future__ import annotations

import pandas as pd
import pytest
import scipy.sparse as sp

import stage4_weights as wts

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
        {"gene_id": B, "disease_id": "EFO_TEST", "datatype_id": "genetic_association",
         "datasource_id": "gwas_credible_sets", "score": 0.7},
        # Decoy: B's known_drug evidence is higher than its genetic_association score.
        {"gene_id": B, "disease_id": "EFO_TEST", "datatype_id": "known_drug",
         "datasource_id": "clinical_precedence", "score": 0.95},
        # C has *only* excluded-datatype evidence.
        {"gene_id": C, "disease_id": "EFO_TEST", "datatype_id": "affected_pathway",
         "datasource_id": "reactome", "score": 0.99},
        # D has no rows at all.
    ]
)


def test_node_weights_reflect_only_allowed_datatypes():
    scores = wts.compute_genetic_evidence_scores(OT_DISEASE_SUBSET, wts.EVIDENCE_DATATYPES, GENE_INDEX)
    assert scores[A] == pytest.approx(0.9)  # max of allowed rows, not the 0.99 decoy
    assert scores[B] == pytest.approx(0.7)
    assert scores[C] == pytest.approx(0.0)
    assert scores[D] == pytest.approx(0.0)

    # The filter is genuinely datatype-driven, not hardcoded: flip the
    # allowlist to affected_pathway only and the picture inverts.
    scores_pathway_only = wts.compute_genetic_evidence_scores(OT_DISEASE_SUBSET, ["affected_pathway"], GENE_INDEX)
    assert scores_pathway_only[A] == pytest.approx(0.99)
    assert scores_pathway_only[B] == pytest.approx(0.0)
    assert scores_pathway_only[C] == pytest.approx(0.99)
    assert scores_pathway_only[D] == pytest.approx(0.0)


def test_gene_with_only_affected_pathway_evidence_gets_zero_weight():
    scores = wts.compute_genetic_evidence_scores(OT_DISEASE_SUBSET, wts.EVIDENCE_DATATYPES, GENE_INDEX)
    assert scores[C] == 0.0


def test_known_drug_evidence_is_never_counted():
    # Only genetic_association counts: known_drug would favour genes already drugged for the disease.
    assert wts.EVIDENCE_DATATYPES == ["genetic_association"]
    known_drug_only = pd.DataFrame([{"gene_id": D, "disease_id": "EFO_TEST", "datatype_id": "known_drug",
                                     "datasource_id": "clinical_precedence", "score": 0.8}])
    scores = wts.compute_genetic_evidence_scores(known_drug_only, wts.EVIDENCE_DATATYPES, GENE_INDEX)
    assert scores[D] == 0.0


# --------------------------------------------------------------------------
# optional user node / edge weights
# --------------------------------------------------------------------------

Z = f"ENSG{99:011d}"  # not a graph gene


def _chain_graph() -> sp.csr_matrix:
    # A<->B, B->C, C->D, all weight 1 (Stage 3's unweighted graph).
    rows, cols = [0, 1, 1, 2], [1, 0, 2, 3]
    return sp.csr_matrix(([1.0] * 4, (rows, cols)), shape=(4, 4))


NODE_WEIGHTS = pd.DataFrame({"gene_id": [A, C, Z], "weight": [0.5, 2.0, 7.0]})  # B, D default
EDGE_WEIGHTS = pd.DataFrame({
    "gene_a": [B, D, A, A],
    "gene_b": [A, C, D, Z],
    "weight": [0.4, 3.0, 9.0, 5.0],  # (B, A) covers A->B and B->A; (A, D) is no edge; Z is no gene
})


def _weights(matrix):
    i = {gene: n for n, gene in enumerate(GENE_INDEX)}
    dense = matrix.toarray()
    return {(a, b): dense[i[a], i[b]] for a, b in [(A, B), (B, A), (B, C), (C, D)]}


def test_without_user_weights_every_edge_is_one():
    weight, stats = wts.apply_user_weights(_chain_graph() * 5.0, GENE_INDEX, None, None, "mean", 1.0, 1.0)
    # Reset to Stage 3's structure: an earlier weighting never compounds.
    assert set(weight.data) == {1.0} and weight.nnz == 4
    assert stats["n_node_weights_used"] == 0 and stats["n_edge_weights_used"] == 0


def test_node_weights_combine_endpoints():
    weight, stats = wts.apply_user_weights(_chain_graph(), GENE_INDEX, NODE_WEIGHTS, None, "mean", 1.0, 1.0)
    assert _weights(weight) == pytest.approx({(A, B): 0.75, (B, A): 0.75, (B, C): 1.5, (C, D): 1.5})
    assert stats["n_node_weights_used"] == 2 and stats["n_node_weights_not_in_graph"] == 1

    weight, _ = wts.apply_user_weights(_chain_graph(), GENE_INDEX, NODE_WEIGHTS, None, "product", 1.0, 1.0)
    assert _weights(weight) == pytest.approx({(A, B): 0.5, (B, A): 0.5, (B, C): 2.0, (C, D): 2.0})

    # A gene missing from the file weighs --default-node-weight.
    weight, _ = wts.apply_user_weights(_chain_graph(), GENE_INDEX, NODE_WEIGHTS, None, "mean", 0.1, 1.0)
    assert _weights(weight)[(B, C)] == pytest.approx((0.1 + 2.0) / 2)


def test_edge_weights_match_either_direction_and_never_add_edges():
    weight, stats = wts.apply_user_weights(_chain_graph(), GENE_INDEX, None, EDGE_WEIGHTS, "mean", 1.0, 0.2)
    assert weight.nnz == 4  # (A, D) is not added
    assert _weights(weight) == pytest.approx({(A, B): 0.4, (B, A): 0.4, (B, C): 0.2, (C, D): 3.0})
    assert stats["n_edge_weights_used"] == 2 and stats["n_edge_weights_not_in_graph"] == 2


def test_node_and_edge_weights_multiply():
    weight, _ = wts.apply_user_weights(_chain_graph(), GENE_INDEX, NODE_WEIGHTS, EDGE_WEIGHTS, "mean", 1.0, 1.0)
    assert _weights(weight) == pytest.approx({(A, B): 0.75 * 0.4, (B, A): 0.75 * 0.4, (B, C): 1.5, (C, D): 1.5 * 3.0})


def test_edge_weights_listing_a_pair_twice_exit():
    twice = pd.DataFrame({"gene_a": [A, B], "gene_b": [B, A], "weight": [0.4, 0.6]})
    with pytest.raises(SystemExit):
        wts.apply_user_weights(_chain_graph(), GENE_INDEX, None, twice, "mean", 1.0, 1.0)


@pytest.mark.parametrize("content", [
    "gene_id\tweight\nENSG00000000001\t0\n",          # weight must be > 0
    "gene_id\tweight\nENSG00000000001\t-1\n",
    "gene_id\tweight\nGENE01\t0.5\n",                 # Ensembl IDs only
    "gene\tweight\nENSG00000000001\t0.5\n",           # wrong column name
    "gene_id\tweight\nENSG00000000001\t0.5\nENSG00000000001\t0.7\n",  # duplicate gene
])
def test_invalid_node_weights_file_exits(tmp_path, content):
    path = tmp_path / "node_weights.tsv"
    path.write_text(content)
    with pytest.raises(SystemExit):
        wts.read_weights_file(path, wts.UserNodeWeightsSchema, "--node-weights")
