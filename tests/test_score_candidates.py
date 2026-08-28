"""pytest tests for score_candidates.py (Stage 5: candidate scoring).

Hand-built graph, calling the module's functions directly, plus one
full ``run()`` round-trip through a written graph directory for the
output-contract checks (columns present/absent, schema, target excluded).

Synthetic graph — unit edge weights everywhere, so every hop costs 1.0 and
shortest-path distances are hand-computable:

    T <-> A -> B -> C        (T<->A is the "short cycle"; A->B->C the chain)
    T -> D                   (D is a direct neighbour of T)

Directed edges: T->A, A->T, A->B, B->C, T->D.

Distances from T: A=1, B=2, C=3, D=1  (all reachable).

Undirected common neighbours with T (T's undirected neighbours = {A, D}):
- A: neighbours {T, B}          -> common {} -> 0
- B: neighbours {A, C}          -> common {A} -> 1
- C: neighbours {B}             -> common {} -> 0
- D: neighbours {T}             -> common {} -> 0

topology_score = (1 + common_neighbours) / (1 + distance):
- A: 1/2   = 0.5
- B: 2/3   ~ 0.6667
- C: 1/4   = 0.25
- D: 1/2   = 0.5

So the topology ranking is B > A == D > C.
"""

from __future__ import annotations

import json

import pytest
import scipy.sparse as sp

import score_candidates as sc

T, A, B, C, D = (f"ENSG{n:011d}" for n in (1, 2, 3, 4, 5))
GENE_INDEX = [T, A, B, C, D]
EDGES = {(T, A): 1.0, (A, T): 1.0, (A, B): 1.0, (B, C): 1.0, (T, D): 1.0}


def _weight_matrix() -> sp.csr_matrix:
    idx = {g: i for i, g in enumerate(GENE_INDEX)}
    rows, cols, data = [], [], []
    for (u, v), w in EDGES.items():
        rows.append(idx[u])
        cols.append(idx[v])
        data.append(w)
    return sp.csr_matrix((data, (rows, cols)), shape=(len(GENE_INDEX), len(GENE_INDEX)))


def _write_graph_dir(path, *, node_weights: dict | None = None) -> None:
    weight = _weight_matrix()
    sp.save_npz(path / "graph_weight.npz", weight)
    sp.save_npz(path / "graph_sign.npz", sp.csr_matrix(weight.shape, dtype="int8"))
    (path / "graph_gene_index.json").write_text(json.dumps(GENE_INDEX))
    metadata = {"seed": 0, "resolved_target": {"gene_id": T, "input": "T", "symbol": "T"}}
    if node_weights is not None:
        metadata["genetic_evidence_score"] = node_weights
    (path / "graph_metadata.json").write_text(json.dumps(metadata))


# --------------------------------------------------------------------------
# topology
# --------------------------------------------------------------------------


def test_topology_scores_match_hand_computed():
    graph = sc.build_digraph(_weight_matrix(), GENE_INDEX)
    scores = sc.topology_scores(graph, T)

    assert scores[A] == pytest.approx(0.5)
    assert scores[B] == pytest.approx(2.0 / 3.0)
    assert scores[C] == pytest.approx(0.25)
    assert scores[D] == pytest.approx(0.5)
    assert T not in scores  # target is never its own candidate
    # B (shares neighbour A with T) outranks the equidistant-but-lonelier D.
    assert scores[B] > scores[D] > scores[C]


def test_unreachable_node_scores_zero():
    # An isolated 6th node with no edges at all — reachable from nothing.
    import numpy as np

    E = f"ENSG{6:011d}"
    idx = GENE_INDEX + [E]
    dense = np.zeros((6, 6))
    dense[:5, :5] = _weight_matrix().toarray()
    graph = sc.build_digraph(sp.csr_matrix(dense), idx)
    scores = sc.topology_scores(graph, T)
    assert scores[E] == 0.0


# --------------------------------------------------------------------------
# rwr
# --------------------------------------------------------------------------


def test_rwr_is_a_probability_distribution_and_favours_neighbours():
    graph = sc.build_digraph(_weight_matrix(), GENE_INDEX)
    scores = sc.rwr_scores(
        graph, GENE_INDEX, T, metadata={}, restart_prob=0.5, node_weight_mix=0.5
    )
    assert T not in scores
    assert all(v >= 0 for v in scores.values())
    # candidates + the target's retained restart mass sum to ~1
    assert 0.0 < sum(scores.values()) < 1.0
    # a direct neighbour of the restart node beats the tail of the chain
    assert scores[A] > scores[C]


def test_rwr_node_weight_blend_shifts_mass_toward_high_evidence_gene():
    graph = sc.build_digraph(_weight_matrix(), GENE_INDEX)
    plain = sc.rwr_scores(graph, GENE_INDEX, T, metadata={}, restart_prob=0.5, node_weight_mix=0.5)
    # C is far from T but carries all the genetic evidence.
    blended = sc.rwr_scores(
        graph, GENE_INDEX, T,
        metadata={"genetic_evidence_score": {C: 1.0, A: 0.0, B: 0.0, D: 0.0, T: 0.0}},
        restart_prob=0.5, node_weight_mix=0.5,
    )
    assert blended[C] > plain[C]


# --------------------------------------------------------------------------
# output contract (run() round-trip)
# --------------------------------------------------------------------------


def _run(tmp_path, method):
    args = sc.parse_args(["--graph-dir", str(tmp_path), "--method", method])
    out = sc.run(args)
    import pandas as pd

    return pd.read_csv(out, sep="\t")


def test_default_method_emits_only_topology_score(tmp_path):
    _write_graph_dir(tmp_path)
    df = _run(tmp_path, "topology")
    assert list(df.columns) == ["gene", "topology_score"]
    assert T not in set(df["gene"])  # target excluded
    assert len(df) == 4
    assert df.iloc[0]["gene"] == B  # ranked by topology_score, B on top


def test_topology_rwr_emits_both_columns_ranked_by_rwr(tmp_path):
    _write_graph_dir(tmp_path)
    df = _run(tmp_path, "topology,rwr")
    assert list(df.columns) == ["gene", "topology_score", "rwr_score"]
    assert df["rwr_score"].is_monotonic_decreasing


def test_rwr_only_still_emits_topology_score(tmp_path):
    _write_graph_dir(tmp_path)
    df = _run(tmp_path, "rwr")
    assert "topology_score" in df.columns
    assert "rwr_score" in df.columns


def test_unknown_method_exits(tmp_path):
    _write_graph_dir(tmp_path)
    with pytest.raises(SystemExit):
        _run(tmp_path, "topology,bfw")
