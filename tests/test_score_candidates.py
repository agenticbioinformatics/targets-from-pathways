"""pytest tests for score_candidates.py (Stage 5: candidate scoring).

Two hand-built graphs, both with unit edge weights so every hop costs 1.0
and shortest-path distances are hand-computable:

- a 5-node fixture (below) driving the function-level topology/RWR checks;
- a separate 6-node fixture with a longer directed chain
  (``_write_e2e_graph_dir``) driving the two end-to-end ``run()`` tests
  that assert the exact ``topology_score`` values, the absent-vs-present
  ``rwr_score`` column, and that RWR is a proper probability distribution.

5-node fixture — unit edge weights everywhere, so every hop costs 1.0 and
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

import networkx as nx
import pandas as pd
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


def _run(tmp_path, method=None):
    """Run Stage 5 end to end. ``method=None`` omits ``--method`` entirely,
    exercising the argparse default."""
    argv = ["--graph-dir", str(tmp_path)]
    if method is not None:
        argv += ["--method", method]
    out = sc.run(sc.parse_args(argv))
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


# --------------------------------------------------------------------------
# End-to-end: short cycle + longer chain, hand-computed distances
# --------------------------------------------------------------------------

# A self-contained 6-node graph, unit edge weights (every hop costs 1.0):
#
#   t <-> a           short 2-cycle (t->a, a->t)
#   t  -> d           d: a near neighbour of t (dangling)
#   a  -> b -> c -> e long directed chain
#
# Distances from t: a=1, d=1, b=2, c=3, e=4.
#
# Undirected neighbours: t={a,d}, a={t,b}, b={a,c}, c={b,e}, e={c}, d={t}.
# Common neighbours with t (t's undirected neighbours = {a, d}):
#   a->0, b->1 (shares a), c->0, e->0, d->0
# topology_score = (1 + common) / (1 + distance):
#   a=1/2, b=2/3, c=1/4, e=1/5, d=1/2
_E2E = [f"ENSG{n:011d}" for n in (101, 102, 103, 104, 105, 106)]
_E2E_T, _E2E_A, _E2E_B, _E2E_C, _E2E_E, _E2E_D = _E2E
_E2E_EDGES = {
    (_E2E_T, _E2E_A): 1.0, (_E2E_A, _E2E_T): 1.0,   # short cycle
    (_E2E_T, _E2E_D): 1.0,                            # near neighbour
    (_E2E_A, _E2E_B): 1.0, (_E2E_B, _E2E_C): 1.0, (_E2E_C, _E2E_E): 1.0,  # long chain
}


def _write_e2e_graph_dir(path) -> None:
    idx = {g: i for i, g in enumerate(_E2E)}
    rows = [idx[u] for u, _ in _E2E_EDGES]
    cols = [idx[v] for _, v in _E2E_EDGES]
    weight = sp.csr_matrix((list(_E2E_EDGES.values()), (rows, cols)), shape=(6, 6))
    sp.save_npz(path / "graph_weight.npz", weight)
    sp.save_npz(path / "graph_sign.npz", sp.csr_matrix((6, 6), dtype="int8"))
    (path / "graph_gene_index.json").write_text(json.dumps(_E2E))
    (path / "graph_metadata.json").write_text(
        json.dumps({"seed": 0, "resolved_target": {"gene_id": _E2E_T, "input": "t", "symbol": "t"}})
    )


def test_default_run_emits_exact_topology_scores_only(tmp_path):
    _write_e2e_graph_dir(tmp_path)

    df = _run(tmp_path)  # no --method at all -> argparse default

    assert list(df.columns) == ["gene", "topology_score"]
    assert "rwr_score" not in df.columns
    assert _E2E_T not in set(df["gene"])  # target is not its own candidate

    got = dict(zip(df["gene"], df["topology_score"]))
    assert got == pytest.approx(
        {_E2E_A: 1 / 2, _E2E_B: 2 / 3, _E2E_C: 1 / 4, _E2E_E: 1 / 5, _E2E_D: 1 / 2}
    )
    # target's own graph neighbours (a, d) outrank the far chain nodes (c, e)
    assert min(got[_E2E_A], got[_E2E_D]) > max(got[_E2E_C], got[_E2E_E])


def test_topology_rwr_run_emits_both_and_rwr_is_a_distribution(tmp_path):
    _write_e2e_graph_dir(tmp_path)

    df = _run(tmp_path, "topology,rwr")
    assert list(df.columns) == ["gene", "topology_score", "rwr_score"]

    rwr = dict(zip(df["gene"], df["rwr_score"]))
    assert all(v >= 0 for v in rwr.values())
    # The emitted column omits the target, so it is a sub-distribution;
    # reconstruct the full stationary distribution the same way score_candidates
    # does and check *that* is a proper probability distribution.
    graph = sc.build_digraph(sp.load_npz(tmp_path / "graph_weight.npz").tocsr(), _E2E)
    personalization = sc._restart_distribution({}, _E2E, _E2E_T, node_weight_mix=0.5)
    full = nx.pagerank(graph, alpha=0.5, personalization=personalization, weight="weight")
    assert sum(full.values()) == pytest.approx(1.0)
    assert all(v >= 0 for v in full.values())
    for gene, score in rwr.items():
        assert score == pytest.approx(full[gene])

    # target's own graph neighbours (a, d) outrank the far chain nodes (c, e)
    assert min(rwr[_E2E_A], rwr[_E2E_D]) > max(rwr[_E2E_C], rwr[_E2E_E])
