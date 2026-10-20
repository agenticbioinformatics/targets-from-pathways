"""pytest tests for stage5_score_candidates.py (Stage 5: propagation from seed genes).

A hand-built 6-node graph with unit edge weights, so every hop costs 1.0
and shortest-path distances are hand-computable. Two seeds, S1 and S2:

    S1 <-> A -> B -> C        (S1<->A is a short cycle; A->B->C a chain)
    S1 -> D <- S2             (D is a direct neighbour of both seeds)

Directed edges: S1->A, A->S1, A->B, B->C, S1->D, S2->D.

Distances from the nearest seed: A=1, D=1, B=2, C=3 (all reachable).

Undirected seed neighbours:
- A: neighbours {S1, B}  -> 1 seed
- D: neighbours {S1, S2} -> 2 seeds
- B: neighbours {A, C}   -> 0
- C: neighbours {B}      -> 0

topology_score = (1 + seed neighbours) / (1 + distance):
- D: 3/2 = 1.5
- A: 2/2 = 1.0
- B: 1/3 ~ 0.3333
- C: 1/4 = 0.25

So the topology ranking is D > A > B > C, and seeds are never candidates.
"""

from __future__ import annotations

import json

import networkx as nx
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import stage5_score_candidates as sc

S1, S2, A, B, C, D = (f"ENSG{n:011d}" for n in (1, 2, 3, 4, 5, 6))
GENE_INDEX = [S1, S2, A, B, C, D]
SEEDS = [S1, S2]
EDGES = {(S1, A): 1.0, (A, S1): 1.0, (A, B): 1.0, (B, C): 1.0, (S1, D): 1.0, (S2, D): 1.0}
EXPECTED_TOPOLOGY = {D: 1.5, A: 1.0, B: 1 / 3, C: 0.25}


def _weight_matrix() -> sp.csr_matrix:
    idx = {g: i for i, g in enumerate(GENE_INDEX)}
    rows, cols, data = [], [], []
    for (u, v), w in EDGES.items():
        rows.append(idx[u])
        cols.append(idx[v])
        data.append(w)
    return sp.csr_matrix((data, (rows, cols)), shape=(len(GENE_INDEX), len(GENE_INDEX)))


def _write_graph_dir(path, seeds_tsv: str | None = None) -> None:
    weight = _weight_matrix()
    sp.save_npz(path / "graph_weight.npz", weight)
    sp.save_npz(path / "graph_sign.npz", sp.csr_matrix(weight.shape, dtype="int8"))
    (path / "graph_gene_index.json").write_text(json.dumps(GENE_INDEX))
    (path / "graph_metadata.json").write_text(json.dumps({"seed": 0}))
    if seeds_tsv is None:
        seeds_tsv = f"gene_id\tscore\n{S1}\t0.9\n{S2}\t0.4\n"
    (path / "disease_seed_genes.tsv").write_text(seeds_tsv)


def _run(tmp_path, method=None, extra=()):
    """Run Stage 5 end to end. ``method=None`` omits ``--method`` entirely,
    exercising the argparse default."""
    argv = ["--graph-dir", str(tmp_path), *extra]
    if method is not None:
        argv += ["--method", method]
    out = sc.run(sc.parse_args(argv))
    return pd.read_csv(out, sep="\t")


# --------------------------------------------------------------------------
# topology
# --------------------------------------------------------------------------


def test_topology_scores_match_hand_computed():
    graph = sc.build_digraph(_weight_matrix(), GENE_INDEX)
    scores = sc.topology_scores(graph, SEEDS)
    assert scores == pytest.approx(EXPECTED_TOPOLOGY)
    assert S1 not in scores and S2 not in scores  # seeds are never candidates


def test_unreachable_node_scores_zero():
    E = f"ENSG{7:011d}"
    dense = np.zeros((7, 7))
    dense[:6, :6] = _weight_matrix().toarray()
    graph = sc.build_digraph(sp.csr_matrix(dense), GENE_INDEX + [E])
    assert sc.topology_scores(graph, SEEDS)[E] == 0.0


# --------------------------------------------------------------------------
# rwr
# --------------------------------------------------------------------------


def test_rwr_restarts_on_seeds_and_favours_their_neighbours():
    graph = sc.build_digraph(_weight_matrix(), GENE_INDEX)
    scores = sc.rwr_scores(graph, SEEDS, restart_prob=0.5)
    assert S1 not in scores and S2 not in scores
    assert all(v >= 0 for v in scores.values())
    # candidates + the seeds' retained mass sum to 1, so the candidates alone sum to < 1
    assert 0.0 < sum(scores.values()) < 1.0
    # D (fed by both seeds) and A (next to S1) beat the tail of the chain
    assert min(scores[D], scores[A]) > scores[C]


def test_rwr_full_distribution_sums_to_one():
    graph = sc.build_digraph(_weight_matrix(), GENE_INDEX)
    restart = {S1: 0.5, S2: 0.5}
    full = nx.pagerank(graph, alpha=0.5, personalization=restart, weight="weight")
    assert sum(full.values()) == pytest.approx(1.0)
    scores = sc.rwr_scores(graph, SEEDS, restart_prob=0.5)
    for gene, score in scores.items():
        assert score == pytest.approx(full[gene])


def test_rwr_rejects_bad_restart_prob():
    graph = sc.build_digraph(_weight_matrix(), GENE_INDEX)
    with pytest.raises(SystemExit):
        sc.rwr_scores(graph, SEEDS, restart_prob=1.0)


# --------------------------------------------------------------------------
# bfwalk (skipped if the bfwalk package is not installed)
# --------------------------------------------------------------------------


def test_bfwalk_follows_edge_direction_and_attenuates():
    pytest.importorskip("BFWalk.bfwalk")
    scores = sc.bfwalk_scores(_weight_matrix(), GENE_INDEX, SEEDS, alpha=0.5, threads=1)
    assert S1 not in scores and S2 not in scores
    assert all(v >= 0 for v in scores.values())
    # the walk spreads S1 -> A -> B -> C, attenuating at each hop
    assert scores[A] > scores[B] > scores[C] > 0
    assert scores[D] > 0


def test_bfwalk_ignores_zero_weight_edges_and_rescales_weights():
    pytest.importorskip("BFWalk.bfwalk")
    weight = (_weight_matrix() * 4.0).tolil()  # weights > 1 are rescaled into ]0, 1]
    weight[GENE_INDEX.index(S1), GENE_INDEX.index(D)] = 0.0  # evidence-free edges carry no flow
    weight[GENE_INDEX.index(S2), GENE_INDEX.index(D)] = 0.0
    scores = sc.bfwalk_scores(weight.tocsr(), GENE_INDEX, SEEDS, alpha=0.5, threads=1)
    assert scores[D] == 0.0
    assert scores[A] > 0


def test_bfwalk_rejects_bad_alpha():
    with pytest.raises(SystemExit):
        sc.bfwalk_scores(_weight_matrix(), GENE_INDEX, SEEDS, alpha=1.5, threads=1)


# --------------------------------------------------------------------------
# run(): seeds input and output contract
# --------------------------------------------------------------------------


def test_default_run_emits_exact_topology_scores_only(tmp_path):
    _write_graph_dir(tmp_path)
    df = _run(tmp_path)  # no --method -> argparse default
    assert list(df.columns) == ["gene", "topology_score"]
    assert set(df["gene"]) == {A, B, C, D}  # seeds excluded
    assert list(df["gene"]) == [D, A, B, C]  # ranked by topology_score
    assert dict(zip(df["gene"], df["topology_score"])) == pytest.approx(EXPECTED_TOPOLOGY)


def test_topology_rwr_emits_both_columns_ranked_by_rwr(tmp_path):
    _write_graph_dir(tmp_path)
    df = _run(tmp_path, "topology,rwr")
    assert list(df.columns) == ["gene", "topology_score", "rwr_score"]
    assert df["rwr_score"].is_monotonic_decreasing


def test_rwr_only_still_emits_topology_score(tmp_path):
    _write_graph_dir(tmp_path)
    df = _run(tmp_path, "rwr")
    assert {"topology_score", "rwr_score"} <= set(df.columns)


def test_topology_rwr_bfwalk_emits_all_columns_ranked_by_bfwalk(tmp_path):
    pytest.importorskip("BFWalk.bfwalk")
    _write_graph_dir(tmp_path)
    df = _run(tmp_path, "topology,rwr,bfwalk")
    assert list(df.columns) == ["gene", "topology_score", "rwr_score", "bfwalk_score"]
    assert df["bfwalk_score"].is_monotonic_decreasing


def test_seed_list_without_score_and_seeds_outside_graph(tmp_path):
    # A hand-made seed list: gene_id only, plus one seed that is not a graph node (dropped).
    outside = f"ENSG{99:011d}"
    _write_graph_dir(tmp_path, seeds_tsv=f"gene_id\n{S1}\n{outside}\n")
    seeds_path = tmp_path / "my_seeds.tsv"
    seeds_path.write_text(f"gene_id\n{S2}\n")
    # --seed-genes overrides the default disease_seed_genes.tsv
    df = _run(tmp_path, extra=["--seed-genes", str(seeds_path)])
    assert S2 not in set(df["gene"]) and S1 in set(df["gene"])
    # the default file: S1 is the only seed in the graph, so S2 becomes a candidate
    df = _run(tmp_path)
    assert S1 not in set(df["gene"]) and S2 in set(df["gene"])


def test_no_seed_in_graph_exits(tmp_path):
    _write_graph_dir(tmp_path, seeds_tsv=f"gene_id\nENSG{99:011d}\n")
    with pytest.raises(SystemExit):
        _run(tmp_path)


def test_missing_seed_file_exits(tmp_path):
    _write_graph_dir(tmp_path)
    (tmp_path / "disease_seed_genes.tsv").unlink()
    with pytest.raises(SystemExit):
        _run(tmp_path)


def test_unknown_method_exits(tmp_path):
    _write_graph_dir(tmp_path)
    with pytest.raises(SystemExit):
        _run(tmp_path, "topology,bfw")
