"""pytest tests for validation/leave_one_out.py and validation/compare_scores.py on a hand-built catalogue.

STAR: a hub H linked both ways to the leaves S1, S2, S3 and X, plus an isolated gene Y. Seeds S1, S2, S3, and Z,
which is not in the graph (so it cannot be held out). Holding out S1 leaves the seeds S2 and S3: the candidates
are H, X, Y and S1. H is closest to the seeds; S1 and X are both leaves of H, so they tie; Y gets nothing. S1's
mid-rank is 1 + 1 (H) + 0.5 (the tie with X) = 2.5 of 4 for both methods, so its AUROC is (4 - 2.5) / 3 = 0.5,
and the same holds for S2 and S3 by symmetry.
ONE: a single seed in the graph - holding it out leaves nothing to propagate from, so it is skipped.

Open Targets (for compare_scores.py): S1 has literature 0.9 and a decoy genetic-association row 0.99, which must
be ignored; X has literature 0.2.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).parent.parent / "validation"))
import compare_scores as cs  # noqa: E402
import leave_one_out as loo  # noqa: E402
import stage3_build_graph as bg  # noqa: E402

H, S1, S2, S3, X, Y, Z = (f"ENSG{n:011d}" for n in (1, 2, 3, 4, 5, 6, 7))


def _write_disease(catalogue, disease, graph, seeds):
    disease_dir = catalogue / f"{disease}_disease_graph"
    disease_dir.mkdir(parents=True)
    weight, sign, gene_index = bg.graph_to_sparse(graph)
    sp.save_npz(disease_dir / "graph_weight.npz", weight)
    sp.save_npz(disease_dir / "graph_sign.npz", sign)
    (disease_dir / "graph_gene_index.json").write_text(json.dumps(gene_index))
    (disease_dir / "graph_metadata.json").write_text("{}")
    pd.DataFrame({"gene_id": list(seeds), "score": list(seeds.values())}).to_csv(
        disease_dir / "disease_seed_genes.tsv", sep="\t", index=False)


@pytest.fixture(scope="module")
def results(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("loo")
    catalogue = tmp / "disease_graph_catalogue"
    star = nx.DiGraph()
    star.add_nodes_from([Y])
    for leaf in (S1, S2, S3, X):
        star.add_edge(H, leaf, sign=0)
        star.add_edge(leaf, H, sign=0)
    _write_disease(catalogue, "STAR", star, {S1: 0.9, S2: 0.5, S3: 0.3, Z: 0.8})
    one = nx.DiGraph()
    one.add_edge(H, X, sign=1)
    _write_disease(catalogue, "ONE", one, {H: 0.7})

    out_dir = tmp / "leave_one_out"
    assert loo.main(["--catalogue", str(catalogue), "--out-dir", str(out_dir), "--n-random", "50"]) == 0

    association_dir = tmp / "associations"
    association_dir.mkdir()
    pd.DataFrame([
        {"targetId": S1, "diseaseId": "STAR", "aggregationType": "datasourceId", "aggregationValue": "europepmc",
         "associationScore": 0.9},
        {"targetId": S1, "diseaseId": "STAR", "aggregationType": "datasourceId", "aggregationValue": "eva",
         "associationScore": 0.99},
        {"targetId": X, "diseaseId": "STAR", "aggregationType": "datasourceId", "aggregationValue": "europepmc",
         "associationScore": 0.2},
    ]).to_parquet(association_dir / "part-0.parquet", index=False)
    assert cs.main(["--leave-one-out-dir", str(out_dir), "--association-dir", str(association_dir),
                    "--n-random", "50"]) == 0
    return out_dir


def test_one_directory_per_disease_and_held_out_seed(results):
    star = results / "STAR_leave_one_out"
    assert sorted(p.name for p in star.iterdir() if p.is_dir()) == sorted(["full", S1, S2, S3])  # not Z
    for name in ["candidate_scores.tsv", "result.json", "seed.log"]:
        assert (star / S1 / name).exists(), name
    for name in ["scores_left_out_details.tsv", "metrics.tsv", "leave_one_out.log", "roc_curve.png",
                 "compare_scores.log", "compare_scores_details.tsv", "roc_curve_vs_open_targets.png"]:
        assert (star / name).exists(), name
    for name in ["leave_one_out_summary.tsv", "leave_one_out.log", "roc_curve_all_diseases.png",
                 "compare_scores_summary.tsv", "compare_scores.log", "roc_curve_vs_open_targets_all_diseases.png"]:
        assert (results / name).exists(), name
    # ONE has a single seed in its graph: skipped, so it has no results.
    assert not (results / "ONE_leave_one_out" / "scores_left_out_details.tsv").exists()


def test_held_out_seed_is_ranked_among_the_candidates(results):
    star = results / "STAR_leave_one_out"
    candidates = pd.read_csv(star / S1 / "candidate_scores.tsv", sep="\t", index_col="gene")
    assert set(candidates.index) == {H, X, Y, S1}  # the remaining seeds S2, S3 are not candidates
    result = json.loads((star / S1 / "result.json").read_text())
    assert result["status"] == "ranked" and result["n_seeds"] == 2 and result["n_candidates"] == 4
    assert result["bfwalk_rank"] == 2.5 and result["rwr_rank"] == 2.5
    details = pd.read_csv(star / "scores_left_out_details.tsv", sep="\t").set_index("gene")
    assert details.loc[S1, "genetic_evidence"] == 0.9
    assert (details["bfwalk_rank"] == 2.5).all() and (details["rwr_rank"] == 2.5).all()


def test_metrics_and_summary(results):
    summary = pd.read_csv(results / "leave_one_out_summary.tsv", sep="\t").set_index("method")
    assert list(summary["diseaseId"].unique()) == ["STAR"]
    assert summary.loc["bfwalk", "auroc"] == pytest.approx(0.5)
    assert summary.loc["rwr", "auroc"] == pytest.approx(0.5)
    assert summary.loc["bfwalk", "recall@10"] == 1.0  # 4 candidates: every rank is within the top 10


def test_compare_scores_ranks_the_open_targets_baseline_without_genetic_evidence(results):
    details = pd.read_csv(results / "STAR_leave_one_out" / "compare_scores_details.tsv", sep="\t").set_index("gene")
    assert details.loc[S1, "ot_score"] == 0.9  # literature, not the 0.99 genetic-association row
    assert details.loc[S1, "ot_rank"] == 1.0   # 0.9 beats X (0.2), H and Y (0)
    assert details.loc[S2, "ot_rank"] == 3.0   # 0 ties with H and Y, below X: 1 + 1 + 0.5 * 2
    summary = pd.read_csv(results / "compare_scores_summary.tsv", sep="\t")
    assert set(summary["method"]) == {"bfwalk", "rwr"}
    assert set(summary["verdict"]) <= {"better", "worse", "no_difference"}


def test_roc_area_is_the_mean_auroc(tmp_path):
    ranks, n_candidates = [1, 3, 5], [5, 5, 5]
    expected = loo.ranking_metrics(ranks, n_candidates, [])["auroc"]
    # The plotted curve's area: 1 - mean((rank - 1) / (n - 1)).
    assert expected == pytest.approx(1 - sum((r - 1) / 4 for r in ranks) / 3)
    loo.plot_roc({"bfwalk": ranks}, n_candidates, tmp_path / "roc.png", "test", np.random.default_rng(0))
    assert (tmp_path / "roc.png").stat().st_size > 0
