"""pytest tests for validation/leave_one_out.py and validation/compare_scores.py on a hand-built catalogue.

STAR: a hub H linked both ways to the leaves S1, S2, S3 and X, plus an isolated gene Y. Seeds S1, S2, S3, and Z,
which is not in the graph (so it cannot be held out). Holding out S1 leaves the seeds S2 and S3: the candidates
are H, X, Y and S1. H is closest to the seeds; S1 and X are both leaves of H, so they tie; Y gets nothing. S1's
mid-rank is 1 + 1 (H) + 0.5 (the tie with X) = 2.5 of 4 for both methods, so its AUROC is (4 - 2.5) / 3 = 0.5,
and the same holds for S2 and S3 by symmetry.
ONE: a single seed in the graph - holding it out leaves nothing to propagate from, so it is skipped.

Open Targets (for compare_scores.py): S1 has literature 0.9 and a decoy genetic-association row 0.99, which must
be ignored; X has literature 0.2. Open Targets lists STAR as MONDO_0000001, with STAR only as a cross-reference
(as 26.06 does with many older EFO IDs), so compare_scores.py must translate the ID.
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
STAR_OT_ID = "MONDO_0000001"


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


def _star_graph():
    star = nx.DiGraph()
    star.add_nodes_from([Y])
    for leaf in (S1, S2, S3, X):
        star.add_edge(H, leaf, sign=0)
        star.add_edge(leaf, H, sign=0)
    return star


@pytest.fixture(scope="module")
def results(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("loo")
    catalogue = tmp / "disease_graph_catalogue"
    _write_disease(catalogue, "STAR", _star_graph(), {S1: 0.9, S2: 0.5, S3: 0.3, Z: 0.8})
    one = nx.DiGraph()
    one.add_edge(H, X, sign=1)
    _write_disease(catalogue, "ONE", one, {H: 0.7})

    out_dir = tmp / "leave_one_out"
    assert loo.main(["--catalogue", str(catalogue), "--out-dir", str(out_dir)]) == 0

    association_dir = tmp / "associations"
    association_dir.mkdir()
    disease_file = tmp / "disease.parquet"
    pd.DataFrame({"id": [STAR_OT_ID, "MONDO_0000002"], "name": ["star", "other"],
                  "dbXRefs": [["STAR", "DOID:1"], None]}).to_parquet(disease_file, index=False)
    pd.DataFrame([
        {"targetId": S1, "diseaseId": STAR_OT_ID, "aggregationType": "datasourceId", "aggregationValue": "europepmc",
         "associationScore": 0.9},
        {"targetId": S1, "diseaseId": STAR_OT_ID, "aggregationType": "datasourceId", "aggregationValue": "eva",
         "associationScore": 0.99},
        {"targetId": X, "diseaseId": STAR_OT_ID, "aggregationType": "datasourceId", "aggregationValue": "europepmc",
         "associationScore": 0.2},
    ]).to_parquet(association_dir / "part-0.parquet", index=False)
    assert cs.main(["--leave-one-out-dir", str(out_dir), "--association-dir", str(association_dir),
                    "--disease-file", str(disease_file)]) == 0
    return out_dir


def test_one_directory_per_disease_and_held_out_seed(results):
    star = results / "STAR_leave_one_out"
    assert sorted(p.name for p in star.iterdir() if p.is_dir()) == sorted(["full", S1, S2, S3])  # not Z
    assert (star / "cache_key.json").exists()
    for name in ["candidate_scores.tsv", "result.json", "seed.log"]:
        assert (star / S1 / name).exists(), name
    for name in ["scores_left_out_details.tsv", "metrics.tsv", "leave_one_out.log", "roc_curve.png",
                 "compare_scores.log", "compare_scores_details.tsv", "roc_curve_vs_open_targets.png"]:
        assert (star / name).exists(), name
    for name in ["leave_one_out_summary.tsv", "leave_one_out.log", "roc_curve_all_diseases.png",
                 "roc_curves_per_disease.png", "roc_curves_per_disease_bfwalk.png", "roc_curves_per_disease_rwr.png",
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
    assert summary.loc["bfwalk", "auroc_random_expected"] == 0.5
    assert summary.loc["bfwalk", "recall@10_random_expected"] == 1.0
    assert not any(column.endswith("p_vs_random") for column in summary.columns)


def test_compare_scores_ranks_the_open_targets_baseline_without_genetic_evidence(results):
    details = pd.read_csv(results / "STAR_leave_one_out" / "compare_scores_details.tsv", sep="\t").set_index("gene")
    assert details.loc[S1, "ot_score"] == 0.9  # literature, not the 0.99 genetic-association row
    assert details.loc[S1, "ot_rank"] == 1.0   # 0.9 beats X (0.2), H and Y (0)
    assert details.loc[S2, "ot_rank"] == 3.0   # 0 ties with H and Y, below X: 1 + 1 + 0.5 * 2
    summary = pd.read_csv(results / "compare_scores_summary.tsv", sep="\t")
    assert set(summary["method"]) == {"bfwalk", "rwr"}
    assert set(summary["ot_diseaseId"]) == {STAR_OT_ID}
    assert set(summary["verdict"]) <= {"better", "worse", "no_difference"}


def test_roc_area_is_the_mean_auroc(tmp_path):
    ranks, n_candidates = [1, 3, 5], [5, 5, 5]
    expected = loo.ranking_metrics(ranks, n_candidates, [])["auroc"]
    # The plotted curve's area: 1 - mean((rank - 1) / (n - 1)).
    assert expected == pytest.approx(1 - sum((r - 1) / 4 for r in ranks) / 3)
    loo.plot_roc({"bfwalk": ranks}, n_candidates, tmp_path / "roc.png", "test")
    assert (tmp_path / "roc.png").stat().st_size > 0


def test_disease_ids_are_translated_by_cross_reference(tmp_path):
    disease_file = tmp_path / "disease.parquet"
    pd.DataFrame({"id": ["MONDO_0005265", "EFO_0000001"], "dbXRefs": [["EFO:0003767", "DOID:0050589"], None]}
                 ).to_parquet(disease_file, index=False)
    mapped = cs.map_disease_ids(["EFO_0003767", "EFO_0000001", "EFO_9999999"], disease_file)
    assert mapped == {"EFO_0003767": "MONDO_0005265", "EFO_0000001": "EFO_0000001", "EFO_9999999": None}


def test_unchanged_inputs_are_not_scored_again(tmp_path):
    """A second run on the same files and parameters reuses the scores; a changed input file, a changed
    parameter or --recompute scores again."""
    catalogue, out_dir = tmp_path / "catalogue", tmp_path / "out"
    _write_disease(catalogue, "STAR", _star_graph(), {S1: 0.9, S2: 0.5, S3: 0.3})
    scores_file = out_dir / "STAR_leave_one_out" / S1 / "candidate_scores.tsv"
    log_file = out_dir / "STAR_leave_one_out" / "leave_one_out.log"
    arguments = ["--catalogue", str(catalogue), "--out-dir", str(out_dir)]

    def run(*extra):
        before = scores_file.stat().st_mtime_ns if scores_file.exists() else None
        assert loo.main(arguments + list(extra)) == 0
        return scores_file.stat().st_mtime_ns != before, "Cached:" in log_file.read_text()

    assert run() == (True, False)
    assert run() == (False, True)
    assert (out_dir / "STAR_leave_one_out" / "metrics.tsv").exists()  # metrics recomputed from the cache
    assert run("--restart-prob", "0.3") == (True, False)
    assert run("--restart-prob", "0.3", "--recompute") == (True, False)
    seeds = catalogue / "STAR_disease_graph" / "disease_seed_genes.tsv"
    seeds.write_text(seeds.read_text().replace("0.3", "0.4"))
    assert run("--restart-prob", "0.3") == (True, False)
    assert run("--restart-prob", "0.3") == (False, True)


@pytest.mark.parametrize("n_diseases", [2, loo.MAX_LEGEND_DISEASES + 1])
def test_per_disease_roc_plots(tmp_path, n_diseases):
    """Each disease's curve has its own AUROC, the pooled one the mean over all held-out seeds; with more
    diseases than the legend lists, the plot is still drawn."""
    pooled = pd.DataFrame({
        "diseaseId": [f"D{d}" for d in range(n_diseases) for _ in range(3)],
        "n_candidates": 5,
        "bfwalk_rank": [1, 2, 3] * n_diseases,
        "rwr_rank": [3, 4, 5] * n_diseases,
    })
    assert loo.roc_curve([1, 2, 3], [5, 5, 5])[2] == pytest.approx(loo.ranking_metrics([1, 2, 3], [5] * 3, [])["auroc"])
    loo.plot_roc_per_disease(pooled, loo.METHODS, tmp_path / "both.png", "test")
    loo.plot_roc_per_disease(pooled, ["rwr"], tmp_path / "rwr.png", "test")
    assert (tmp_path / "both.png").stat().st_size > 0 and (tmp_path / "rwr.png").stat().st_size > 0
