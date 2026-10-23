"""pytest tests for validation/leave_one_out.py and validation/compare_scores.py on a hand-built catalogue.

STAR: a hub H linked both ways to the leaves S1, S2, S3 and X, plus an isolated gene Y. Seeds S1, S2, S3, and Z,
which is not in the graph (so it cannot be held out). Holding out S1 leaves the seeds S2 and S3: the candidates
are H, X, Y and S1. H is closest to the seeds; S1 and X are both leaves of H, so they tie; Y gets nothing. S1's
mid-rank is 1 + 1 (H) + 0.5 (the tie with X) = 2.5 of 4 for both methods, so its AUROC is (4 - 2.5) / 3 = 0.5,
and the same holds for S2 and S3 by symmetry. Both methods agree, so every rank aggregation keeps that order and
S1 is again 2.5 of 4. AUPRC: all three held-out genes are recovered at depth 2.5, where the top-2.5 lists hold
3 x 2.5 candidates, so precision 0.4 at recall 1; the prevalence is 3 / 12 = 0.25.
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
    for name in ["scores_left_out_details.tsv", "metrics.tsv", "leave_one_out.log", "roc_pr_curve.png",
                 "compare_scores.log", "compare_scores_details.tsv", "roc_curve_vs_open_targets.png"]:
        assert (star / name).exists(), name
    for name in ["leave_one_out_summary.tsv", "leave_one_out.log", "roc_curve_all_diseases.png",
                 "roc_pr_curve_combined_all_diseases.png", "combined_score_tests.tsv",
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
    for combination in loo.COMBINATIONS:
        assert (details[f"{combination}_rank"] == 2.5).all() and result[f"{combination}_rank"] == 2.5


def test_metrics_and_summary(results):
    summary = pd.read_csv(results / "leave_one_out_summary.tsv", sep="\t").set_index("ranking")
    assert list(summary.index) == loo.RANKINGS + ["random"]
    assert list(summary["diseaseId"].unique()) == ["STAR"]
    assert (summary["n_held_out"] == 3).all() and (summary["n_candidates"] == 4).all()
    for ranking in loo.RANKINGS:
        assert summary.loc[ranking, "auroc"] == pytest.approx(0.5), ranking
        assert summary.loc[ranking, "auprc"] == pytest.approx(0.4), ranking
        assert summary.loc[ranking, "auprc_fold_enrichment"] == pytest.approx(0.4 / 0.25), ranking
    # Depth 20 holds all 4 candidates of each run: precision@20 is the prevalence, no enrichment.
    assert summary.loc["bfwalk", "precision@20"] == pytest.approx(0.25)
    assert summary.loc["bfwalk", "precision@20_fold_enrichment"] == pytest.approx(1.0)
    assert summary.loc["bfwalk", "recall@10"] == 1.0  # 4 candidates: every rank is within the top 10
    random = summary.loc["random"]
    assert random["auroc"] == 0.5 and random["recall@10"] == 1.0 and random["precision@20_fold_enrichment"] == 1.0
    # With 3 held-out genes a random ranking's AUPRC is well above the prevalence (see random_auprc).
    assert random["auprc"] > 0.25 and random["auprc_fold_enrichment"] == pytest.approx(random["auprc"] / 0.25)
    metrics = pd.read_csv(results / "STAR_leave_one_out" / "metrics.tsv", sep="\t", index_col="ranking")
    assert list(metrics.index) == loo.RANKINGS + ["random"]
    tests = pd.read_csv(results / "combined_score_tests.tsv", sep="\t")
    assert set(tests["diseaseId"]) == {"all", "STAR"}
    assert list(tests.columns) == ["diseaseId", "measure", "combination", "n_held_out", "mean", "mean_bfwalk",
                                   "mean_rwr", "mean_random", "p_vs_bfwalk", "p_vs_rwr", "p_vs_random"]
    # Every combination ranks each held-out gene exactly as BFWalk and RWR do: nothing to test.
    assert tests["p_vs_bfwalk"].isna().all() and tests["p_vs_rwr"].isna().all()


def test_compare_scores_ranks_the_open_targets_baseline_without_genetic_evidence(results):
    details = pd.read_csv(results / "STAR_leave_one_out" / "compare_scores_details.tsv", sep="\t").set_index("gene")
    assert details.loc[S1, "ot_score"] == 0.9  # literature, not the 0.99 genetic-association row
    assert details.loc[S1, "ot_rank"] == 1.0   # 0.9 beats X (0.2), H and Y (0)
    assert details.loc[S2, "ot_rank"] == 3.0   # 0 ties with H and Y, below X: 1 + 1 + 0.5 * 2
    summary = pd.read_csv(results / "compare_scores_summary.tsv", sep="\t")
    assert set(summary["method"]) == {"bfwalk", "rwr"}
    assert set(summary["ot_diseaseId"]) == {STAR_OT_ID}
    assert set(summary["verdict"]) <= {"better", "worse", "no_difference"}


def test_roc_area_is_the_mean_auroc():
    ranks, n_candidates = [1, 3, 5], [5, 5, 5]
    expected = loo.metrics_table({"x": ranks}, n_candidates, []).loc["x", "auroc"]
    # The plotted curve's area: 1 - mean((rank - 1) / (n - 1)).
    assert expected == pytest.approx(1 - sum((r - 1) / 4 for r in ranks) / 3)
    assert loo.roc_curve(ranks, n_candidates)[2] == pytest.approx(expected)


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


def test_combinations():
    """G is second for both methods, each of which puts a different gene first: bfwalk ranks A 1, G 2, C 3, B 4;
    rwr B 1, G 2, C 3, A 4.
    - rrf: A and B 1/61 + 1/64 = 0.03202, G 2/62 = 0.03226, C 2/63: G first.
    - borda: A and B 3 + 0, G 2 + 2, C 1 + 1: G first.
    - rra: normalised ranks A (0.25, 1): rho = min(1 - 0.75^2, 1^2) = 0.4375; G (0.5, 0.5): min(0.75, 0.25) = 0.25;
      C (0.75, 0.75): 0.5625. G first.
    - union (better rank): A and B 1, G 2, C 3: A and B first.
    - intersection (worse rank): G 2, C 3, A and B 4: G first, A and B last."""
    scores = pd.DataFrame({"bfwalk_score": [0.4, 0.1, 0.2, 0.3], "rwr_score": [1, 400, 200, 300]},
                          index=pd.Index(["A", "B", "C", "G"], name="gene"))
    combined = loo.combined_scores(scores)
    assert set(combined) == set(loo.COMBINATIONS)
    assert combined["rrf"].to_dict() == pytest.approx({"A": 1 / 61 + 1 / 64, "B": 1 / 61 + 1 / 64, "C": 2 / 63, "G": 2 / 62})
    assert combined["borda"].to_dict() == {"A": 3.0, "B": 3.0, "C": 2.0, "G": 4.0}
    assert (-combined["rra"]).to_dict() == pytest.approx({"A": 0.4375, "B": 0.4375, "C": 0.5625, "G": 0.25})
    for combination in ["rrf", "borda", "rra", "intersection"]:
        assert loo.mid_rank(combined[combination], "G") == 1.0, combination
    assert loo.mid_rank(combined["rrf"], "C") == 4.0
    assert loo.mid_rank(combined["union"], "A") == 1.5 and loo.mid_rank(combined["union"], "G") == 3.0
    assert loo.mid_rank(combined["intersection"], "A") == 3.5


def test_combined_score_tests():
    """Every combination ranks every held-out gene first, BFWalk and RWR last: significantly better than both,
    by per-gene AUROC (1 vs 0) and average precision vs random (100 / harmonic(100) vs 1 / harmonic(100))."""
    n = 20
    ranked = pd.DataFrame({"n_candidates": [100] * n, "bfwalk_rank": [100.0] * n, "rwr_rank": [100.0] * n,
                           **{f"{combination}_rank": [1.0] * n for combination in loo.COMBINATIONS}})
    rows = {(row["measure"], row["combination"]): row for row in loo.combined_score_tests(ranked, np.random.default_rng(0))}
    assert set(rows) == {(measure, combination) for measure in loo.PER_GENE for combination in loo.COMBINATIONS}
    for row in rows.values():
        assert row["n_held_out"] == n
        assert row["p_vs_bfwalk"] < 1e-5 and row["p_vs_rwr"] < 1e-5 and row["p_vs_random"] < 1e-4
    assert (rows[("auroc", "rrf")]["mean"], rows[("auroc", "rrf")]["mean_bfwalk"]) == (1.0, 0.0)
    harmonic = sum(1 / i for i in range(1, 101))
    assert rows[("average_precision_vs_random", "rrf")]["mean"] == pytest.approx(100 / harmonic)
    assert rows[("average_precision_vs_random", "rrf")]["mean_bfwalk"] == pytest.approx(1 / harmonic)
    # Identical rankings: nothing to test.
    same = ranked.assign(rrf_rank=ranked["bfwalk_rank"])
    assert np.isnan(loo.combined_score_tests(same, np.random.default_rng(0))[0]["p_vs_bfwalk"])


def test_average_precision_vs_random_is_one_for_a_random_ranking():
    rng = np.random.default_rng(0)
    n_candidates = np.array([50] * 20000 + [3000] * 20000)
    assert np.mean(loo.average_precision_vs_random(loo.random_ranks(n_candidates, rng), n_candidates)) == \
        pytest.approx(1.0, abs=0.05)


def test_pr_curve():
    """Ranks 1, 3, 5 of 5 candidates each: at depth 1, 1 of 3 genes among 3 candidates; at 3, 2 among 9; at 5,
    3 among 15. AUPRC = average precision, 1/3 x (1/3 + 2/9 + 1/5), not the trapezoidal area. A run with fewer
    candidates than the depth adds only its own."""
    recall, precision, auprc = loo.pr_curve([1, 3, 5], [5, 5, 5])
    assert list(recall) == pytest.approx([1 / 3, 2 / 3, 1])
    assert list(precision) == pytest.approx([1 / 3, 2 / 9, 1 / 5])
    assert auprc == pytest.approx((1 / 3 + 2 / 9 + 1 / 5) / 3)
    assert loo.pr_curve([1, 5], [2, 10])[1][-1] == pytest.approx(2 / (2 + 5))
    assert loo.prevalence([5, 5, 5]) == pytest.approx(0.2)


def test_pr_curve_enrichment_and_precision_at_k():
    """Two runs, 2 and 10 candidates (prevalence 1/2 and 1/10), ranks 1 and 5. At depth 1: 1 held-out gene where
    random expects 1/2 + 1/10 = 0.6; at depth 5: 2 where random expects 1 + 5/10 = 1.5. Within one disease (the
    same n), the enrichment is just precision x n."""
    recall, enrichment, area = loo.pr_curve([1, 5], [2, 10], enrichment=True)
    assert list(enrichment) == pytest.approx([1 / 0.6, 2 / 1.5])
    assert area == pytest.approx(0.5 / 0.6 + 0.5 * 2 / 1.5)
    assert loo.pr_curve([1, 3, 5], [5, 5, 5], enrichment=True)[2] == pytest.approx(5 * loo.pr_curve([1, 3, 5], [5] * 3)[2])
    assert loo.precision_at_k([1, 3, 5], [5, 5, 5], 3) == pytest.approx(2 / 9)
    assert loo.precision_at_k([1, 3, 5], [5, 5, 5], 3, enrichment=True) == pytest.approx(2 / (3 * 3 / 5))
    random = loo.metrics_table({}, [5, 5, 5], [3], follow_up_k=3).loc["random"]
    assert random["precision@3"] == pytest.approx(0.2) and random["precision@3_fold_enrichment"] == 1.0


def test_random_auprc():
    """One held-out gene among n candidates: a random ranking's AUPRC is E[1 / rank] = harmonic(n) / n, more than
    twice the prevalence 1 / n for n = 4; with many held-out genes it approaches the prevalence."""
    rng = np.random.default_rng(0)
    assert loo.random_auprc([4], rng) == pytest.approx((1 + 1 / 2 + 1 / 3 + 1 / 4) / 4, abs=0.06)
    assert loo.random_auprc([1000] * 2000, rng) == pytest.approx(1 / 1000, rel=0.05)


def test_top_k_size_matches_a_direct_sum():
    n_candidates = np.array([3.0, 10.0, 10.0, 50.0])
    for k in [1, 2.5, 3, 7, 10, 60]:
        assert loo.top_k_size(k, n_candidates)[0] == pytest.approx(np.minimum(k, n_candidates).sum())
        assert loo.top_k_size(k, n_candidates, True)[0] == pytest.approx((np.minimum(k, n_candidates) / n_candidates).sum())


def test_plots(tmp_path):
    ranks = {ranking: [1, 3, 5] for ranking in loo.RANKINGS}
    loo.plot_roc_pr(ranks, [5, 5, 5], tmp_path / "disease.png", "test", follow_up_k=2)
    loo.plot_roc_pr(ranks, [5, 5, 5], tmp_path / "combined.png", "test", follow_up_k=2, enrichment=True)
    loo.plot_roc({"bfwalk": [1, 3, 5], "open_targets": [2, 2, 2]}, [5, 5, 5], tmp_path / "roc.png", "test")
    for name in ["disease.png", "combined.png", "roc.png"]:
        assert (tmp_path / name).stat().st_size > 0


def test_cached_results_get_up_to_date_combinations(tmp_path):
    """On a cache hit, the combinations' ranks are recomputed from candidate_scores.tsv, without scoring again:
    missing ones (cached before they existed) and stale ones (cached under an earlier definition)."""
    catalogue, out_dir = tmp_path / "catalogue", tmp_path / "out"
    _write_disease(catalogue, "STAR", _star_graph(), {S1: 0.9, S2: 0.5, S3: 0.3})
    arguments = ["--catalogue", str(catalogue), "--out-dir", str(out_dir)]
    assert loo.main(arguments) == 0
    details_file = out_dir / "STAR_leave_one_out" / "scores_left_out_details.tsv"
    pd.read_csv(details_file, sep="\t").drop(columns="rra_rank").assign(rrf_rank=99.0).to_csv(details_file, sep="\t", index=False)
    assert loo.main(arguments) == 0
    assert "Cached:" in (out_dir / "STAR_leave_one_out" / "leave_one_out.log").read_text()
    details = pd.read_csv(details_file, sep="\t")
    assert (details["rra_rank"] == 2.5).all() and (details["rrf_rank"] == 2.5).all()
