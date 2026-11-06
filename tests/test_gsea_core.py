"""Tests for pipeline/gsea_core.py, the GSEA shared by Stage 2 and the validation catalogue.

tests/fixtures/gsea_reference.json holds ES, p-values, NES, draw counts and FDR computed on a synthetic input by
the code DisPathCatalogue v2 ran: jdhayhurst/blitzgsea v0.0.1's enrichment_score and get_peak_size_adv, and
DisPathCatalogue's permutation_stats and bh_fdr. The port has to reproduce them exactly.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import gsea_core as g

REFERENCE = json.loads((Path(__file__).parent / "fixtures" / "gsea_reference.json").read_text())


@pytest.fixture(scope="module")
def reference_result():
    return g.run_gsea(REFERENCE["scores"], REFERENCE["library"]).set_index("set_id")


def test_reproduces_blitzgsea_fork_and_catalogue_statistics(reference_result):
    expected = pd.DataFrame(REFERENCE["expected"]).T
    assert sorted(reference_result.index) == sorted(expected.index)
    for column in ["es", "pval", "nes", "fdr"]:
        np.testing.assert_allclose(reference_result[column], expected.loc[reference_result.index, column].astype(float),
                                   rtol=1e-12, atol=1e-12, err_msg=column)
    for column in ["n_perm", "geneset_size"]:
        assert (reference_result[column] == expected.loc[reference_result.index, column].astype(int)).all(), column


def test_null_sampler_matches_blitzgsea():
    background = g.ordered_background(REFERENCE["library"])
    _, values = g.build_signature(REFERENCE["scores"], background)
    abs_signature = np.abs(values - values.mean())
    np.testing.assert_allclose(g.null_es(abs_signature, 7, 50, 123), REFERENCE["null_k7_n50_seed123"], rtol=0, atol=1e-7)


def test_enrichment_score_by_hand():
    # Weights 4, 3, 2, 1 (already sorted); hits at positions 0 and 2: running sum +4/6, -1/2, +2/6, -1/2.
    es = g.enrichment_score(np.array([4.0, 3.0, 2.0, 1.0]), [0, 2])
    assert es == pytest.approx(4 / 6)
    # Hits at the bottom: the running sum only falls until the hits, so ES is negative.
    assert g.enrichment_score(np.array([4.0, 3.0, 2.0, 1.0]), [2, 3]) == pytest.approx(-1.0)


def test_top_sets_are_significant_and_bottom_sets_have_negative_es(reference_result):
    top = reference_result.loc[["S00", "S01", "S02", "S03"]]
    assert (top["es"] > 0).all() and (top["fdr"] < 0.05).all()
    assert (reference_result.loc[["S04", "S05"], "es"] < 0).all()
    selected = g.significant_sets(reference_result.reset_index(), 0.05)
    assert set(selected["set_id"]) >= {"S00", "S01", "S02", "S03"}
    assert (selected["es"] > 0).all() and selected["fdr"].is_monotonic_increasing


def test_deterministic():
    first = g.run_gsea(REFERENCE["scores"], REFERENCE["library"])
    second = g.run_gsea(REFERENCE["scores"], REFERENCE["library"])
    pd.testing.assert_frame_equal(first, second)


def test_tie_order_comes_from_order_seed():
    """Genes tied at 0 are ordered by the seeded background permutation, not alphabetically."""
    background = g.ordered_background(REFERENCE["library"], order_seed=0)
    assert background != sorted(background)
    assert background == g.ordered_background(REFERENCE["library"], order_seed=0)
    assert background != g.ordered_background(REFERENCE["library"], order_seed=1)
    genes, values = g.build_signature(REFERENCE["scores"], background)
    zero_block = [gene for gene, value in zip(genes, values) if value == 0]
    assert zero_block == [gene for gene in background if REFERENCE["scores"].get(gene, 0) == 0]


def test_null_cache_reuses_and_does_not_resample(reference_result, monkeypatch):
    cache = {}
    filled = g.run_gsea(REFERENCE["scores"], REFERENCE["library"], null_cache=cache).set_index("set_id")
    pd.testing.assert_frame_equal(filled, reference_result)
    assert set(cache) == set(reference_result["geneset_size"])

    def no_sampling(*args, **kwargs):
        raise AssertionError("a cached size was sampled again")

    monkeypatch.setattr(g, "null_es", no_sampling)
    reused = g.run_gsea(REFERENCE["scores"], REFERENCE["library"], null_cache=cache).set_index("set_id")
    pd.testing.assert_frame_equal(reused, reference_result)


def test_removing_a_top_gene_weakens_its_sets_with_the_cached_null(reference_result):
    """What leave-one-out with a refitted GSEA does: zero one gene's score and reuse the full run's null."""
    cache = {}
    g.run_gsea(REFERENCE["scores"], REFERENCE["library"], null_cache=cache)
    top_gene = max(REFERENCE["scores"], key=REFERENCE["scores"].get)
    sets_with_it = [s for s, members in REFERENCE["library"].items() if top_gene in members]
    assert sets_with_it
    scores = {gene: score for gene, score in REFERENCE["scores"].items() if gene != top_gene}
    refit = g.run_gsea(scores, REFERENCE["library"], null_cache=cache).set_index("set_id")
    for set_id in sets_with_it:
        if reference_result.loc[set_id, "es"] > 0:
            assert refit.loc[set_id, "es"] < reference_result.loc[set_id, "es"]


def test_bh_fdr():
    p = np.array([0.01, 0.04, 0.03, np.nan, 0.5])
    q = g.bh_fdr(p)
    assert np.isnan(q[3])
    # Sorted p 0.01, 0.03, 0.04, 0.5 (m = 4): p * m / rank = 0.04, 0.06, 0.0533, 0.5, then the running minimum
    # from the largest down.
    np.testing.assert_allclose(q[[0, 2, 1, 4]], [0.04, 0.16 / 3, 0.16 / 3, 0.5])


def test_size_limits_and_empty_result():
    library = {"tiny": ["A", "B", "C"], "ok": ["A", "B", "C", "D", "E"], "other": ["F", "G", "H", "I", "J"]}
    result = g.run_gsea({"A": 1.0, "B": 0.5}, library)
    assert list(result["set_id"]) == ["ok", "other"]
    # A set that is the whole background has no misses, so it can't be tested.
    assert g.run_gsea({"A": 1.0}, {"all": ["A", "B", "C", "D", "E"]}).empty
    empty = g.run_gsea({"A": 1.0}, {"tiny": ["A", "B"]})
    assert empty.empty and list(empty.columns) == g.RESULT_COLUMNS
