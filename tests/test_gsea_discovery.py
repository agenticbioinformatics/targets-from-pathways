"""Tests for stage2_gsea_discovery.py (Stage 2: GSEA disease-pathway discovery).

Synthetic DataFrames, calling Stage 2's functions directly; manifest plumbing is covered by the ingest tests and
the end-to-end test. GSEA itself (gsea_core.py) is deterministic and checked against the DisPathCatalogue v2
code in test_gsea_core.py, so these tests can assert exact outcomes.

Fixture: genes G1 (strongest genetic evidence) to G25, scored linearly, plus Z1-Z40 with no evidence (score 0 in
the ranking, as every background gene without evidence).
- strong_pathway = G1-G8: the top of the ranking.
- borderline_pathway = G11-G15: p < 0.05 on its own, but not after BH-FDR over the tested sets.
- bottom_pathway = Z1-Z8: no evidence at all, so its ES is negative.
- noise1-noise8: two scored genes and three without evidence each, scattered.
"""

import pandas as pd
import pytest

import gsea_core
import stage2_gsea_discovery as gd

FDR_THRESHOLD = 0.05

GENES = [f"G{i}" for i in range(1, 26)]
NO_EVIDENCE = [f"Z{i}" for i in range(1, 41)]
_SCORES = {gid: round(1.0 - i * (0.98 / 24), 4) for i, gid in enumerate(GENES)}
STRONG_PATHWAY = GENES[0:8]
BORDERLINE_PATHWAY = GENES[10:15]
BOTTOM_PATHWAY = NO_EVIDENCE[0:8]
NOISE_SETS = {f"noise{j}": [GENES[(j * 3) % 25], GENES[(j * 7 + 2) % 25], NO_EVIDENCE[j], NO_EVIDENCE[j + 10],
                            NO_EVIDENCE[j + 20]] for j in range(1, 9)}


def _gene_sets_df(specs, names=None):
    names = names or {}
    rows = [{"set_id": set_id, "set_name": names.get(set_id, set_id), "source_db": "reactome", "gene_id": gene_id}
            for set_id, gene_ids in specs.items() for gene_id in gene_ids]
    return pd.DataFrame(rows)


def _all_gene_sets():
    return _gene_sets_df({"strong_pathway": STRONG_PATHWAY, "borderline_pathway": BORDERLINE_PATHWAY,
                          "bottom_pathway": BOTTOM_PATHWAY, **NOISE_SETS},
                         {"strong_pathway": "Strongly Enriched Pathway", "borderline_pathway": "Borderline Pathway"})


def _ot_rows(scores, datasource="eva", datatype="genetic_association"):
    return [{"gene_id": gene, "disease_id": "EFO_TEST", "datatype_id": datatype, "datasource_id": datasource,
             "score": score} for gene, score in scores.items()]


def _signature():
    return gd.build_signature(pd.DataFrame(_ot_rows(_SCORES)))


@pytest.fixture(scope="module")
def gsea_result():
    return gd.run_gsea(_signature(), _all_gene_sets(), seed=0)


@pytest.fixture(scope="module")
def significant(gsea_result):
    return gd.filter_results(gsea_result, FDR_THRESHOLD)


def test_deliberately_enriched_pathway_is_significant(gsea_result, significant):
    strong = gsea_result.set_index("set_id").loc["strong_pathway"]
    assert strong.es > 0 and strong.fdr < FDR_THRESHOLD
    assert list(significant["set_id"]) == ["strong_pathway"]
    assert significant.iloc[0]["gene_set"] == "Strongly Enriched Pathway"


def test_fdr_correction_makes_borderline_pathway_nonsignificant(gsea_result, significant):
    borderline = gsea_result.set_index("set_id").loc["borderline_pathway"]
    assert borderline.es > 0 and borderline.pval < 0.05
    assert borderline.fdr > FDR_THRESHOLD
    assert "borderline_pathway" not in set(significant["set_id"])


def test_every_set_is_tested_and_genes_without_evidence_are_in_the_ranking(gsea_result):
    assert set(gsea_result["set_id"]) == set(_all_gene_sets()["set_id"])
    # The bottom pathway's genes have no evidence; they are ranked (score 0), so its ES is negative.
    assert gsea_result.set_index("set_id").loc["bottom_pathway", "es"] < 0
    assert list(gsea_result.columns) == gd.GSEA_RESULT_COLUMNS


def test_negative_es_is_never_disease_relevant():
    result = pd.DataFrame({"set_id": ["up", "down"], "gene_set": ["up", "down"], "source_db": "reactome",
                           "es": [0.5, -0.9], "nes": [1.5, -2.0], "pval": [0.001, 0.0001], "fdr": [0.002, 0.0002],
                           "n_perm": 1000, "geneset_size": 10})
    assert list(gd.filter_results(result, 0.05)["set_id"]) == ["up"]


def test_matches_gsea_core_directly(gsea_result):
    """Stage 2 is a thin layer over gsea_core: same scores and library give the same statistics."""
    library = {set_id: sorted(genes) for set_id, genes in _all_gene_sets().groupby("set_id")["gene_id"]}
    direct = gsea_core.run_gsea(_SCORES, library, seed=0, order_seed=0).set_index("set_id")
    stage2 = gsea_result.set_index("set_id")
    pd.testing.assert_series_equal(stage2["es"], direct.loc[stage2.index, "es"])
    pd.testing.assert_series_equal(stage2["pval"], direct.loc[stage2.index, "pval"])


def test_build_signature_ignores_non_genetic_association_rows():
    rows = _ot_rows({"G1": 0.20}) + _ot_rows({"G1": 0.99}, datasource="reactome", datatype="affected_pathway")
    assert gd.build_signature(pd.DataFrame(rows)).set_index("gene").loc["G1", "score"] == 0.20


def test_build_signature_aggregates_datasources_like_open_targets():
    """Two genetic datasources, 0.8 and 0.4: (0.8 / 1 + 0.4 / 4) / (1 + 1 / 4) = 0.72, as Open Targets'
    genetic_association datatype score (DisPathCatalogue's aggregation). The old max would give 0.8."""
    rows = _ot_rows({"G1": 0.4}, datasource="gene_burden") + _ot_rows({"G1": 0.8}, datasource="eva")
    assert gd.build_signature(pd.DataFrame(rows)).set_index("gene").loc["G1", "score"] == pytest.approx(0.72)
    assert gd.datatype_harmonic_sum([0.6]) == pytest.approx(0.6)


def test_build_signature_drops_zero_scores():
    signature = gd.build_signature(pd.DataFrame(_ot_rows({"G1": 0.5, "G2": 0.0})))
    assert list(signature["gene"]) == ["G1"]


def _assert_nothing_tested(result, caplog, reason):
    assert result.empty
    assert list(result.columns) == gd.GSEA_RESULT_COLUMNS
    assert gd.filter_results(result, FDR_THRESHOLD).empty
    assert reason in caplog.text and "disease_pathways.tsv is written empty" in caplog.text


def test_too_few_genes_with_evidence_skips_gsea(caplog, monkeypatch):
    def must_not_run(*args, **kwargs):
        raise AssertionError("GSEA must not run")
    monkeypatch.setattr(gd.gsea_core, "run_gsea", must_not_run)
    result = gd.run_gsea(_signature().head(gsea_core.MIN_SIZE - 1), _all_gene_sets(), seed=0)
    _assert_nothing_tested(result, caplog, f"Only {gsea_core.MIN_SIZE - 1} gene(s) with genetic evidence")


def test_genes_outside_every_gene_set_do_not_count(caplog):
    outside = gd.build_signature(pd.DataFrame(_ot_rows({f"X{i}": 0.9 for i in range(20)})))
    result = gd.run_gsea(outside, _all_gene_sets(), seed=0)
    _assert_nothing_tested(result, caplog, "Only 0 gene(s) with genetic evidence are in a gene set")


def test_no_gene_set_of_testable_size_skips_gsea(caplog):
    tiny = _gene_sets_df({"tiny": GENES[:3]})
    signature = _signature()
    result = gd.run_gsea(signature, tiny, seed=0, min_genes=1)
    _assert_nothing_tested(result, caplog, "GSEA tested nothing")
