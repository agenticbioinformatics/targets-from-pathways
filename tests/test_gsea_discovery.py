"""pytest tests for stage2_gsea_discovery.py (Stage 2: GSEA disease-pathway discovery).

Everything here is a synthetic, hand-built DataFrame — no manifest.json or
parquet round-trip — calling stage2_gsea_discovery.py's functions directly.
Resolving artifacts out of a real manifest is Stage 1 plumbing, already
covered by tests/test_ingest_*.py and stage2_gsea_discovery.py's own
``_resolve_artifact_path``; what needs covering here is the GSEA/statistics
logic itself.

blitzgsea is *not* perfectly deterministic run-to-run even with a fixed
``seed`` (observed: its internal calibration pool occasionally lands on one
of two nearby outcomes, presumably worker-scheduling-order dependent).
Every assertion below therefore checks a threshold relationship (p-value or
FDR against a cutoff), never an exact float — and the fixture's margins
were chosen empirically, from several real runs, to hold comfortably on
either side of that jitter.

25-gene, 9-pathway shared fixture (``_all_gene_sets`` / ``_signature``):
genes G1 (highest genetic_association score) through G25 (lowest), ranked
linearly.

- ``strong_pathway`` = G1-G8: a clean, unambiguous top-of-ranking cluster —
  the "deliberately enriched" pathway (req. 1).
- ``borderline_pathway`` = G17-G21: a smaller, less extreme cluster,
  individually significant at the nominal p-value threshold on its own, but
  not after BH-FDR is applied across the whole tested library (req. 2) —
  this is the whole point of testing pathways together rather than one at a
  time.
- ``dup_ancestor`` = G9-G14 (6 genes) / ``dup_child`` = G9-G13 (5 genes):
  Jaccard(ancestor, child) = 5/6 = 0.833 > 0.7, so collapsing must drop the
  larger ancestor and keep the smaller, more specific child (req. 4).
- ``noise1``..``noise6``: five genes each, scattered across the ranking
  (no coherent cluster) — non-significant filler that both pads the tested
  library (so BH-FDR in req. 2 has something to correct against) and
  exercises "tested but excluded from output" (req. 3).
"""

from __future__ import annotations

import pandas as pd
import pytest

import stage2_gsea_discovery as gd

PVAL_THRESHOLD = 0.05
FDR_THRESHOLD = 0.1

GENES = [f"G{i}" for i in range(1, 26)]  # G1..G25
_SCORES = {gid: round(1.0 - i * (0.98 / 24), 4) for i, gid in enumerate(GENES)}

STRONG_PATHWAY = GENES[0:8]  # G1-G8
BORDERLINE_PATHWAY = GENES[16:21]  # G17-G21
DUP_ANCESTOR = GENES[8:14]  # G9-G14 (6 genes)
DUP_CHILD = GENES[8:13]  # G9-G13 (5 genes) -- subset of DUP_ANCESTOR
NOISE_SETS = {
    "noise1": [GENES[2], GENES[11], GENES[17], GENES[22], GENES[6]],
    "noise2": [GENES[4], GENES[13], GENES[20], GENES[8], GENES[24]],
    "noise3": [GENES[0], GENES[15], GENES[18], GENES[5], GENES[21]],
    "noise4": [GENES[3], GENES[9], GENES[16], GENES[23], GENES[12]],
    "noise5": [GENES[7], GENES[14], GENES[19], GENES[1], GENES[10]],
    "noise6": [GENES[5], GENES[11], GENES[22], GENES[2], GENES[17]],
}


def _gene_sets_df(specs: dict[str, list[str]], names: dict[str, str] | None = None) -> pd.DataFrame:
    names = names or {}
    rows = [
        {"set_id": set_id, "set_name": names.get(set_id, set_id), "source_db": "reactome", "gene_id": gene_id}
        for set_id, gene_ids in specs.items()
        for gene_id in gene_ids
    ]
    return pd.DataFrame(rows)


def _all_gene_sets() -> pd.DataFrame:
    specs = {
        "strong_pathway": STRONG_PATHWAY,
        "borderline_pathway": BORDERLINE_PATHWAY,
        "dup_ancestor": DUP_ANCESTOR,
        "dup_child": DUP_CHILD,
        **NOISE_SETS,
    }
    names = {
        "strong_pathway": "Strongly Enriched Pathway",
        "borderline_pathway": "Borderline Pathway",
        "dup_ancestor": "Ancestor Broad Pathway",
        "dup_child": "Child Specific Pathway",
    }
    return _gene_sets_df(specs, names)


def _signature() -> pd.DataFrame:
    """Build the ranking signature the same way stage2_gsea_discovery.py does in
    production: via build_signature() over an ot_disease_subset-shaped
    DataFrame, not a hand-built (gene, score) table directly."""
    ot_disease_subset = pd.DataFrame(
        [
            {
                "gene_id": gid,
                "disease_id": "EFO_TEST",
                "datatype_id": "genetic_association",
                "datasource_id": "eva",
                "score": _SCORES[gid],
            }
            for gid in GENES
        ]
    )
    return gd.build_signature(ot_disease_subset)


@pytest.fixture(scope="module")
def gsea_result() -> pd.DataFrame:
    """Collapse + run GSEA once for the shared fixture; every threshold
    assertion below reads from this same run."""
    collapsed = gd.collapse_near_duplicate_sets(_all_gene_sets())
    result, _sets_by_id = gd.run_gsea(_signature(), collapsed, seed=0)
    return result


@pytest.fixture(scope="module")
def significant(gsea_result: pd.DataFrame) -> pd.DataFrame:
    return gd.filter_results(gsea_result, PVAL_THRESHOLD, FDR_THRESHOLD)


# ==========================================================================
# 1. The deliberately-enriched pathway is returned with p-value below threshold.
# ==========================================================================


def test_deliberately_enriched_pathway_is_significant(gsea_result, significant):
    strong = gsea_result.set_index("set_id").loc["strong_pathway"]
    assert strong.pval < PVAL_THRESHOLD
    assert strong.fdr < FDR_THRESHOLD
    assert "strong_pathway" in set(significant.set_id)


# ==========================================================================
# 2. BH-FDR correction is actually applied: an individually-significant
#    pathway becomes non-significant once corrected across the tested library.
# ==========================================================================


def test_fdr_correction_makes_borderline_pathway_nonsignificant(gsea_result, significant):
    borderline = gsea_result.set_index("set_id").loc["borderline_pathway"]
    # Individually significant by the raw p-value alone...
    assert borderline.pval < PVAL_THRESHOLD
    # ...but BH-FDR, corrected across all 8 pathways actually tested, pushes
    # it back out -- this is only possible because run_gsea tests it
    # alongside the rest of the library rather than in isolation.
    assert borderline.fdr > FDR_THRESHOLD
    assert "borderline_pathway" not in set(significant.set_id)
    # The positive control from test 1 must still survive: the fixture
    # isn't just failing everything.
    assert "strong_pathway" in set(significant.set_id)


# ==========================================================================
# 3. Pathways not meeting threshold are excluded from output.
# ==========================================================================


def test_noise_pathways_are_tested_but_excluded(gsea_result, significant):
    noise_ids = set(NOISE_SETS)
    assert noise_ids <= set(gsea_result.set_id), "noise pathways should still be tested"
    assert not (noise_ids & set(significant.set_id)), "noise pathways should never pass the significance filter"


def test_filter_results_only_returns_rows_meeting_both_thresholds(significant):
    assert (significant["pval"] <= PVAL_THRESHOLD).all()
    assert (significant["fdr"] <= FDR_THRESHOLD).all()


# ==========================================================================
# 4. Two near-identical gene sets (Jaccard > 0.7) are collapsed to one before
#    testing.
# ==========================================================================


def test_collapse_near_duplicate_sets_jaccard_math():
    """A minimal, hand-verifiable Jaccard case, independent of the shared
    fixture: "big" (10 genes) and "small" (its first 8 genes) have
    Jaccard = 8/10 = 0.8 > 0.7, so "big" must be dropped and the smaller
    "small" kept; "unrelated" shares nothing and must survive untouched."""
    gene_sets = _gene_sets_df(
        {
            "big": [f"X{i}" for i in range(10)],
            "small": [f"X{i}" for i in range(8)],
            "unrelated": [f"Y{i}" for i in range(5)],
        }
    )
    collapsed = gd.collapse_near_duplicate_sets(gene_sets)
    assert set(collapsed.set_id) == {"small", "unrelated"}


def test_collapse_drops_ancestor_keeps_child_in_full_fixture():
    collapsed = gd.collapse_near_duplicate_sets(_all_gene_sets())
    kept = set(collapsed.set_id)
    assert "dup_child" in kept
    assert "dup_ancestor" not in kept
    # Nothing else was collapsed as collateral damage.
    assert kept == set(_all_gene_sets().set_id) - {"dup_ancestor"}


def test_collapsed_ancestor_never_reaches_gsea(gsea_result):
    """Confirms the collapse happens *before* testing, not just before
    output filtering: the dropped ancestor must be absent even from
    gsea_result (every tested pathway, unfiltered by significance)."""
    assert "dup_ancestor" not in set(gsea_result.set_id)
    assert "dup_child" in set(gsea_result.set_id)


# ==========================================================================
# Bonus: build_signature()'s genetic_association-only restriction, which the
# module docstring calls out as load-bearing for avoiding circularity.
# ==========================================================================


def test_build_signature_ignores_non_genetic_association_rows():
    ot_disease_subset = pd.DataFrame(
        [
            {"gene_id": "G1", "disease_id": "EFO_TEST", "datatype_id": "genetic_association",
             "datasource_id": "eva", "score": 0.20},
            # Decoy: a much higher score under a different datatype must be ignored.
            {"gene_id": "G1", "disease_id": "EFO_TEST", "datatype_id": "affected_pathway",
             "datasource_id": "reactome", "score": 0.99},
        ]
    )
    signature = gd.build_signature(ot_disease_subset)
    assert signature.set_index("gene").loc["G1", "score"] == 0.20
