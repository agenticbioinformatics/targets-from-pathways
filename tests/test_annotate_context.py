"""pytest tests for stage6_annotate_context.py (Stage 6: contextual annotation).

Hand-built inputs, calling the module's functions directly, plus full
``run()`` round-trips through a written run directory (candidate_scores.tsv +
graph_metadata.json + manifest.json + a tiny Open Targets `target` parquet).

Stage 6 only annotates: genetic evidence, tractability and safety are added
for reading, and candidates keep Stage 5's ranking (bfwalk_score, else
rwr_score). There is no composite score.
"""

from __future__ import annotations

import json

import pandas as pd
import pandera.errors
import pytest

import stage6_annotate_context as ac
from stage0_schemas import AnnotatedCandidatesSchema, CandidateScoresSchema

A = "ENSG00000000001"
B = "ENSG00000000002"
CAND = "ENSG00000000003"  # a third gene, absent from the target parquet fixture


# --------------------------------------------------------------------------
# tractability / safety derivation
# --------------------------------------------------------------------------


def testtractability_bucket():
    assert ac.tractability_bucket([{"id": "Approved Drug", "modality": "SM", "value": True}]) == "clinical"
    assert ac.tractability_bucket(
        [{"id": "Phase 1 Clinical", "modality": "SM", "value": True}]
    ) == "clinical"
    assert ac.tractability_bucket(
        [{"id": "Structure with Ligand", "modality": "SM", "value": True}]
    ) == "discovery"
    # a flag that is present but False does not count
    assert ac.tractability_bucket([{"id": "Druggable Family", "modality": "SM", "value": False}]) == "unknown"
    assert ac.tractability_bucket([]) == "unknown"
    assert ac.tractability_bucket(None) == "unknown"


def test_safety_flag_and_count():
    assert ac.safety_flag([{"event": "x"}, {"event": "y"}]) == ("has_liabilities", 2)
    # no annotation at all -> unknown, never a numeric/"safe" default
    assert ac.safety_flag([]) == ("unknown", 0)
    assert ac.safety_flag(None) == ("unknown", 0)


# --------------------------------------------------------------------------
# schemas
# --------------------------------------------------------------------------


def test_scores_without_any_propagation_method_are_rejected():
    # Stage 5 always runs rwr and/or bfwalk, so a scores table with neither is malformed.
    with pytest.raises(pandera.errors.SchemaError, match="rwr_score"):
        CandidateScoresSchema.validate(pd.DataFrame({"gene": [A]}))
    annotated = pd.DataFrame({
        "gene": [A], "tractability": ["unknown"], "safety": ["unknown"], "n_safety_liabilities": [0],
    })
    with pytest.raises(pandera.errors.SchemaError, match="rwr_score"):
        AnnotatedCandidatesSchema.validate(annotated)
    AnnotatedCandidatesSchema.validate(annotated.assign(rwr_score=0.1))


# --------------------------------------------------------------------------
# build_annotations: absent gene -> unknown/unknown, never penalised
# --------------------------------------------------------------------------


def _write_target_parquet(path):
    df = pd.DataFrame(
        {
            "id": [A, B],
            "tractability": [
                [{"modality": "SM", "id": "Approved Drug", "value": True}],
                [{"modality": "SM", "id": "Structure with Ligand", "value": True}],
            ],
            # A: two liabilities; B: present but NO safety annotation (empty list)
            "safetyLiabilities": [
                [{"event": "hepatotoxicity", "eventId": "EFO_0000001"}],
                [],
            ],
        }
    )
    df.to_parquet(path, index=False)


def test_build_annotations_absent_gene_is_unknown_not_penalised(tmp_path):
    pq = tmp_path / "target.parquet"
    _write_target_parquet(pq)

    ann = ac.build_annotations([pq], [A, B, CAND]).set_index("gene")

    assert ann.loc[A, "tractability"] == "clinical"
    assert ann.loc[A, "safety"] == "has_liabilities"
    assert ann.loc[A, "n_safety_liabilities"] == 1

    # B: present in the parquet but empty safetyLiabilities -> unknown, not "safe"
    assert ann.loc[B, "safety"] == "unknown"
    assert ann.loc[B, "n_safety_liabilities"] == 0

    # CAND: absent from the parquet entirely -> unknown / unknown / 0, no numeric default
    assert ann.loc[CAND, "tractability"] == "unknown"
    assert ann.loc[CAND, "safety"] == "unknown"
    assert ann.loc[CAND, "n_safety_liabilities"] == 0


# --------------------------------------------------------------------------
# end-to-end run()
# --------------------------------------------------------------------------


def _write_run_dir(path, *, with_bfwalk=True, with_stage4=True):
    # Deliberately not in ranking order, so run() must restore Stage 5's ranking itself.
    cols = {"gene": [A, B, CAND], "rwr_score": [0.4, 0.2, 0.1]}
    if with_bfwalk:
        cols["bfwalk_score"] = [1.0, 3.0, 0.5]
    pd.DataFrame(cols).to_csv(path / "candidate_scores.tsv", sep="\t", index=False)

    meta = {"seed": 0}
    if with_stage4:
        meta["genetic_evidence_score"] = {A: 0.8, B: 0.5, CAND: 0.1}
    (path / "graph_metadata.json").write_text(json.dumps(meta))

    _write_target_parquet(path / "ot_target_subset.parquet")
    manifest = {
        "run_id": "r", "git_commit": "u", "created_at": "2026-08-28T00:00:00Z", "seed": 0,
        "disease": {"efo_id": "EFO_0000001", "name": "d", "n_associated_genes": 3},
        "sources": [
            {"db": "opentargets", "version": "26.06", "files": [
                {"path": str(path / "ot_target_subset.parquet"), "sha256": "a" * 64, "bytes": 1}
            ]}
        ],
        "cli_parameters": {}, "output_artifacts": [],
        "coverage_report": {}, "scale_report": {
            "gene_set_size_distribution_before_cap": [], "gene_set_size_distribution_after_cap": [],
            "sets_retained": 0, "interaction_counts_by_sign": {},
            "interaction_counts_by_evidence_type": {},
        },
    }
    (path / "manifest.json").write_text(json.dumps(manifest))


def test_run_end_to_end_columns_ranking_and_unknown_gene(tmp_path):
    _write_run_dir(tmp_path)
    out = ac.run(ac.parse_args(["--manifest", str(tmp_path / "manifest.json")]))
    df = pd.read_csv(out, sep="\t")

    assert list(df.columns) == [
        "gene", "rwr_score", "bfwalk_score", "genetic_evidence_score",
        "tractability", "safety", "n_safety_liabilities",
    ]
    # Stage 5's ranking (bfwalk over rwr), whatever order the candidates file is in
    assert list(df["gene"]) == [B, A, CAND]
    assert dict(zip(df["gene"], df["genetic_evidence_score"])) == pytest.approx({A: 0.8, B: 0.5, CAND: 0.1})

    row = df.set_index("gene").loc[CAND]
    assert row["tractability"] == "unknown" and row["safety"] == "unknown"
    assert row["n_safety_liabilities"] == 0


def test_run_rwr_only_without_stage4_ranks_by_rwr(tmp_path):
    _write_run_dir(tmp_path, with_bfwalk=False, with_stage4=False)
    out = ac.run(ac.parse_args(["--manifest", str(tmp_path / "manifest.json")]))
    df = pd.read_csv(out, sep="\t")

    assert list(df.columns) == ["gene", "rwr_score", "tractability", "safety", "n_safety_liabilities"]
    assert list(df["gene"]) == [A, B, CAND]


def test_annotations_never_change_the_ranking(tmp_path):
    """A has a clinical tractability flag and strong genetic evidence, CAND has neither; ranked by
    bfwalk_score they stay in Stage 5's order (B, A, CAND) - no composite reorders them."""
    _write_run_dir(tmp_path)
    df = pd.read_csv(ac.run(ac.parse_args(["--manifest", str(tmp_path / "manifest.json")])), sep="\t")
    assert df["bfwalk_score"].is_monotonic_decreasing
    assert df.set_index("gene").loc[A, "tractability"] == "clinical"
    assert list(df["gene"]).index(A) > list(df["gene"]).index(B)
