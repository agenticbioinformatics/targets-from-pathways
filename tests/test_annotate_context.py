"""pytest tests for annotate_context.py (Stage 6: contextual annotation +
composite score).

Hand-built inputs, calling the module's functions directly, plus one full
``run()`` round-trip through a written run directory (candidate_scores.tsv +
graph_metadata.json + manifest.json + a tiny Open Targets `target` parquet).

Composite fixture (two candidates, so min-max maps them to the endpoints):

    gene  topology_score  rwr_score  genetic_evidence_score  tractability  safety
    A     1.0             0.4        0.8                     clinical      has_liabilities
    B     3.0             0.2        0.5                     unknown       unknown

Normalised components: topology = min-max(topology_score) -> A=0, B=1;
rwr = min-max(rwr_score) -> A=1, B=0; genetic_evidence as-is -> A=0.8,
B=0.5; tractability {clinical:1, unknown:0}; safety {has_liabilities:0,
unknown:0.5}.

With weights topology:0.4, rwr:0.3, genetic_evidence:0.2, tractability:0.05,
safety:0.05 (sum 1):
    composite A = .3*1 + .2*.8 + .05*1            = 0.510
    composite B = .4*1 + .2*.5 + .05*.5           = 0.525
so B edges out A.
"""

from __future__ import annotations

import json

import pandas as pd
import pytest

import annotate_context as ac

A = "ENSG00000000001"
B = "ENSG00000000002"
CAND = "ENSG00000000003"  # a third gene, absent from the target parquet fixture


# --------------------------------------------------------------------------
# tractability / safety derivation
# --------------------------------------------------------------------------


def test_bucket_tractability():
    assert ac._bucket_tractability([{"id": "Approved Drug", "modality": "SM", "value": True}]) == "clinical"
    assert ac._bucket_tractability(
        [{"id": "Phase 1 Clinical", "modality": "SM", "value": True}]
    ) == "clinical"
    assert ac._bucket_tractability(
        [{"id": "Structure with Ligand", "modality": "SM", "value": True}]
    ) == "discovery"
    # a flag that is present but False does not count
    assert ac._bucket_tractability([{"id": "Druggable Family", "modality": "SM", "value": False}]) == "unknown"
    assert ac._bucket_tractability([]) == "unknown"
    assert ac._bucket_tractability(None) == "unknown"


def test_safety_flag_and_count():
    assert ac._safety([{"event": "x"}, {"event": "y"}]) == ("has_liabilities", 2)
    # no annotation at all -> unknown, never a numeric/"safe" default
    assert ac._safety([]) == ("unknown", 0)
    assert ac._safety(None) == ("unknown", 0)


# --------------------------------------------------------------------------
# weights
# --------------------------------------------------------------------------


def test_parse_weights_ok_and_rejects_junk():
    w = ac.parse_weights("topology:0.5,safety:0.5")
    assert w == {"topology": 0.5, "safety": 0.5}
    with pytest.raises(SystemExit):
        ac.parse_weights("bogus:0.5")
    with pytest.raises(SystemExit):
        ac.parse_weights("topology:-1")
    with pytest.raises(SystemExit):
        ac.parse_weights("topology:0")  # nothing positive


# --------------------------------------------------------------------------
# composite score
# --------------------------------------------------------------------------


def _composite_fixture() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene": [A, B],
            "topology_score": [1.0, 3.0],
            "rwr_score": [0.4, 0.2],
            "genetic_evidence_score": [0.8, 0.5],
            "tractability": ["clinical", "unknown"],
            "safety": ["has_liabilities", "unknown"],
        }
    )


def test_composite_matches_hand_computed():
    df = _composite_fixture()
    weights = {"topology": 0.4, "rwr": 0.3, "genetic_evidence": 0.2, "tractability": 0.05, "safety": 0.05}
    composite, breakdown, fractions = ac.compute_composite(df, weights)

    # weights already sum to 1, so the renormalised fractions equal them
    assert fractions == pytest.approx(weights)
    assert composite.tolist() == pytest.approx([0.510, 0.525])
    # breakdown is each component's *weighted contribution*; the terms sum to composite
    assert breakdown.iloc[0] == "topology=0.000|rwr=0.300|genetic_evidence=0.160|tractability=0.050|safety=0.000"
    contrib_sum = sum(float(kv.split("=")[1]) for kv in breakdown.iloc[0].split("|"))
    assert contrib_sum == pytest.approx(0.510)


def test_composite_drops_unavailable_components_and_renormalises():
    df = _composite_fixture().drop(columns=["rwr_score", "genetic_evidence_score"])
    # rwr + genetic_evidence are weighted but have no column -> dropped, divisor = 0.2 + 0.3
    composite, breakdown, used = ac.compute_composite(
        df, {"topology": 0.2, "rwr": 0.3, "genetic_evidence": 0.4, "safety": 0.3}
    )
    assert set(used) == {"topology", "safety"}
    # A: (.2*0 + .3*0)/.5 = 0 ; B: (.2*1 + .3*.5)/.5 = 0.7
    assert composite.tolist() == pytest.approx([0.0, 0.7])
    assert "rwr" not in breakdown.iloc[0] and "genetic_evidence" not in breakdown.iloc[0]


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


def _write_run_dir(path, *, with_rwr=True, with_stage4=True):
    cols = {"gene": [A, B, CAND], "topology_score": [1.0, 3.0, 0.5]}
    if with_rwr:
        cols["rwr_score"] = [0.4, 0.2, 0.1]
    pd.DataFrame(cols).to_csv(path / "candidate_scores.tsv", sep="\t", index=False)

    meta = {"seed": 0, "resolved_target": {"gene_id": "ENSG00000000099", "input": "t", "symbol": "t"}}
    if with_stage4:
        meta["genetic_evidence_score"] = {A: 0.8, B: 0.5, CAND: 0.1}
    (path / "graph_metadata.json").write_text(json.dumps(meta))

    _write_target_parquet(path / "ot_target_subset.parquet")
    manifest = {
        "run_id": "r", "git_commit": "u", "created_at": "2026-08-28T00:00:00Z", "seed": 0,
        "resolved_target": {"input": "t", "gene_id": "ENSG00000000099", "symbol": "t"},
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
            "interaction_counts_by_evidence_type": {}, "projected_comembership_edge_count": 0,
        },
    }
    (path / "manifest.json").write_text(json.dumps(manifest))


def test_run_end_to_end_columns_ranking_and_unknown_gene(tmp_path):
    _write_run_dir(tmp_path)
    out = ac.run(ac.parse_args(["--manifest", str(tmp_path / "manifest.json")]))
    df = pd.read_csv(out, sep="\t")

    assert list(df.columns) == [
        "gene", "topology_score", "rwr_score", "genetic_evidence_score",
        "tractability", "safety", "n_safety_liabilities",
        "composite_score", "composite_breakdown", "composite_weights",
    ]
    # composite_breakdown terms sum to composite_score (per-component contributions)
    for _, r in df.iterrows():
        parts = sum(float(kv.split("=")[1]) for kv in r["composite_breakdown"].split("|"))
        assert parts == pytest.approx(r["composite_score"], abs=2e-3)
    # composite_weights is the same self-contained string on every row
    assert df["composite_weights"].nunique() == 1
    assert df["composite_score"].is_monotonic_decreasing
    assert df["composite_score"].between(0.0, 1.0).all()

    row = df.set_index("gene").loc[CAND]
    assert row["tractability"] == "unknown" and row["safety"] == "unknown"
    assert row["n_safety_liabilities"] == 0


def test_run_without_rwr_or_stage4_drops_those_components(tmp_path):
    _write_run_dir(tmp_path, with_rwr=False, with_stage4=False)
    out = ac.run(ac.parse_args(["--manifest", str(tmp_path / "manifest.json")]))
    df = pd.read_csv(out, sep="\t")

    assert "rwr_score" not in df.columns
    assert "genetic_evidence_score" not in df.columns
    assert "composite_score" in df.columns
    assert df["composite_score"].between(0.0, 1.0).all()
    assert "rwr" not in df["composite_breakdown"].iloc[0]


# --------------------------------------------------------------------------
# Fixture tractability/safety table -> hand-computed composite, end to end
# --------------------------------------------------------------------------
#
# Three candidates E1/E2/E3, chosen so min-max lands them on 0 / 0.5 / 1:
#
#   gene  topology  rwr   genetic_evidence   in target table?
#   E1    1.0       0.3   0.9                yes: clinical flag, 1 liability
#   E2    2.0       0.2   0.6                yes: discovery flag, EMPTY liabilities
#   E3    3.0       0.1   0.3                NO (absent entirely)
#
# Normalised components:
#   topology (min-max)         E1=0.0  E2=0.5  E3=1.0
#   rwr      (min-max)          E1=1.0  E2=0.5  E3=0.0
#   genetic_evidence (as-is)    E1=0.9  E2=0.6  E3=0.3
#   tractability {clin1,disc.5,unk0}  E1=1.0  E2=0.5  E3=0.0
#   safety {has_liab 0, unknown .5}   E1=0.0  E2=0.5  E3=0.5
#
# weights topology:0.4 rwr:0.2 genetic_evidence:0.2 tractability:0.1 safety:0.1 (sum 1):
#   E1 = .2*1 + .2*.9 + .1*1            = 0.48
#   E2 = .4*.5 + .2*.5 + .2*.6 + .1*.5 + .1*.5 = 0.52
#   E3 = .4*1 + .2*.3 + .1*.5           = 0.51
# -> ranking E2, E3, E1.

E1, E2, E3 = "ENSG00000000011", "ENSG00000000012", "ENSG00000000013"
_WEIGHTS = "topology:0.4,rwr:0.2,genetic_evidence:0.2,tractability:0.1,safety:0.1"


def _write_fixture_run_dir(path):
    pd.DataFrame(
        {"gene": [E1, E2, E3], "topology_score": [1.0, 2.0, 3.0], "rwr_score": [0.3, 0.2, 0.1]}
    ).to_csv(path / "candidate_scores.tsv", sep="\t", index=False)

    (path / "graph_metadata.json").write_text(
        json.dumps(
            {
                "seed": 0,
                "resolved_target": {"gene_id": "ENSG00000000099", "input": "t", "symbol": "t"},
                "genetic_evidence_score": {E1: 0.9, E2: 0.6, E3: 0.3},
            }
        )
    )

    # small tractability/safety table: E1 and E2 only; E2 has no liabilities
    # annotation at all; E3 is deliberately absent.
    pd.DataFrame(
        {
            "id": [E1, E2],
            "tractability": [
                [{"modality": "SM", "id": "Approved Drug", "value": True}],
                [{"modality": "SM", "id": "Structure with Ligand", "value": True}],
            ],
            "safetyLiabilities": [[{"event": "hepatotoxicity", "eventId": "EFO_0000001"}], []],
        }
    ).to_parquet(path / "ot_target_subset.parquet", index=False)

    (path / "manifest.json").write_text(
        json.dumps(
            {
                "run_id": "r", "git_commit": "u", "created_at": "2026-08-28T00:00:00Z", "seed": 0,
                "resolved_target": {"input": "t", "gene_id": "ENSG00000000099", "symbol": "t"},
                "disease": {"efo_id": "EFO_0000001", "name": "d", "n_associated_genes": 3},
                "sources": [
                    {"db": "opentargets", "version": "26.06", "files": [
                        {"path": str(path / "ot_target_subset.parquet"), "sha256": "a" * 64, "bytes": 1}
                    ]}
                ],
                "cli_parameters": {}, "output_artifacts": [], "coverage_report": {},
                "scale_report": {
                    "gene_set_size_distribution_before_cap": [], "gene_set_size_distribution_after_cap": [],
                    "sets_retained": 0, "interaction_counts_by_sign": {},
                    "interaction_counts_by_evidence_type": {}, "projected_comembership_edge_count": 0,
                },
            }
        )
    )


def test_fixture_table_hand_computed_composite_and_unknown_gene(tmp_path):
    _write_fixture_run_dir(tmp_path)
    out = ac.run(ac.parse_args(["--manifest", str(tmp_path / "manifest.json"), "--weights", _WEIGHTS]))
    df = pd.read_csv(out, sep="\t").set_index("gene")

    # composite matches the hand-computation above, exactly
    assert df.loc[E1, "composite_score"] == pytest.approx(0.48)
    assert df.loc[E2, "composite_score"] == pytest.approx(0.52)
    assert df.loc[E3, "composite_score"] == pytest.approx(0.51)
    assert list(df.sort_values("composite_score", ascending=False).index) == [E2, E3, E1]

    # E2 has NO safety annotation (empty list) -> "unknown", never a numeric
    # default and never "safe"
    assert df.loc[E2, "safety"] == "unknown"
    assert df.loc[E2, "n_safety_liabilities"] == 0
    assert df.loc[E2, "tractability"] == "discovery"

    # E3 absent from the table entirely -> unknown / unknown / 0
    assert df.loc[E3, "safety"] == "unknown"
    assert df.loc[E3, "tractability"] == "unknown"
    assert df.loc[E3, "n_safety_liabilities"] == 0

    # E1 does have a liability -> flagged, counted
    assert df.loc[E1, "safety"] == "has_liabilities"
    assert df.loc[E1, "n_safety_liabilities"] == 1
    assert df.loc[E1, "tractability"] == "clinical"
