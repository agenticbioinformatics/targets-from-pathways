"""pytest tests for stage7_report.py (Stage 7: HTML report).

`evidence_trace` and `benchmark_summary` are tested as pure functions on
hand-built inputs; `build_report` is exercised end to end against the
checked-in example run (read-only — the HTML is written to tmp_path).
"""

from __future__ import annotations

import json

import pandas as pd
import pytest

import stage7_report as rep

EXAMPLE_RUN = rep.REPO_DIR / "example_data" / "example_run"


# --------------------------------------------------------------------------
# evidence_trace
# --------------------------------------------------------------------------


def _mini_run():
    """T is the target. Pathways: P1 = {T, C} (target-containing),
    P2 = {C, X} (disease-relevant, no target). Graph nodes: T, C, X.
    Interaction C->T (in graph) and C->Y (Y not in graph -> dropped)."""
    T, C, X, Y = (f"ENSG{n:011d}" for n in (1, 2, 3, 4))
    gene_sets = pd.DataFrame(
        [
            {"set_id": "P1", "set_name": "Target Pathway", "gene_id": T},
            {"set_id": "P1", "set_name": "Target Pathway", "gene_id": C},
            {"set_id": "P2", "set_name": "Disease Pathway", "gene_id": C},
            {"set_id": "P2", "set_name": "Disease Pathway", "gene_id": X},
        ]
    )
    return {
        "manifest": {"resolved_target": {"gene_id": T}},
        "symbols": {T: "TGT", C: "CAND", X: "XG", Y: "YG"},
        "graph_index": {T, C, X},
        "admitted_datatypes": {"genetic_association", "known_drug"},
        "stage4_ran": True,
        "gene_sets": gene_sets,
        "disease_pathways": pd.DataFrame([{"set_id": "P2"}]),
        "ot": pd.DataFrame(
            [
                {"gene_id": C, "datatype_id": "genetic_association", "datasource_id": "eva", "score": 0.7},
                {"gene_id": C, "datatype_id": "genetic_association", "datasource_id": "gwas", "score": 0.9},
                {"gene_id": C, "datatype_id": "affected_pathway", "datasource_id": "reactome", "score": 0.99},
            ]
        ),
        "interactions": pd.DataFrame(
            [
                {"gene_a": C, "gene_b": T, "directed": True, "sign": 1,
                 "evidence_type": "curated", "confidence": 0.8},
                {"gene_a": C, "gene_b": Y, "directed": True, "sign": -1,
                 "evidence_type": "curated", "confidence": 0.5},
            ]
        ),
    }, (T, C, X, Y)


def test_evidence_trace_shared_pathways_and_scoping():
    run, (T, C, X, Y) = _mini_run()
    tr = rep.evidence_trace(C, run)

    assert tr["shared_pathways"] == ["Target Pathway"]        # P1 has both C and T
    assert tr["member_pathways"] == ["Disease Pathway"]       # P2 has C, not T

    # datatypes: max over datasources, admitted flagged, disallowed still shown
    ga = next(d for d in tr["datatypes"] if d["datatype"] == "genetic_association")
    assert ga["score"] == pytest.approx(0.9)
    assert ga["feeds"] is True
    assert set(ga["datasources"]) == {"eva", "gwas"}
    ap = next(d for d in tr["datatypes"] if d["datatype"] == "affected_pathway")
    assert ap["feeds"] is False

    # only the interaction whose other endpoint is a graph node survives
    assert len(tr["interactions"]) == 1
    it = tr["interactions"][0]
    assert (it["a"], it["b"]) == ("CAND", "TGT")
    assert it["to_target"] is True
    assert it["sign"] == "activating (+)"


def test_evidence_trace_feeds_false_when_stage4_absent():
    run, (T, C, X, Y) = _mini_run()
    run["stage4_ran"] = False
    tr = rep.evidence_trace(C, run)
    assert all(d["feeds"] is False for d in tr["datatypes"])


# --------------------------------------------------------------------------
# benchmark_summary fallback
# --------------------------------------------------------------------------


def test_benchmark_summary_prefers_real_then_example_then_none(tmp_path):
    text, label = rep.benchmark_summary(tmp_path)
    assert text is None and label == "none"

    (tmp_path / "benchmark_summary.example.txt").write_text("example numbers")
    text, label = rep.benchmark_summary(tmp_path)
    assert text == "example numbers" and "EXAMPLE" in label

    (tmp_path / "benchmark_summary.txt").write_text("real numbers")
    text, label = rep.benchmark_summary(tmp_path)
    assert text == "real numbers" and label == "benchmark_summary.txt"


# --------------------------------------------------------------------------
# build_report end to end (against the checked-in example run)
# --------------------------------------------------------------------------


@pytest.mark.skipif(
    not (EXAMPLE_RUN / "candidates_annotated.tsv").exists(),
    reason="example_data/example_run not populated — run example_data/build_gsea_example.py + Stages 3-6",
)
def test_build_report_is_self_contained_and_symbol_labelled(tmp_path):
    out = rep.build_report(EXAMPLE_RUN, tmp_path / "report.html")
    html = out.read_text()

    # single self-contained file: no external resources
    assert "http://" not in html and "https://" not in html
    assert "<script>" in html and "<style>" in html

    annotated = pd.read_csv(EXAMPLE_RUN / "candidates_annotated.tsv", sep="\t")
    genes = pd.read_parquet(EXAMPLE_RUN / "genes.parquet")
    symbols = dict(zip(genes["gene_id"], genes["symbol"]))
    assert set(symbols.values()) != set(symbols.keys())  # symbols really differ from Ensembl ids

    tbody = html.split("<tbody", 1)[1].split("</tbody>", 1)[0]
    assert tbody.count("<tr>") == len(annotated)          # one table row per candidate
    assert html.count("<details ") == len(annotated)      # one detail block per candidate

    for gid in annotated["gene"]:
        sym = symbols[gid]
        assert f'id="g-{gid}"' in html                    # anchored detail block
        # the ranked table shows the SYMBOL (as a link), with the Ensembl id in its own cell
        assert f">{sym}</a>" in tbody
        assert f'class="txt ens">{gid}<' in tbody

    # benchmark panel embedded (the example placeholder, clearly labelled)
    assert "Benchmark validation" in html
    assert "EXAMPLE" in html
