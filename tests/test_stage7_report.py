"""pytest tests for stage7_report.py (Stage 7: HTML report).

`evidence_trace` is tested as a pure function on
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
    """T is the only seed gene. Disease-relevant pathways: P1 = {T, C} (contains a seed),
    P2 = {C, X} (no seed). P3 = {T, C} contains a seed but is not disease-relevant, so it was
    not used to build the graph and is not shown. Graph nodes: T, C, X.
    Interaction C->T (in graph) and C->Y (Y not in graph -> dropped)."""
    T, C, X, Y = (f"ENSG{n:011d}" for n in (1, 2, 3, 4))
    gene_sets = pd.DataFrame(
        [
            {"set_id": "P1", "set_name": "Seed Pathway", "gene_id": T},
            {"set_id": "P1", "set_name": "Seed Pathway", "gene_id": C},
            {"set_id": "P2", "set_name": "Disease Pathway", "gene_id": C},
            {"set_id": "P2", "set_name": "Disease Pathway", "gene_id": X},
            {"set_id": "P3", "set_name": "Other Seed Pathway", "gene_id": T},
            {"set_id": "P3", "set_name": "Other Seed Pathway", "gene_id": C},
        ]
    )
    return {
        "manifest": {},
        "seeds": {T},
        "symbols": {T: "TGT", C: "CAND", X: "XG", Y: "YG"},
        "graph_index": {T, C, X},
        "gene_sets": gene_sets,
        "disease_pathways": pd.DataFrame([{"set_id": "P1"}, {"set_id": "P2"}]),
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

    assert tr["shared_pathways"] == [("Seed Pathway", 1)]     # P1 has C and one seed (T); P3 is not disease-relevant
    assert tr["member_pathways"] == ["Disease Pathway"]       # P2 has C but no seed

    # datatypes: max over datasources, every datatype shown, strongest first
    assert [d["datatype"] for d in tr["datatypes"]] == ["affected_pathway", "genetic_association"]
    ga = tr["datatypes"][1]
    assert ga["score"] == pytest.approx(0.9)
    assert set(ga["datasources"]) == {"eva", "gwas"}
    assert "feeds" not in ga  # the report shows no genetic-evidence score, so nothing "feeds" it

    # only the interaction whose other endpoint is a graph node survives
    assert len(tr["interactions"]) == 1
    it = tr["interactions"][0]
    assert (it["a"], it["b"]) == ("CAND", "TGT")
    assert it["to_seed"] is True
    assert it["sign"] == "activating (+)"


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
    assert tbody.count("<tr") == len(annotated)          # one table row per candidate
    # genetic evidence is not in the report, although candidates_annotated.tsv has the column
    assert "genetic_evidence_score" in annotated.columns
    thead = html.split("<thead>", 1)[1].split("</thead>", 1)[0]
    assert "genetic ev." not in thead and "genetic-evidence score" not in html
    assert html.count('class="trace-entry"') == len(annotated)  # one trace per candidate
    # The traces live in the side panel, which is closed until a symbol is clicked, not below the table.
    panel = html.split('<aside id="tracepanel"', 1)[1].split("</aside>", 1)[0]
    assert panel.startswith(" hidden")
    assert panel.count('class="trace-entry"') == len(annotated)
    assert html.index("</table>") < html.index('<aside id="tracepanel"')
    assert "function openGene" in html and "function closeGene" in html

    for gid in annotated["gene"]:
        sym = symbols[gid]
        assert f'id="g-{gid}"' in html                    # anchored detail block
        # the ranked table shows the SYMBOL (as a link), with the Ensembl id in its own cell
        assert f">{sym}</a>" in tbody
        assert f'class="txt ens">{gid}<' in tbody


@pytest.mark.skipif(
    not (EXAMPLE_RUN / "candidates_annotated.tsv").exists(),
    reason="example_data/example_run not populated",
)
def test_build_report_traces_only_top_candidates(tmp_path):
    html = rep.build_report(EXAMPLE_RUN, tmp_path / "report.html", n_traces=2).read_text()
    annotated = pd.read_csv(EXAMPLE_RUN / "candidates_annotated.tsv", sep="\t")
    tbody = html.split("<tbody", 1)[1].split("</tbody>", 1)[0]
    assert tbody.count("<tr") == len(annotated)          # the table still lists every candidate
    assert html.count('class="trace-entry"') == 2
    for rank, gid in enumerate(annotated["gene"], start=1):
        has_trace = f'id="g-{gid}"' in html
        assert has_trace == (rank <= 2)
        assert (f"openGene('g-{gid}')" in tbody) == has_trace  # no link to a trace that is not there
        assert (f'data-trace="g-{gid}"' in tbody) == has_trace  # the row highlighted while its trace is open
    assert f"top 2 of {len(annotated)} candidates" in html


# --------------------------------------------------------------------------
# graceful degradation: missing optional artifacts
# --------------------------------------------------------------------------


@pytest.mark.skipif(
    not (EXAMPLE_RUN / "candidates_annotated.tsv").exists(),
    reason="example_data/example_run not populated",
)
def test_partial_run_dir_degrades_gracefully(tmp_path):
    """Only manifest.json + candidates_annotated.tsv present, and the table
    has bfwalk_score but no rwr_score or genetic_evidence_score (a default,
    bfwalk-only Stage 5 on a Stage 3 graph). The report still renders: the
    missing columns are omitted, every trace section
    says which artifact is missing, and a header note lists them."""
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "manifest.json").write_text((EXAMPLE_RUN / "manifest.json").read_text())
    ann = pd.read_csv(EXAMPLE_RUN / "candidates_annotated.tsv", sep="\t").drop(
        columns=["rwr_score", "genetic_evidence_score"], errors="ignore"
    )
    assert "bfwalk_score" in ann.columns
    ann.to_csv(run_dir / "candidates_annotated.tsv", sep="\t", index=False)

    html = rep.build_report(run_dir, run_dir / "report.html").read_text()

    thead = html.split("<thead>", 1)[1].split("</thead>", 1)[0]
    assert ">rwr<" not in thead and "genetic ev." not in thead   # optional columns dropped
    assert ">composite<" not in thead and ">bfwalk<" in thead    # the score that ran is shown
    assert html.count('class="trace-entry"') == len(ann)  # every candidate still rendered

    assert "Partial run:" in html
    for artifact in ("gene_sets.parquet", "interactions.parquet", "ot_disease_subset.parquet"):
        assert f"{artifact} not in the run directory" in html

    assert "http://" not in html and "https://" not in html      # still self-contained


def test_load_run_requires_only_manifest_and_annotated(tmp_path):
    """A directory with neither required artifact fails loudly; with both,
    it loads (every other artifact optional)."""
    with pytest.raises(SystemExit):
        rep.load_run(tmp_path)

    (tmp_path / "manifest.json").write_text((EXAMPLE_RUN / "manifest.json").read_text())
    (tmp_path / "candidates_annotated.tsv").write_text(
        (EXAMPLE_RUN / "candidates_annotated.tsv").read_text()
    )
    run = rep.load_run(tmp_path)
    assert run["gene_sets"] is None and run["graph_index"] is None
    assert "gene_sets.parquet" in run["missing"]
