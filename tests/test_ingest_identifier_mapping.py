"""Identifier mapping (README.md's Stage 1 testing requirement):

a clean HGNC symbol maps to the right Ensembl ID; a deprecated symbol maps
via the synonym field; an unmappable symbol is recorded as dropped and
counted correctly in coverage_report.json.

Gene universe (tests/fixtures/opentargets/target.parquet): TP53 (approved
symbol, clean case), PTGS2 with synonym "COX2" (deprecated/alias-symbol
case — COX2 is not PTGS2's approved symbol), and FAKEGENE123 (referenced by
the GMT/FI fixtures but absent from target.parquet — the unmappable case).
"""

import json
from pathlib import Path

from conftest import make_args
from stage1_ingest import _build_symbol_index, build_genes_table, resolve_symbol, run_ingest

FIXTURE_TARGET_PARQUET = Path(__file__).parent / "fixtures" / "opentargets" / "target.parquet"


def test_resolve_symbol_clean_synonym_and_unmappable():
    genes_df = build_genes_table([FIXTURE_TARGET_PARQUET])
    index = _build_symbol_index(genes_df)

    gene_id, reason = resolve_symbol("TP53", index)
    assert gene_id == "ENSG00000141510"
    assert reason == "approved_symbol"

    gene_id, reason = resolve_symbol("COX2", index)
    assert gene_id == "ENSG00000073756"  # PTGS2's gene_id, reached via the synonym field
    assert reason == "synonym"

    gene_id, reason = resolve_symbol("FAKEGENE123", index)
    assert gene_id is None
    assert reason == "unmapped"


def test_unmappable_symbol_recorded_in_coverage_report(tmp_path, patched_acquisition):
    out_dir = tmp_path / "out"
    run_ingest(make_args(out_dir=out_dir, data_dir=tmp_path / "data"))

    coverage = json.loads((out_dir / "coverage_report.json").read_text())
    gmt_entry = coverage["reactome_gmt"]

    assert "FAKEGENE123" in gmt_entry["example_dropped_ids"]
    assert gmt_entry["dropped_count"] == gmt_entry["raw_gene_count"] - gmt_entry["mapped_count"]
    # TP53, COX2/PTGS2, and BRCA1 all map cleanly across the two fixture
    # pathways, so FAKEGENE123 is the only drop.
    assert gmt_entry["dropped_count"] == 1

    fi_entry = coverage["reactome_fi"]
    assert "FAKEGENE123" in fi_entry["example_dropped_ids"]
    assert fi_entry["dropped_count"] >= 1
