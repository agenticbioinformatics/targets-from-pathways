"""Schema conformance (README.md's Stage 1 testing requirement):

every emitted artifact validates against its Stage 0 schema, and a source
file with a corrupt row fails loudly rather than being silently skipped.
"""

from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

import stage1_ingest as ingest
from conftest import make_args
from stage1_ingest import run_ingest
from stage0_schemas import GeneSetsSchema, GenesSchema, InteractionsSchema, OTAssociationsSchema, assert_foreign_key

FIXTURES_DIR = Path(__file__).parent / "fixtures"
TARGET_FIXTURE = FIXTURES_DIR / "opentargets" / "target.parquet"


def test_every_emitted_artifact_validates_against_its_stage0_schema(tmp_path, patched_acquisition):
    out_dir = tmp_path / "out"
    run_ingest(make_args(out_dir=out_dir, data_dir=tmp_path / "data"))

    # Deliberately a bare pd.read_parquet, not the in-process DataFrame
    # run_ingest already validated before writing — this is "validate on
    # read" for whatever the next stage would actually do.
    genes = pd.read_parquet(out_dir / "genes.parquet")
    gene_sets = pd.read_parquet(out_dir / "gene_sets.parquet")
    interactions = pd.read_parquet(out_dir / "interactions.parquet")
    ot_subset = pd.read_parquet(out_dir / "ot_disease_subset.parquet")

    GenesSchema.validate(genes)
    GeneSetsSchema.validate(gene_sets)
    InteractionsSchema.validate(interactions)
    OTAssociationsSchema.validate(ot_subset)

    assert_foreign_key(gene_sets, "gene_id", genes)
    assert_foreign_key(interactions, "gene_a", genes)
    assert_foreign_key(interactions, "gene_b", genes)
    assert_foreign_key(ot_subset, "gene_id", genes)

    # Non-empty, so the above validations exercise real rows, not just an
    # empty-frame pass-through.
    assert len(genes) == 3
    assert len(gene_sets) > 0
    assert len(interactions) > 0
    assert len(ot_subset) > 0


def test_malformed_fi_row_fails_loudly(tmp_path, monkeypatch):
    """A ragged row (fewer fields than the header) in the FI file must raise,
    not get silently dropped.

    This is a real, non-obvious trap: pandas' C parser does NOT raise on a
    short row — it pads the missing trailing fields with NaN, which would
    otherwise miss every _FI_DIRECTION_CODES key and vanish with nothing but
    a logged warning (verified: a naive `pd.read_csv` over this exact file
    silently returns a NaN-padded row instead of raising — caught during
    review of this very test, before stage1_ingest.py had an explicit check for it).
    The zip wrapping matters here too: it's what makes this test actually
    reach the row-parsing code, instead of failing earlier and for an
    unrelated reason (an un-zipped fixture makes _extract_single_member's
    zipfile.ZipFile() call raise BadZipFile before any row is ever read).
    """
    import zipfile

    bad_fi_txt = tmp_path / "bad_fi_src.txt"
    bad_fi_txt.write_text(
        "Gene1\tGene2\tAnnotation\tDirection\tScore\n"
        "TP53\tCOX2\tactivate\t->\t1.00\n"
        "TP53\tBRCA1\tcomplex\n"  # missing Direction and Score columns
    )
    bad_fi_zip = tmp_path / "bad_fi.zip"
    with zipfile.ZipFile(bad_fi_zip, "w") as zf:
        zf.write(bad_fi_txt, arcname="FIsInGene_test_with_annotations.txt")

    from conftest import _make_fake_ensure_directory_source, _make_fake_ensure_pinned_file

    monkeypatch.setattr(
        ingest, "_ensure_pinned_file", _make_fake_ensure_pinned_file(tmp_path, overrides={"reactome_fi": bad_fi_zip})
    )
    monkeypatch.setattr(ingest, "_ensure_directory_source", _make_fake_ensure_directory_source())

    out_dir = tmp_path / "out"
    with pytest.raises(SystemExit):
        run_ingest(make_args(out_dir=out_dir, data_dir=tmp_path / "data"))
    assert not (out_dir / "manifest.json").exists()


def test_unrecognized_datasource_fails_loudly(tmp_path, monkeypatch, patched_acquisition):
    """An OT association row with a datasource this module doesn't know how
    to bucket into a datatype must hard-fail (see stage1_ingest.py's
    _DATASOURCE_TO_DATATYPE), not be silently dropped or mis-tagged."""
    bad_assoc = tmp_path / "bad_association.parquet"
    schema = pa.schema(
        [
            ("targetId", pa.string()),
            ("diseaseId", pa.string()),
            ("aggregationType", pa.string()),
            ("aggregationValue", pa.string()),
            ("associationScore", pa.float64()),
        ]
    )
    rows = [
        {
            "targetId": "ENSG00000141510",
            "diseaseId": "EFO_9000001",
            "aggregationType": "datasourceId",
            "aggregationValue": "totally_bogus_datasource",
            "associationScore": 0.5,
        }
    ]
    pq.write_table(pa.Table.from_pylist(rows, schema=schema), bad_assoc)

    from conftest import _make_fake_ensure_directory_source

    monkeypatch.setattr(
        ingest,
        "_ensure_directory_source",
        _make_fake_ensure_directory_source(
            overrides={"opentargets_association_by_datasource_direct": [bad_assoc]}
        ),
    )

    out_dir = tmp_path / "out"
    with pytest.raises(SystemExit):
        run_ingest(make_args(out_dir=out_dir, data_dir=tmp_path / "data"))
    assert not (out_dir / "manifest.json").exists()
