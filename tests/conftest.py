"""Shared fixtures for stage1_ingest.py's test suite.

The tests never touch the network or the real pinned registry: they patch
stage1_ingest.py's two acquisition entry points (`_ensure_pinned_file`,
`_ensure_directory_source`) to hand back the tiny committed fixtures under
tests/fixtures/ instead. Everything downstream of acquisition — extraction,
identifier mapping, gene_sets/interactions construction, coverage/scale
reporting, schema validation, and manifest writing — runs unmodified.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "pipeline"))
import stage1_ingest as ingest  # noqa: E402

FIXTURES_DIR = Path(__file__).parent / "fixtures"

_PINNED_FIXTURES = {
    "reactome_gmt": FIXTURES_DIR / "reactome" / "ReactomePathways.gmt.zip",
    "reactome_ensembl2reactome": FIXTURES_DIR / "reactome" / "Ensembl2Reactome.txt",
    "reactome_fi": FIXTURES_DIR / "reactome" / "FIsInGene_test_with_annotations.txt.zip",
    "opentargets_disease": FIXTURES_DIR / "opentargets" / "disease.parquet",
}
_DIRECTORY_FIXTURES = {
    "opentargets_target": [FIXTURES_DIR / "opentargets" / "target.parquet"],
    "opentargets_association_by_datasource_indirect": [FIXTURES_DIR / "opentargets" / "association.parquet"],
}


def _make_fake_ensure_pinned_file(tmp_path, overrides: dict | None = None):
    """A fake `_ensure_pinned_file` that serves the committed fixtures (zip
    fixtures copied+kept-zipped into tmp_path, so repeated/parallel test runs
    never write into tests/fixtures/), except for any `key` in `overrides`,
    which is served verbatim — used to swap in a single corrupt source file
    while leaving every other source as the normal tiny fixture."""
    overrides = overrides or {}

    def fake_ensure_pinned_file(key, data_dir, no_download):
        if key in overrides:
            return overrides[key]
        src = _PINNED_FIXTURES[key]
        if src.suffix == ".zip":
            dest = tmp_path / src.name
            dest.write_bytes(src.read_bytes())
            return dest
        return src

    return fake_ensure_pinned_file


def _make_fake_ensure_directory_source(overrides: dict | None = None):
    overrides = overrides or {}

    def fake_ensure_directory_source(key, data_dir, no_download):
        return overrides.get(key, _DIRECTORY_FIXTURES[key])

    return fake_ensure_directory_source


@pytest.fixture
def patched_acquisition(monkeypatch, tmp_path):
    """Redirect stage1_ingest.py's acquisition layer at the tiny fixtures."""
    monkeypatch.setattr(ingest, "_ensure_pinned_file", _make_fake_ensure_pinned_file(tmp_path))
    monkeypatch.setattr(ingest, "_ensure_directory_source", _make_fake_ensure_directory_source())


def make_args(out_dir: Path, data_dir: Path, **overrides) -> argparse.Namespace:
    """Build an stage1_ingest.py CLI Namespace without going through argv parsing."""
    defaults = dict(
        target="TP53",
        disease="EFO_9000001",
        data_dir=data_dir,
        out_dir=out_dir,
        pathway_db="reactome",
        no_download=False,
        seed=0,
        min_set_size=1,
        max_set_size=10,
        fi_curated_only=True,
        fi_min_score=None,
        allow_low_coverage=True,
    )
    defaults.update(overrides)
    return argparse.Namespace(**defaults)
