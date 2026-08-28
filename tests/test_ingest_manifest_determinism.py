"""Manifest determinism (README.md's Stage 1 testing requirement):

two runs over identical inputs produce byte-identical manifests apart from
created_at/run_id, and every artifact sha256 recorded in the manifest
matches the file on disk.

run_id is itself a deterministic hash of (target, disease, seed,
cli_parameters) — see stage1_ingest.py's _run_id — and cli_parameters includes
out_dir, so two runs into *different* out-dirs legitimately get different
run_ids even though target/disease/seed are identical. Using two distinct
tmp out-dirs here (the realistic way to run this pipeline twice) is exactly
why run_id has to be excluded alongside created_at, not just timestamp
non-determinism.
"""

import hashlib
import json
from pathlib import Path

from conftest import make_args
from stage1_ingest import run_ingest


def _sha256_file(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _normalize(manifest: dict, out_dir) -> dict:
    manifest = dict(manifest)
    manifest.pop("created_at", None)
    manifest.pop("run_id", None)
    manifest["cli_parameters"] = dict(manifest["cli_parameters"])
    manifest["cli_parameters"]["out_dir"] = "NORMALIZED_OUT_DIR"
    manifest["output_artifacts"] = [
        {**a, "path": a["path"].replace(str(out_dir), "NORMALIZED_OUT_DIR")} for a in manifest["output_artifacts"]
    ]
    return manifest


def test_two_runs_over_identical_inputs_produce_byte_identical_manifests(tmp_path, patched_acquisition):
    data_dir = tmp_path / "data"
    out_dir_a = tmp_path / "out_a"
    out_dir_b = tmp_path / "out_b"

    run_ingest(make_args(out_dir=out_dir_a, data_dir=data_dir))
    run_ingest(make_args(out_dir=out_dir_b, data_dir=data_dir))

    manifest_a = json.loads((out_dir_a / "manifest.json").read_text())
    manifest_b = json.loads((out_dir_b / "manifest.json").read_text())

    assert _normalize(manifest_a, out_dir_a) == _normalize(manifest_b, out_dir_b)
    # Both real, non-empty timestamps/ids, just not required to match.
    assert manifest_a["created_at"] and manifest_b["created_at"]
    assert manifest_a["run_id"] and manifest_b["run_id"]


def test_every_manifest_artifact_sha256_matches_the_file_on_disk(tmp_path, patched_acquisition):
    out_dir = tmp_path / "out"
    run_ingest(make_args(out_dir=out_dir, data_dir=tmp_path / "data"))
    manifest = json.loads((out_dir / "manifest.json").read_text())

    assert manifest["output_artifacts"], "expected at least one recorded output artifact"
    for artifact in manifest["output_artifacts"]:
        path = Path(artifact["path"])
        assert path.exists(), f"artifact recorded in manifest is missing on disk: {path}"
        assert _sha256_file(path) == artifact["sha256"]

    for source in manifest["sources"]:
        for f in source["files"]:
            path = Path(f["path"])
            assert path.exists(), f"source file recorded in manifest is missing on disk: {path}"
            assert _sha256_file(path) == f["sha256"]
            assert path.stat().st_size == f["bytes"]
