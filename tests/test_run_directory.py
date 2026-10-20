"""pytest tests for run-directory integrity across stages, on copies of the checked-in example run.

- ``find_artifact`` (stage0_schemas.py): a stage reads the Stage 1 files next to its --manifest, even when
  the manifest lists a path that exists elsewhere (a copied or moved run directory), and refuses a file
  whose sha256 no longer matches the manifest.
- Benchmark holdout (Stages 2 + 4): a held-out gene is not a seed, and its own Open Targets evidence does
  not feed Stage 4's scores or edge weights - otherwise it would be ranked by the very association the
  benchmark is meant to recover.
"""

from __future__ import annotations

import json
import shutil

import pandas as pd
import pytest
import scipy.sparse as sp

import stage2_gsea_discovery as gsea
import stage3_build_graph as bg
import stage4_genetic_evidence_weights as gew
from stage0_schemas import find_artifact, load_manifest, sha256_of_file

from test_run_pipeline import EXAMPLE_RUN, STAGE1_FILES

HOLDOUT_FILE = EXAMPLE_RUN.parent / "benchmark_holdout_example.txt"
HELD_OUT = "ENSG00000000017"  # OLDG17 in benchmark_holdout_example.txt

pytestmark = pytest.mark.skipif(
    not (EXAMPLE_RUN / "manifest.json").exists(),
    reason="example_data/example_run not populated",
)


def _copy_run(tmp_path, listed_dir=None):
    """Copy the example's Stage 1 files to tmp_path/run. listed_dir: rewrite the manifest's artifact paths
    to point there instead."""
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    for name in STAGE1_FILES:
        shutil.copy(EXAMPLE_RUN / name, run_dir / name)
    if listed_dir is not None:
        manifest = json.loads((run_dir / "manifest.json").read_text())
        for artifact in manifest["output_artifacts"]:
            artifact["path"] = str(listed_dir / artifact["path"].split("/")[-1])
        (run_dir / "manifest.json").write_text(json.dumps(manifest))
    return run_dir


def _replace_artifact(run_dir, name, df):
    """Rewrite a Stage 1 artifact and record its new sha256, as if Stage 1 had written it."""
    path = run_dir / name
    df.to_parquet(path, index=False)
    manifest = json.loads((run_dir / "manifest.json").read_text())
    for artifact in manifest["output_artifacts"]:
        if artifact["path"].split("/")[-1] == name:
            artifact["sha256"] = sha256_of_file(path)
    (run_dir / "manifest.json").write_text(json.dumps(manifest))


# --------------------------------------------------------------------------
# find_artifact
# --------------------------------------------------------------------------


def test_copied_run_reads_its_own_files_not_the_listed_ones(tmp_path):
    # The manifest lists the original run's files, which exist - the copy next to the manifest must win.
    run_dir = _copy_run(tmp_path, listed_dir=EXAMPLE_RUN)
    manifest_path = run_dir / "manifest.json"
    assert find_artifact(load_manifest(manifest_path), manifest_path, "genes.parquet") == run_dir / "genes.parquet"


def test_listed_path_is_the_fallback(tmp_path):
    run_dir = _copy_run(tmp_path, listed_dir=EXAMPLE_RUN)
    (run_dir / "genes.parquet").unlink()
    manifest_path = run_dir / "manifest.json"
    assert find_artifact(load_manifest(manifest_path), manifest_path, "genes.parquet") == EXAMPLE_RUN / "genes.parquet"


def test_artifact_changed_after_stage1_is_refused(tmp_path):
    run_dir = _copy_run(tmp_path)
    ot = pd.read_parquet(run_dir / "ot_disease_subset.parquet")
    ot.iloc[:-1].to_parquet(run_dir / "ot_disease_subset.parquet", index=False)
    manifest_path = run_dir / "manifest.json"
    with pytest.raises(SystemExit):
        find_artifact(load_manifest(manifest_path), manifest_path, "ot_disease_subset.parquet")


# --------------------------------------------------------------------------
# benchmark holdout
# --------------------------------------------------------------------------


def _run_stages_2_to_4(run_dir, holdout=True):
    manifest = str(run_dir / "manifest.json")
    stage2_argv = ["--manifest", manifest]
    if holdout:
        stage2_argv += ["--benchmark-holdout-file", str(HOLDOUT_FILE)]
    gsea.run(gsea.parse_args(stage2_argv))
    bg.run(bg.parse_args(["--manifest", manifest]))
    gew.run(gew.parse_args(["--manifest", manifest]))


def test_holdout_gene_evidence_does_not_feed_stage4(tmp_path):
    run_dir = _copy_run(tmp_path)
    # Give the held-out gene strong genetic evidence: a leak would now be large.
    ot = pd.read_parquet(run_dir / "ot_disease_subset.parquet")
    ot.loc[ot["gene_id"] == HELD_OUT, "score"] = 0.9
    _replace_artifact(run_dir, "ot_disease_subset.parquet", ot)

    _run_stages_2_to_4(run_dir)

    assert HELD_OUT not in set(pd.read_csv(run_dir / "disease_seed_genes.tsv", sep="\t")["gene_id"])
    assert list(pd.read_csv(run_dir / "benchmark_holdout_genes.tsv", sep="\t")["gene_id"]) == [HELD_OUT]

    metadata = json.loads((run_dir / "graph_metadata.json").read_text())
    assert metadata["genetic_evidence_score"][HELD_OUT] == 0.0
    assert metadata["stage4_genetic_evidence"]["holdout_genes"] == [HELD_OUT]

    # Every edge touching the held-out gene is weighted as for a gene without evidence.
    index = json.loads((run_dir / "graph_gene_index.json").read_text())
    scores = metadata["genetic_evidence_score"]
    weight = sp.load_npz(run_dir / "graph_weight.npz").tocoo()
    i = index.index(HELD_OUT)
    for row, col, value in zip(weight.row, weight.col, weight.data):
        if i in (row, col):
            other = index[col] if row == i else index[row]
            assert value == pytest.approx(scores[other] / 2.0)


def test_rerun_without_holdout_clears_the_holdout_file(tmp_path):
    run_dir = _copy_run(tmp_path)
    _run_stages_2_to_4(run_dir, holdout=True)
    _run_stages_2_to_4(run_dir, holdout=False)
    assert pd.read_csv(run_dir / "benchmark_holdout_genes.tsv", sep="\t").empty
    metadata = json.loads((run_dir / "graph_metadata.json").read_text())
    assert metadata["genetic_evidence_score"][HELD_OUT] > 0.0


def test_stage4_rerun_gives_the_same_graph(tmp_path):
    # Stage 4 overwrites graph_weight.npz in place, so running it twice must not compound.
    run_dir = _copy_run(tmp_path)
    _run_stages_2_to_4(run_dir)
    first = sp.load_npz(run_dir / "graph_weight.npz")
    gew.run(gew.parse_args(["--manifest", str(run_dir / "manifest.json")]))
    second = sp.load_npz(run_dir / "graph_weight.npz")
    assert first.nnz == second.nnz
    assert abs(first - second).max() == 0
