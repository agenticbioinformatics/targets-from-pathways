"""Smoke test for run_pipeline.py (the Stages 1-6 orchestrator).

Stage 1 needs real Open Targets / Reactome data, so this uses
``--from-manifest`` on a copy of the checked-in synthetic Stage 1 output
(``example_data/example_run/``) and checks that Stages 2-6 + the report
chain through and produce the expected artifacts. It is a "does it wire
up" check, not a correctness check on rankings.
"""

from __future__ import annotations

import shutil

import pandas as pd
import pytest

import run_pipeline

EXAMPLE_RUN = run_pipeline.PIPELINE_DIR.parent / "example_data" / "example_run"
STAGE1_FILES = [
    "manifest.json", "genes.parquet", "gene_sets.parquet", "interactions.parquet",
    "ot_disease_subset.parquet", "ot_target_subset.parquet",
]


@pytest.mark.skipif(
    not (EXAMPLE_RUN / "manifest.json").exists(),
    reason="example_data/example_run not populated",
)
def test_from_manifest_chains_stages_and_writes_report(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    for name in STAGE1_FILES:
        shutil.copy(EXAMPLE_RUN / name, run_dir / name)

    run_pipeline.main(["--from-manifest", str(run_dir / "manifest.json"), "--report"])

    # every artifact Stages 2-6 + the report are expected to drop
    expected = {
        "stage 2": ["disease_pathways.tsv"],
        "stage 3": ["graph_weight.npz", "graph_sign.npz", "graph_gene_index.json", "graph_metadata.json"],
        "stage 4": ["gene_weights.tsv"],  # + graph_* rewritten in place
        "stage 5": ["candidate_scores.tsv"],
        "stage 6": ["candidates_annotated.tsv"],
        "report": ["report.html"],
    }
    for stage, names in expected.items():
        for name in names:
            assert (run_dir / name).exists(), f"{stage}: {name} not produced"

    annotated = pd.read_csv(run_dir / "candidates_annotated.tsv", sep="\t")
    assert len(annotated) == 16  # 17 graph genes minus the target
    assert {"composite_score", "composite_breakdown", "composite_weights"} <= set(annotated.columns)
    assert annotated["composite_score"].between(0.0, 1.0).all()

    html = (run_dir / "report.html").read_text()
    assert "http://" not in html and "https://" not in html
    assert html.count("<details ") == len(annotated)


def test_missing_required_args_exits():
    with pytest.raises(SystemExit):
        run_pipeline.main(["--target", "PTGS2"])  # no --disease/--data-dir/--out-dir, no --from-manifest
