"""Smoke tests for run_pipeline.py (the Stages 1-6 orchestrator).

Stage 1 needs real Open Targets / Reactome data, so the chaining test uses the session's
``synthetic_full_run`` (conftest.py: ``run_pipeline.py --from-manifest --report`` on the synthetic Stage 1
output) and checks that Stages 2-6 + the report produce the expected artifacts. It is a "does it wire
up" check, not a correctness check on rankings.
"""

from __future__ import annotations

import subprocess
import sys

import pandas as pd
import pytest

import run_pipeline


def test_from_manifest_chains_stages_and_writes_report(synthetic_full_run):
    run_dir = synthetic_full_run
    # every artifact Stages 2-6 + the report are expected to drop
    expected = {
        "stage 2": ["disease_pathways.tsv", "disease_seed_genes.tsv", "benchmark_holdout_genes.tsv"],
        "stage 3": ["graph_weight.npz", "graph_sign.npz", "graph_gene_index.json", "graph_metadata.json"],
        "stage 4": ["genetic_evidence_scores.tsv"],  # + graph files rewritten in place
        "stage 5": ["candidate_scores.tsv"],
        "stage 6": ["candidates_annotated.tsv"],
        "report": ["report.html"],
    }
    for stage, names in expected.items():
        for name in names:
            assert (run_dir / name).exists(), f"{stage}: {name} not produced"

    annotated = pd.read_csv(run_dir / "candidates_annotated.tsv", sep="\t")
    assert len(annotated) == 3  # 19 genes on the 3 GSEA pathways minus 16 seed genes among them
    # default Stage 5 method, and Stage 6 keeps its ranking - there is no composite score
    assert "bfwalk_score" in annotated.columns and annotated["bfwalk_score"].is_monotonic_decreasing
    assert not any(column.startswith("composite") for column in annotated.columns)

    html = (run_dir / "report.html").read_text()
    assert "http://" not in html and "https://" not in html
    assert html.count('class="trace-entry"') == len(annotated)


def test_missing_required_args_exits():
    with pytest.raises(SystemExit):
        run_pipeline.main(["--disease", "MONDO_0005265"])  # no --data-dir/--out-dir, no --from-manifest


def test_stage_failure_names_the_stage_and_command(tmp_path, caplog):
    """A stage exiting non-zero -> a helpful ERROR line (which stage, its
    command, "see output above") and exit code 1, not a bare traceback."""
    (tmp_path / "manifest.json").write_text('{"not": "a valid manifest"}')

    with caplog.at_level("ERROR"):
        with pytest.raises(SystemExit) as exc:
            run_pipeline.main(["--from-manifest", str(tmp_path / "manifest.json")])

    assert exc.value.code == 1
    msg = "\n".join(r.message for r in caplog.records)
    assert "Stage 2" in msg
    assert "stage2_gsea_discovery.py" in msg
    assert "exit code 1" in msg
    assert "command:" in msg


def test_unknown_flags_forward_to_stage1(monkeypatch, tmp_path):
    """An option run_pipeline doesn't define is passed through to
    stage1_ingest.py (so `--allow-low-coverage` etc. need no run_pipeline
    flag). Later stages are not run here — run_stage is stubbed."""
    calls = []
    monkeypatch.setattr(run_pipeline, "run_stage",
                        lambda label, script, args, python_exe: calls.append((script, [str(a) for a in args])))

    with pytest.raises(SystemExit):  # stages are stubbed, so the final artifact check fails — fine
        run_pipeline.main([
            "--disease", "MONDO_0005265",
            "--data-dir", str(tmp_path / "d"), "--out-dir", str(tmp_path / "o"),
            "--allow-low-coverage", "--min-set-size", "5",
        ])

    stage1_script, stage1_args = calls[0]
    assert stage1_script == "stage1_ingest.py"
    assert "--allow-low-coverage" in stage1_args
    assert stage1_args[stage1_args.index("--min-set-size") + 1] == "5"
    # run_pipeline's own flags are NOT forwarded
    assert "--report" not in stage1_args and "--method" not in stage1_args


def test_stage1_flags_with_from_manifest_are_rejected(tmp_path, caplog):
    (tmp_path / "manifest.json").write_text("{}")
    with caplog.at_level("ERROR"):
        with pytest.raises(SystemExit):
            run_pipeline.main(["--from-manifest", str(tmp_path / "manifest.json"), "--allow-low-coverage"])
    assert "allow-low-coverage" in "\n".join(r.message for r in caplog.records)


def test_help_runs_the_installation_check(capsys):
    with pytest.raises(SystemExit) as exc:
        run_pipeline.main(["--help"])
    assert exc.value.code == 0
    assert "all stages' packages are" in capsys.readouterr().out


def test_missing_package_is_reported_before_help(monkeypatch, caplog):
    real_import = run_pipeline.importlib.import_module

    def import_module(name):
        if name == "stage2_gsea_discovery":
            raise ModuleNotFoundError("No module named 'blitzgsea'")
        return real_import(name)

    monkeypatch.setattr(run_pipeline.importlib, "import_module", import_module)
    with caplog.at_level("ERROR"):
        with pytest.raises(SystemExit) as exc:
            run_pipeline.main(["--help"])
    assert exc.value.code == 1
    message = "\n".join(r.message for r in caplog.records)
    assert "stage2_gsea_discovery.py: No module named 'blitzgsea'" in message
    assert "pip install -r requirements.txt" in message


def test_help_works_as_a_subprocess():
    # The real command a user types to check the installation.
    result = subprocess.run([sys.executable, str(run_pipeline.PIPELINE_DIR / "run_pipeline.py"), "--help"],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "--from-manifest" in result.stdout
