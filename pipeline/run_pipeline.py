"""Pipeline orchestrator — chains Stages 1-6 for one (target, disease).

Each stage is run as its own subprocess (the same CLI documented in the
root README), so a stage failure surfaces its own error and stops the run
with a non-zero exit rather than being swallowed.

Two entry modes:

- ``--target X --disease Y --data-dir DIR --out-dir DIR`` — full run,
  starting from Stage 1 (needs the pinned Open Targets / Reactome data in
  ``--data-dir``, or a network connection for Stage 1 to download it).
- ``--from-manifest path/to/manifest.json`` — skip Stage 1 and start from
  an existing Stage 1 output directory. This is the offline path: e.g.
  ``example_data/example_run/`` is a hand-built Stage 1 output
  (``example_data/build_gsea_example.py``) that needs no real data.

Stages 2-6 all write into the manifest's directory, so that directory is
the "run directory" the report (``stage7_report.py``) then reads.

Any command-line option this orchestrator does not recognise is forwarded
verbatim to ``stage1_ingest.py`` — that stage is the one with a rich
config surface. So a full run that needs, say, a lower coverage floor is
just::

    python3 pipeline/run_pipeline.py --target NOD2 --disease MONDO_0005265 \\
        --data-dir data/ --out-dir out/ --report --allow-low-coverage

(``--allow-low-coverage`` is not a ``run_pipeline`` flag; it passes
through.) Forwarding an unknown option together with ``--from-manifest``
is an error, since Stage 1 does not run in that mode.
"""

from __future__ import annotations

import argparse
import logging
import subprocess
import sys
from pathlib import Path

PIPELINE_DIR = Path(__file__).resolve().parent

logger = logging.getLogger("run_pipeline")


def _fail(msg: str) -> None:
    logger.error(msg)
    sys.exit(1)


def _run_stage(label: str, script: str, args: list, *, python_exe: str) -> None:
    cmd = [python_exe, str(PIPELINE_DIR / script), *(str(a) for a in args)]
    logger.info("Running %s …", label)
    # stdout/stderr are NOT captured — the stage streams its own progress and,
    # on failure, its own error (e.g. Stage 1's "Could not resolve --target
    # 'FOO'. Did you mean: ...?") straight to the console.
    result = subprocess.run(cmd)
    if result.returncode != 0:
        _fail(
            f"{label} [{script}] failed with exit code {result.returncode}. "
            f"See its output above for the cause (e.g. an unresolved --target/--disease "
            f"prints Stage 1's suggestions).\n  command: {' '.join(cmd)}"
        )


def run_pipeline(
    *,
    target: str | None = None,
    disease: str | None = None,
    data_dir: Path | None = None,
    out_dir: Path | None = None,
    from_manifest: Path | None = None,
    method: str = "topology,rwr",
    weights: str | None = None,
    stage1_args: list[str] | None = None,
    python_exe: str | None = None,
) -> Path:
    """Run Stages 1-6 (or 2-6 with ``from_manifest``). Returns the path to
    ``candidates_annotated.tsv`` in the run directory. ``stage1_args`` is a
    list of extra CLI tokens forwarded to ``stage1_ingest.py`` (e.g.
    ``["--allow-low-coverage"]``)."""
    python_exe = python_exe or sys.executable
    stage1_args = stage1_args or []

    if from_manifest is not None:
        if stage1_args:
            _fail(
                f"{stage1_args} look like stage1_ingest.py options, but --from-manifest "
                f"skips Stage 1 — drop them or run a full pipeline."
            )
        manifest = Path(from_manifest).resolve()
        if not manifest.exists():
            _fail(f"--from-manifest {manifest} does not exist.")
        run_dir = manifest.parent
        logger.info("Starting from %s (Stage 1 skipped).", manifest)
    else:
        if not (target and disease and data_dir and out_dir):
            _fail(
                "--target, --disease, --data-dir and --out-dir are all required "
                "unless --from-manifest is given."
            )
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        _run_stage(
            "Stage 1 (ingest)",
            "stage1_ingest.py",
            ["--target", target, "--disease", disease, "--data-dir", data_dir,
             "--out-dir", out_dir, *stage1_args],
            python_exe=python_exe,
        )
        manifest = out_dir / "manifest.json"
        run_dir = out_dir

    _run_stage("Stage 2 (disease-pathway discovery)", "stage2_gsea_discovery.py",
               ["--manifest", manifest], python_exe=python_exe)
    _run_stage("Stage 3 (graph construction)", "stage3_build_graph.py",
               ["--manifest", manifest], python_exe=python_exe)
    _run_stage("Stage 4 (genetic-evidence weighting)", "stage4_genetic_evidence_weights.py",
               ["--manifest", manifest], python_exe=python_exe)
    _run_stage("Stage 5 (candidate scoring)", "stage5_score_candidates.py",
               ["--graph-dir", run_dir, "--method", method], python_exe=python_exe)
    stage6_args = ["--manifest", manifest]
    if weights:
        stage6_args += ["--weights", weights]
    _run_stage("Stage 6 (contextual annotation)", "stage6_annotate_context.py",
               stage6_args, python_exe=python_exe)

    annotated = run_dir / "candidates_annotated.tsv"
    if not annotated.exists():
        _fail(f"Pipeline finished but {annotated} is missing.")
    logger.info("Pipeline complete: %s", annotated)
    return annotated


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="run_pipeline.py",
        description="Chain Stages 1-6 for one (target, disease), or Stages 2-6 from an "
        "existing Stage 1 manifest. Any option not listed below is forwarded to "
        "stage1_ingest.py (run `python3 pipeline/stage1_ingest.py --help` for its flags, "
        "e.g. --allow-low-coverage, --no-download, --seed, --min-set-size).",
    )
    p.add_argument("--target", help="Gene symbol or Ensembl ID (Stage 1).")
    p.add_argument("--disease", help="Disease ontology ID, e.g. EFO_0005755 (Stage 1).")
    p.add_argument("--data-dir", type=Path, help="Pinned source cache for Stage 1.")
    p.add_argument("--out-dir", type=Path, help="Run directory Stage 1 writes into.")
    p.add_argument(
        "--from-manifest",
        type=Path,
        help="Skip Stage 1; start from this Stage 1 manifest.json (its directory is the run directory).",
    )
    p.add_argument("--method", default="topology,rwr", help="Stage 5 --method (default: topology,rwr).")
    p.add_argument("--weights", default=None, help="Stage 6 --weights (default: Stage 6's own default).")
    p.add_argument(
        "--report",
        action="store_true",
        help="After the pipeline, also write report.html in the run directory (via stage7_report.py).",
    )
    return p


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    args, stage1_extra = build_arg_parser().parse_known_args(argv)
    annotated = run_pipeline(
        target=args.target,
        disease=args.disease,
        data_dir=args.data_dir,
        out_dir=args.out_dir,
        from_manifest=args.from_manifest,
        method=args.method,
        weights=args.weights,
        stage1_args=stage1_extra,
    )
    if args.report:
        import stage7_report

        report = stage7_report.build_report(annotated.parent, annotated.parent / "report.html")
        logger.info("Report: %s", report)


if __name__ == "__main__":
    main()
