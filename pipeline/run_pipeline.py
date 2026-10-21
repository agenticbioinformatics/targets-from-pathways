"""Run Stages 1-6 for one disease, or Stages 2-6 from an existing Stage 1 manifest (--from-manifest).

Each stage runs as its own subprocess with its normal CLI, so a failing stage prints its own error and
stops the run. All stages write into the manifest's directory, which --report then turns into report.html.
Options this script does not know are passed on to stage1_ingest.py (e.g. --allow-low-coverage).
"""

import argparse
import logging
import subprocess
import sys
from pathlib import Path

from stage0_schemas import fail
from stage7_report import build_report

PIPELINE_DIR = Path(__file__).resolve().parent

logger = logging.getLogger("run_pipeline")


def run_stage(label, script, args, python_exe):
    command = [python_exe, str(PIPELINE_DIR / script)]
    for arg in args:
        command.append(str(arg))
    logger.info("Running %s ...", label)
    # Output is not captured, so the stage's own progress and error messages reach the console.
    result = subprocess.run(command)
    if result.returncode != 0:
        fail(f"{label} [{script}] failed with exit code {result.returncode}. See its output above for the cause.\n"
             f"  command: {' '.join(command)}")


def run_pipeline(disease=None, data_dir=None, out_dir=None, from_manifest=None,
                 method="bfwalk", stage1_args=None, python_exe=None, node_weights=None, edge_weights=None):
    """Return the path of candidates_annotated.tsv. stage1_args are extra options for stage1_ingest.py;
    node_weights / edge_weights are Stage 4's optional user weight files."""
    python_exe = python_exe or sys.executable
    stage1_args = stage1_args or []

    if from_manifest is not None:
        if stage1_args:
            fail(f"{stage1_args} look like stage1_ingest.py options, but --from-manifest skips Stage 1.")
        manifest = Path(from_manifest).resolve()
        if not manifest.exists():
            fail(f"--from-manifest {manifest} does not exist.")
        logger.info("Starting from %s (Stage 1 skipped).", manifest)
    else:
        if not (disease and data_dir and out_dir):
            fail("--disease, --data-dir and --out-dir are all required unless --from-manifest is given.")
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        run_stage("Stage 1 (ingest)", "stage1_ingest.py",
                  ["--disease", disease, "--data-dir", data_dir, "--out-dir", out_dir] + stage1_args,
                  python_exe)
        manifest = Path(out_dir) / "manifest.json"
    run_dir = manifest.parent

    run_stage("Stage 2 (disease-pathway discovery)", "stage2_gsea_discovery.py", ["--manifest", manifest], python_exe)
    run_stage("Stage 3 (graph construction)", "stage3_build_graph.py", ["--manifest", manifest], python_exe)
    stage4_args = ["--manifest", manifest]
    if node_weights is not None:
        stage4_args += ["--node-weights", Path(node_weights).resolve()]
    if edge_weights is not None:
        stage4_args += ["--edge-weights", Path(edge_weights).resolve()]
    run_stage("Stage 4 (optional user weights)", "stage4_genetic_evidence_weights.py", stage4_args, python_exe)
    run_stage("Stage 5 (candidate scoring)", "stage5_score_candidates.py",
              ["--graph-dir", run_dir, "--method", method], python_exe)
    run_stage("Stage 6 (contextual annotation)", "stage6_annotate_context.py", ["--manifest", manifest], python_exe)

    annotated = run_dir / "candidates_annotated.tsv"
    if not annotated.exists():
        fail(f"Pipeline finished but {annotated} is missing.")
    logger.info("Pipeline complete: %s", annotated)
    return annotated


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    parser = argparse.ArgumentParser(
        prog="run_pipeline.py",
        description="Chain Stages 1-6 for one disease, or Stages 2-6 from an existing Stage 1 manifest. "
        "Any option not listed below is forwarded to stage1_ingest.py (see `python3 pipeline/stage1_ingest.py "
        "--help`, e.g. --allow-low-coverage, --no-download, --seed, --min-set-size).",
    )
    parser.add_argument("--disease", help="Disease ontology ID, e.g. MONDO_0005265 (Stage 1).")
    parser.add_argument("--data-dir", type=Path, help="Source data cache for Stage 1.")
    parser.add_argument("--out-dir", type=Path, help="Run directory Stage 1 writes into.")
    parser.add_argument("--from-manifest", type=Path,
                        help="Skip Stage 1 and start from this manifest.json (its directory is the run directory).")
    parser.add_argument("--method", default="bfwalk",
                        help="Stage 5 --method: bfwalk and/or rwr, comma-separated (default: bfwalk).")
    parser.add_argument("--node-weights", type=Path,
                        help="Stage 4 --node-weights: optional TSV of gene_id, weight to weight the graph with.")
    parser.add_argument("--edge-weights", type=Path,
                        help="Stage 4 --edge-weights: optional TSV of gene_a, gene_b, weight to weight the graph with.")
    parser.add_argument("--report", action="store_true", help="Also write report.html in the run directory (Stage 7).")
    args, stage1_args = parser.parse_known_args(argv)

    annotated = run_pipeline(disease=args.disease, data_dir=args.data_dir, out_dir=args.out_dir,
                             from_manifest=args.from_manifest, method=args.method,
                             stage1_args=stage1_args, node_weights=args.node_weights,
                             edge_weights=args.edge_weights)
    if args.report:
        build_report(annotated.parent, annotated.parent / "report.html")


if __name__ == "__main__":
    main()
