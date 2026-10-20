"""Stage 6 - contextual annotation.

Adds to each Stage 5 candidate its Stage 4 genetic_evidence_score, an Open Targets tractability bucket
(clinical / discovery / unknown) and a safety flag (has_liabilities / unknown - no annotation is not proof of
safety). The annotations are for reading only: candidates keep Stage 5's ranking (bfwalk_score, else
rwr_score). Writes candidates_annotated.tsv. No rows are filtered and no tissue-expression data is used.
"""

import argparse
import json
import logging
from pathlib import Path

import pandas as pd
import pyarrow.compute as pa_compute
import pyarrow.dataset as pa_dataset

from stage0_schemas import (
    ENSEMBL_GENE_ID_PATTERN,
    AnnotatedCandidatesSchema,
    CandidateScoresSchema,
    fail,
    load_manifest,
    ranking_column,
)

logger = logging.getLogger("annotate_context")

CLINICAL_TRACTABILITY_IDS = {
    "approved drug",
    "phase 1 clinical",
    "phase 2 clinical",
    "phase 3 clinical",
    "advanced clinical",
    "advanced clinical trials",
}
OUTPUT_COLUMNS = [
    "gene", "rwr_score", "bfwalk_score", "genetic_evidence_score",
    "tractability", "safety", "n_safety_liabilities",
]


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="stage6_annotate_context.py",
        description="Stage 6 - contextual annotation. Joins an Open Targets tractability bucket and a safety "
        "flag onto Stage 5's candidates, keeping Stage 5's ranking. Reads "
        "candidate_scores.tsv + graph_metadata.json + the OT target parquet (via --manifest); writes "
        "candidates_annotated.tsv.",
    )
    parser.add_argument("--manifest", required=True, type=Path, help="Path to Stage 1's manifest.json.")
    parser.add_argument("--candidates", type=Path, default=None, help="Stage 5's candidate_scores.tsv. Default: next to --manifest.")
    parser.add_argument("--graph-dir", type=Path, default=None,
                        help="Directory with Stage 4's graph_metadata.json. Default: --manifest's directory.")
    parser.add_argument("--target-parquet", type=Path, default=None,
                        help="Open Targets target parquet file or directory. Default: taken from manifest.sources.")
    parser.add_argument("--out", type=Path, default=None, help="Output TSV. Default: <candidates dir>/candidates_annotated.tsv.")
    return parser.parse_args(argv)



def find_next_to_manifest(path, manifest_path):
    # Paths in the manifest may be relative to another working directory.
    if path.exists():
        return path
    next_to_manifest = manifest_path.parent / path.name
    if next_to_manifest.exists():
        return next_to_manifest
    fail(f"{path.name} not found at {path} (also checked {next_to_manifest}).")


def find_target_parquet_files(manifest, manifest_path, override):
    if override is not None:
        if override.is_dir():
            parts = sorted(override.glob("*.parquet"))
            if not parts:
                fail(f"--target-parquet {override} is a directory with no .parquet files.")
            return parts
        if not override.exists():
            fail(f"--target-parquet {override} does not exist.")
        return [override]

    # Open Targets "target" is a directory of part files, so match "target" anywhere in the path.
    paths = []
    for source in manifest.sources:
        if source.db != "opentargets":
            continue
        for source_file in source.files:
            lower_path = source_file.path.lower()
            if lower_path.endswith(".parquet") and "target" in lower_path and "disease" not in lower_path \
                    and "association" not in lower_path:
                paths.append(find_next_to_manifest(Path(source_file.path), manifest_path))
    if not paths:
        fail("manifest.sources lists no Open Targets target parquet file - pass --target-parquet.")
    return paths


def load_genetic_evidence(graph_dir):
    metadata_path = graph_dir / "graph_metadata.json"
    if not metadata_path.exists():
        logger.warning("No graph_metadata.json in %s - 'genetic_evidence' component unavailable.", graph_dir)
        return None
    scores = json.loads(metadata_path.read_text()).get("genetic_evidence_score")
    if not scores:
        logger.warning("graph_metadata.json has no genetic_evidence_score (Stage 4 not run?).")
        return None
    return scores


def list_length(value):
    # Parquet lists arrive as numpy arrays; a missing value may be None or a NaN float.
    try:
        return len(value)
    except TypeError:
        return 0


def tractability_bucket(entries):
    if list_length(entries) == 0:
        return "unknown"
    bucket = "unknown"
    for entry in entries:
        if bool(entry["value"]):
            if str(entry["id"]).lower() in CLINICAL_TRACTABILITY_IDS:
                return "clinical"
            bucket = "discovery"
    return bucket


def safety_flag(liabilities):
    n_liabilities = list_length(liabilities)
    if n_liabilities > 0:
        return "has_liabilities", n_liabilities
    return "unknown", 0


def build_annotations(target_parquet_paths, genes):
    """One row per gene; genes missing from the target table get unknown / unknown / 0."""
    paths = []
    for path in target_parquet_paths:
        paths.append(str(path))
    target_table = pa_dataset.dataset(paths, format="parquet").to_table(
        columns=["id", "tractability", "safetyLiabilities"],
        filter=pa_compute.field("id").isin(list(genes)),
    ).to_pandas()

    rows = []
    annotated_genes = set()
    for row in target_table.itertuples(index=False):
        safety, n_liabilities = safety_flag(row.safetyLiabilities)
        rows.append([row.id, tractability_bucket(row.tractability), safety, n_liabilities])
        annotated_genes.add(row.id)
    for gene in sorted(set(genes) - annotated_genes):
        rows.append([gene, "unknown", "unknown", 0])
    return pd.DataFrame(rows, columns=["gene", "tractability", "safety", "n_safety_liabilities"])


def run(args):
    manifest = load_manifest(args.manifest)
    candidates_path = find_next_to_manifest(args.candidates or args.manifest.parent / "candidate_scores.tsv", args.manifest)
    graph_dir = args.graph_dir or args.manifest.parent
    target_parquet_files = find_target_parquet_files(manifest, args.manifest, args.target_parquet)

    candidates = CandidateScoresSchema.validate(pd.read_csv(candidates_path, sep="\t"))
    for gene in candidates["gene"]:
        if not ENSEMBL_GENE_ID_PATTERN.match(gene):
            fail(f"candidate_scores.tsv has a non-Ensembl gene id: {gene}")
    logger.info("Annotating %d candidates from %s.", len(candidates), candidates_path)

    evidence_scores = load_genetic_evidence(graph_dir)
    if evidence_scores is not None:
        candidates["genetic_evidence_score"] = candidates["gene"].map(evidence_scores).astype(float)
        n_missing = int(candidates["genetic_evidence_score"].isna().sum())
        if n_missing:
            logger.warning("%d candidate(s) missing from genetic_evidence_score - filled with 0.0.", n_missing)
            candidates["genetic_evidence_score"] = candidates["genetic_evidence_score"].fillna(0.0)

    annotations = build_annotations(target_parquet_files, candidates["gene"].tolist())
    candidates = candidates.merge(annotations, on="gene", how="left")
    logger.info("OT annotation: tractability %d/%d, safety liabilities %d/%d.",
                (candidates["tractability"] != "unknown").sum(), len(candidates),
                (candidates["safety"] != "unknown").sum(), len(candidates))

    # Stage 5's ranking, restored in case --candidates points at a re-sorted file.
    rank_column = ranking_column(candidates.columns)
    candidates = candidates.sort_values(rank_column, ascending=False, kind="stable").reset_index(drop=True)

    columns = []
    for column in OUTPUT_COLUMNS:
        if column in candidates.columns:
            columns.append(column)
    candidates = AnnotatedCandidatesSchema.validate(candidates[columns])

    out_path = args.out or candidates_path.parent / "candidates_annotated.tsv"
    candidates.to_csv(out_path, sep="\t", index=False)
    logger.info("Stage 6 complete (ranked by %s): %s", rank_column, out_path)
    return out_path


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    run(parse_args(argv))


if __name__ == "__main__":
    main()
