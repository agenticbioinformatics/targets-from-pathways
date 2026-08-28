"""Stage 2 — Disease-pathway discovery via GSEA.

Consumes only a Stage 1 manifest (``--manifest``): resolves the
``gene_sets.parquet`` and ``ot_disease_subset.parquet`` artifacts it points
to, and the resolved target gene, from the manifest. Never opens a raw
source file — everything it needs was already normalized and canonicalized
by ``stage1_ingest.py``.

Two design choices are load-bearing for the pipeline's orthogonality
discipline (README.md, "Correcting Stage 2 for annotation bias" and Iter-2
of the review changelog) and are not incidental implementation detail:

1. The GSEA ranking signature is built from ``ot_disease_subset``'s
   ``genetic_association`` datatype rows only — never Open Targets'
   aggregated/overall association score, which is itself partly
   pathway-derived (it folds in the ``affected_pathway`` datatype). Ranking
   disease genes by the aggregated score and then testing them against
   Reactome pathways would validate pathway evidence against itself; Stage 4
   depends on the same non-pathway-datatype restriction later for the same
   reason.
2. Gene sets are collapsed on near-duplication (Jaccard > 0.7, keeping the
   smaller/more specific set) before testing, in addition to Stage 1's size
   cap. ``ReactomePathways.gmt`` encodes every level of Reactome's pathway
   hierarchy as a separate set, so a pathway and its near-identical parent
   both enter the library; testing both violates the independence
   Benjamini-Hochberg assumes and makes the correction anticonservative.

Neither choice resolves README.md's open decision #2 ("correcting Stage 2
for annotation bias") — BH-FDR over the collapsed sets remains option C, a
placeholder, and well-studied genes will still test significant in more
pathways than genuinely disease-specific ones. ``run_gsea`` returns every
tested pathway's stats before any significance filtering specifically so a
per-pathway specificity score (one of the other two options) can be added as
an extra column between ``run_gsea`` and ``filter_results`` without
reshaping either function.

The output artifact (``disease_pathways.tsv``) is a stable, documented
inter-stage contract, not an incidental side effect: it is validated against
``schemas.DiseasePathwaysSchema`` before being written, exactly like Stage
1's parquet artifacts, and Stage 3 (graph construction) and Stage 7
(benchmark validation) are expected to read it as that schema rather than
re-deriving its column layout from this file's source.

``--benchmark-holdout-file`` exists for Stage 7: a literature-curated
resistance/compensation pair (e.g. EGFR->MET) is only valid ground truth if
its genes were never part of Stage 2's own disease-associated seed set —
otherwise a benchmark gene's presence in a discovered pathway would be
circular (of course a seed gene's own pathway looks disease-relevant). The
holdout file lists genes one per line (``#`` comments and blank lines
skipped) and is resolved to Ensembl gene IDs via Stage 1's genes.parquet —
via the exact same ``resolve_target``/symbol-index logic ``stage1_ingest.py`` uses
for ``--target`` — never by string-matching against ``ot_disease_subset``,
which carries no symbol column to match against in the first place.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from collections import defaultdict
from pathlib import Path

import blitzgsea as blitz
import pandas as pd

from stage1_ingest import _build_symbol_index, resolve_target
from stage0_schemas import DiseasePathwaysSchema, GeneSetsSchema, GenesSchema, Manifest, OTAssociationsSchema

logger = logging.getLogger("gsea_discovery")

# README.md's Stage 2 open-decision note calls this threshold a placeholder;
# it is fixed here (not a CLI flag) because the task that scoped this script
# fixed it, not because it is settled.
JACCARD_COLLAPSE_THRESHOLD = 0.7

GENETIC_ASSOCIATION_DATATYPE = "genetic_association"

# DiseasePathwaysSchema is the single source of truth for both the output
# TSV's column set and their order.
OUTPUT_COLUMNS = list(DiseasePathwaysSchema.columns)


# ==========================================================================
# CLI
# ==========================================================================


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="stage2_gsea_discovery.py",
        description="Stage 2 — disease-pathway discovery. FDR-corrected GSEA of the disease's "
        "genetic-association genes against Reactome gene sets (size-capped, redundancy-collapsed). "
        "Reads gene_sets.parquet + ot_disease_subset.parquet via --manifest; writes "
        "disease_pathways.tsv next to it.",
    )
    p.add_argument("--manifest", required=True, type=Path, help="Path to Stage 1's manifest.json.")
    p.add_argument("--pval-threshold", type=float, default=0.05, help="Nominal p-value cutoff.")
    p.add_argument("--fdr-threshold", type=float, default=0.1, help="BH-FDR cutoff.")
    p.add_argument(
        "--benchmark-holdout-file",
        type=Path,
        default=None,
        help="Optional file of gene symbols/Ensembl IDs (one per line, '#' comments allowed) to "
        "resolve via Stage 1's genes.parquet and exclude from the disease-associated seed gene "
        "set before running GSEA — for Stage 7's benchmark validation, to keep known "
        "resistance-pair genes out of Stage 2's own evidence.",
    )
    return p


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def _fail(msg: str) -> None:
    logger.error(msg)
    sys.exit(1)


# ==========================================================================
# Manifest-driven artifact loading
# ==========================================================================


def load_manifest(manifest_path: Path) -> Manifest:
    return Manifest.model_validate(json.loads(manifest_path.read_text()))


def _resolve_artifact_path(manifest: Manifest, manifest_path: Path, filename: str) -> Path:
    """Find a Stage 1 output artifact by filename and resolve it to a real path.

    ``output_artifacts[i].path`` is whatever ``--out-dir`` Stage 1 was run
    with, which may be relative to a working directory that no longer
    matches this process's cwd — so a relative path is retried next to
    ``manifest.json`` itself before giving up.
    """
    for artifact in manifest.output_artifacts:
        if Path(artifact.path).name != filename:
            continue
        as_given = Path(artifact.path)
        if as_given.exists():
            return as_given
        alongside_manifest = manifest_path.parent / as_given.name
        if alongside_manifest.exists():
            return alongside_manifest
        _fail(
            f"Manifest artifact {filename!r} listed at {as_given} does not exist "
            f"(also checked {alongside_manifest})."
        )
    _fail(f"Manifest {manifest_path} lists no output artifact named {filename!r}.")


def load_gene_sets(manifest: Manifest, manifest_path: Path) -> pd.DataFrame:
    path = _resolve_artifact_path(manifest, manifest_path, "gene_sets.parquet")
    return GeneSetsSchema.validate(pd.read_parquet(path))


def load_disease_subset(manifest: Manifest, manifest_path: Path) -> pd.DataFrame:
    path = _resolve_artifact_path(manifest, manifest_path, "ot_disease_subset.parquet")
    return OTAssociationsSchema.validate(pd.read_parquet(path))


def load_genes(manifest: Manifest, manifest_path: Path) -> pd.DataFrame:
    path = _resolve_artifact_path(manifest, manifest_path, "genes.parquet")
    return GenesSchema.validate(pd.read_parquet(path))


# ==========================================================================
# Benchmark holdout (--benchmark-holdout-file)
# ==========================================================================


def resolve_benchmark_holdout(holdout_path: Path, genes_df: pd.DataFrame) -> set[str]:
    """Resolve a benchmark holdout file to a set of Ensembl gene IDs.

    Reuses stage1_ingest.py's ``resolve_target`` (symbol, synonym, or bare Ensembl
    ID, same as ``--target``) against Stage 1's genes.parquet — the ID
    authority — rather than matching raw holdout strings against anything in
    ot_disease_subset, which is already Ensembl-keyed and carries no symbol
    column to match against. An unresolvable entry fails the run loudly:
    a silently-dropped holdout gene would silently reintroduce the exact
    circularity this flag exists to prevent.
    """
    if not holdout_path.exists():
        _fail(f"--benchmark-holdout-file {holdout_path} does not exist.")

    index = _build_symbol_index(genes_df)
    holdout_ids: set[str] = set()
    unresolved: list[str] = []
    for line in holdout_path.read_text().splitlines():
        raw = line.split("#", 1)[0].strip()
        if not raw:
            continue
        gene_id, _symbol = resolve_target(raw, genes_df, index)
        if gene_id is None:
            unresolved.append(raw)
        else:
            holdout_ids.add(gene_id)
    if unresolved:
        _fail(
            f"--benchmark-holdout-file {holdout_path}: could not resolve "
            f"{len(unresolved)} entr{'y' if len(unresolved) == 1 else 'ies'} to a gene in "
            f"Open Targets: {unresolved}."
        )
    return holdout_ids


def apply_benchmark_holdout(ot_disease_subset: pd.DataFrame, holdout_gene_ids: set[str]) -> pd.DataFrame:
    excluded = set(ot_disease_subset["gene_id"]) & holdout_gene_ids
    if excluded:
        logger.info("Benchmark holdout: excluded %d gene(s) from the seed set.", len(excluded))
    return ot_disease_subset[~ot_disease_subset["gene_id"].isin(holdout_gene_ids)]


# ==========================================================================
# Signature (ranked disease-association scores)
# ==========================================================================


def build_signature(ot_disease_subset: pd.DataFrame) -> pd.DataFrame:
    """Build the preranked-GSEA signature from genetic_association evidence only.

    A gene can carry multiple genetic_association rows (one per datasource,
    e.g. eva, gene2phenotype, gwas_credible_sets) — these are collapsed to
    one ranking value per gene by taking the max, i.e. its single strongest
    piece of independent genetic evidence, not Open Targets' own
    cross-datasource aggregate.
    """
    genetic = ot_disease_subset[ot_disease_subset["datatype_id"] == GENETIC_ASSOCIATION_DATATYPE]
    if genetic.empty:
        _fail(
            "ot_disease_subset.parquet has no genetic_association rows — "
            "cannot build a GSEA ranking signature."
        )
    signature = genetic.groupby("gene_id", as_index=False)["score"].max()
    # Association scores are one-directional (0-1, "more evidence" only), unlike
    # a typical up/down expression signature — blitzgsea's ES is computed on
    # |value| after centering, so this still ranks genes correctly by evidence
    # strength; the sign of its ES/NES reflects clustering at the top vs.
    # bottom of that ranking, not up- vs down-regulation.
    return signature.rename(columns={"gene_id": "gene"})[["gene", "score"]]


# ==========================================================================
# Near-duplicate set collapsing
# ==========================================================================


def collapse_near_duplicate_sets(
    gene_sets: pd.DataFrame, jaccard_threshold: float = JACCARD_COLLAPSE_THRESHOLD
) -> pd.DataFrame:
    """Drop sets whose Jaccard similarity to an already-kept smaller set exceeds
    ``jaccard_threshold``, keeping the smaller/more specific set in each pair.

    Processes sets smallest-first so the first member of any near-duplicate
    cluster encountered is always the most specific one, then checks each
    later (larger-or-equal) candidate only against genes it shares with
    already-kept sets — an inverted-gene-index lookup, not an all-pairs scan,
    since Reactome's full hierarchy can carry thousands of sets.
    """
    sets_by_id: dict[str, frozenset[str]] = {
        set_id: frozenset(genes) for set_id, genes in gene_sets.groupby("set_id")["gene_id"]
    }
    order = sorted(sets_by_id, key=lambda set_id: (len(sets_by_id[set_id]), set_id))

    kept_ids: list[str] = []
    gene_to_kept_sets: dict[str, set[str]] = defaultdict(set)
    n_dropped = 0

    for set_id in order:
        genes = sets_by_id[set_id]
        candidate_ids: set[str] = set()
        for gene in genes:
            candidate_ids |= gene_to_kept_sets.get(gene, set())

        is_near_duplicate = False
        for candidate_id in candidate_ids:
            candidate_genes = sets_by_id[candidate_id]
            union_size = len(genes | candidate_genes)
            jaccard = len(genes & candidate_genes) / union_size if union_size else 0.0
            if jaccard > jaccard_threshold:
                is_near_duplicate = True
                break

        if is_near_duplicate:
            n_dropped += 1
            continue
        kept_ids.append(set_id)
        for gene in genes:
            gene_to_kept_sets[gene].add(set_id)

    if n_dropped:
        logger.info("Collapsed %d near-duplicate gene set(s) (Jaccard > %.2f).", n_dropped, jaccard_threshold)
    return gene_sets[gene_sets["set_id"].isin(kept_ids)]


# ==========================================================================
# GSEA
# ==========================================================================


def run_gsea(
    signature: pd.DataFrame, gene_sets: pd.DataFrame, seed: int
) -> tuple[pd.DataFrame, dict[str, frozenset[str]]]:
    """Run preranked GSEA (blitzgsea) of the signature against the collapsed
    gene sets and return every tested pathway's stats, unfiltered.

    Returns one row per pathway actually tested (i.e. after blitzgsea's own
    overlap-with-signature min/max size filtering), with columns ``set_id``,
    ``gene_set`` (set_name), ``source_db``, ``pval``, ``fdr``,
    ``contains_target`` — deliberately unfiltered by significance, so a
    later specificity-score column and its own threshold can be inserted
    here without reshaping ``filter_results``.
    """
    library = {set_id: list(genes) for set_id, genes in gene_sets.groupby("set_id")["gene_id"]}
    result = blitz.gsea(signature, library, seed=seed, verbose=False, progress=False)
    if result.empty:
        _fail("blitzgsea returned no tested pathways (no gene set had enough signature overlap).")

    result = result.reset_index().rename(columns={"Term": "set_id"})
    set_meta = gene_sets.drop_duplicates("set_id").set_index("set_id")[["set_name", "source_db"]]
    result = result.join(set_meta, on="set_id").rename(columns={"set_name": "gene_set"})

    sets_by_id = {set_id: frozenset(genes) for set_id, genes in gene_sets.groupby("set_id")["gene_id"]}
    return result, sets_by_id


def add_contains_target(result: pd.DataFrame, sets_by_id: dict[str, frozenset[str]], target_gene_id: str) -> pd.DataFrame:
    result = result.copy()
    result["contains_target"] = result["set_id"].map(lambda set_id: target_gene_id in sets_by_id[set_id])
    return result


def filter_results(result: pd.DataFrame, pval_threshold: float, fdr_threshold: float) -> pd.DataFrame:
    """Significance filter, applied last against named columns only — an
    extension point (see module docstring) for a future specificity-score
    threshold alongside pval_threshold/fdr_threshold."""
    mask = (result["pval"] <= pval_threshold) & (result["fdr"] <= fdr_threshold)
    return result.loc[mask].sort_values(["fdr", "pval"])


# ==========================================================================
# Orchestration
# ==========================================================================


def run(args: argparse.Namespace) -> Path:
    manifest = load_manifest(args.manifest)
    gene_sets = load_gene_sets(manifest, args.manifest)
    ot_disease_subset = load_disease_subset(manifest, args.manifest)
    target_gene_id = manifest.resolved_target.gene_id

    if args.benchmark_holdout_file is not None:
        genes_df = load_genes(manifest, args.manifest)
        holdout_gene_ids = resolve_benchmark_holdout(args.benchmark_holdout_file, genes_df)
        ot_disease_subset = apply_benchmark_holdout(ot_disease_subset, holdout_gene_ids)

    signature = build_signature(ot_disease_subset)
    n_sets_before = gene_sets["set_id"].nunique()
    collapsed_gene_sets = collapse_near_duplicate_sets(gene_sets)

    result, sets_by_id = run_gsea(signature, collapsed_gene_sets, seed=manifest.seed)
    result = add_contains_target(result, sets_by_id, target_gene_id)
    significant = filter_results(result, args.pval_threshold, args.fdr_threshold)
    logger.info(
        "GSEA: %d seed genes, %d/%d gene sets tested, %d significant.",
        len(signature), len(result), n_sets_before, len(significant),
    )

    out_path = args.manifest.parent / "disease_pathways.tsv"
    to_write = significant[OUTPUT_COLUMNS].copy()
    to_write[["pval", "fdr"]] = to_write[["pval", "fdr"]].apply(lambda s: s.map(lambda v: float(f"{v:.4g}")))
    to_write = DiseasePathwaysSchema.validate(to_write)
    to_write.to_csv(out_path, sep="\t", index=False)
    logger.info("Stage 2 complete: %s", out_path)
    return out_path


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    args = parse_args(argv)
    run(args)


if __name__ == "__main__":
    main()
