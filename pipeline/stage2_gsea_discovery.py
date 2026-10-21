"""Stage 2 - disease-pathway discovery with GSEA.

Runs preranked GSEA of the disease's seed genes - the genes with an Open Targets genetic_association score
> 0 for the disease - against the Reactome gene sets listed in Stage 1's manifest and writes the significant
pathways to disease_pathways.tsv next to the manifest. Stage 3 builds the graph from those pathways.
The same seed genes (minus any benchmark holdout) are written to disease_seed_genes.tsv; Stage 5 propagates
from them. The holdout itself is written to benchmark_holdout_genes.tsv (empty without
--benchmark-holdout-file), so Stage 4 can ignore its evidence.
"""

import argparse
import logging
from pathlib import Path

import blitzgsea
import pandas as pd

from stage0_schemas import (
    DiseasePathwaysSchema,
    GeneSetsSchema,
    GenesSchema,
    HoldoutGenesSchema,
    OTAssociationsSchema,
    SeedGenesSchema,
    fail,
    find_artifact,
    load_manifest,
)
from stage1_ingest import build_symbol_index, resolve_gene

logger = logging.getLogger("gsea_discovery")

# Placeholder threshold, not tuned.
JACCARD_COLLAPSE_THRESHOLD = 0.7
# blitzgsea tests only gene sets with at least this many signature genes (its min_size default).
BLITZGSEA_MIN_SIZE = 5
NO_GSEA_CONSEQUENCE = ("disease_pathways.tsv is written empty, so Stage 3 has no pathway to build the graph "
                       "from and stops.")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="stage2_gsea_discovery.py",
        description="Stage 2 - disease-pathway discovery. FDR-corrected GSEA of the disease's "
        "genetic-association genes against Reactome gene sets (size-capped, redundancy-collapsed). "
        "Reads gene_sets.parquet + ot_disease_subset.parquet via --manifest; writes "
        "disease_pathways.tsv, disease_seed_genes.tsv and benchmark_holdout_genes.tsv next to it.",
    )
    parser.add_argument("--manifest", required=True, type=Path, help="Path to Stage 1's manifest.json.")
    parser.add_argument("--pval-threshold", type=float, default=0.05, help="Nominal p-value cutoff.")
    parser.add_argument("--fdr-threshold", type=float, default=0.1, help="BH-FDR cutoff.")
    parser.add_argument(
        "--benchmark-holdout-file", type=Path, default=None,
        help="Optional file of gene symbols/Ensembl IDs (one per line, '#' comments allowed) to exclude "
        "from the disease seed genes, so benchmark genes are not part of Stage 2's own evidence.",
    )
    return parser.parse_args(argv)


def read_holdout_genes(holdout_path, genes_df):
    """Resolve the holdout file to Ensembl gene IDs. An unresolvable entry stops the run, because silently
    dropping it would let a benchmark gene leak into the seed set."""
    if not holdout_path.exists():
        fail(f"--benchmark-holdout-file {holdout_path} does not exist.")
    index = build_symbol_index(genes_df)
    holdout_gene_ids = set()
    unresolved = []
    for line in holdout_path.read_text().splitlines():
        entry = line.split("#")[0].strip()
        if not entry:
            continue
        gene_id, _ = resolve_gene(entry, genes_df, index)
        if gene_id is None:
            unresolved.append(entry)
        else:
            holdout_gene_ids.add(gene_id)
    if unresolved:
        fail(f"--benchmark-holdout-file {holdout_path}: could not resolve {len(unresolved)} entries: {unresolved}.")
    return holdout_gene_ids


def build_signature(ot_disease_subset):
    """Rank the seed genes by genetic_association evidence only. Open Targets' overall score includes pathway
    evidence, and testing that against pathways would be circular. A gene's score is its strongest datasource;
    only genes scoring > 0 are seeds."""
    genetic = ot_disease_subset[ot_disease_subset["datatype_id"] == "genetic_association"]
    signature = genetic.groupby("gene_id", as_index=False)["score"].max()
    signature = signature[signature["score"] > 0]
    if signature.empty:
        fail("ot_disease_subset.parquet has no genetic_association score > 0 - cannot build a GSEA ranking signature.")
    return signature.rename(columns={"gene_id": "gene"})[["gene", "score"]]


def gene_sets_as_dict(gene_sets):
    sets_by_id = {}
    for set_id, genes in gene_sets.groupby("set_id")["gene_id"]:
        sets_by_id[set_id] = set(genes)
    return sets_by_id


def collapse_near_duplicate_sets(gene_sets, jaccard_threshold=JACCARD_COLLAPSE_THRESHOLD):
    """Drop a gene set if it is nearly identical (Jaccard > threshold) to a smaller set already kept.
    Reactome lists every hierarchy level as its own set, and testing near-copies breaks the FDR correction."""
    sets_by_id = gene_sets_as_dict(gene_sets)
    sizes_and_ids = []
    for set_id, genes in sets_by_id.items():
        sizes_and_ids.append((len(genes), set_id))

    kept_ids = []
    # gene -> kept sets containing it, so each set is compared only with kept sets that share a gene.
    kept_sets_by_gene = {}
    for _, set_id in sorted(sizes_and_ids):
        genes = sets_by_id[set_id]
        overlapping_kept_sets = set()
        for gene in genes:
            overlapping_kept_sets |= kept_sets_by_gene.get(gene, set())

        is_near_duplicate = False
        for kept_id in overlapping_kept_sets:
            kept_genes = sets_by_id[kept_id]
            jaccard = len(genes & kept_genes) / len(genes | kept_genes)
            if jaccard > jaccard_threshold:
                is_near_duplicate = True
                break
        if is_near_duplicate:
            continue

        kept_ids.append(set_id)
        for gene in genes:
            if gene not in kept_sets_by_gene:
                kept_sets_by_gene[gene] = set()
            kept_sets_by_gene[gene].add(set_id)

    n_dropped = len(sets_by_id) - len(kept_ids)
    if n_dropped:
        logger.info("Collapsed %d near-duplicate gene set(s) (Jaccard > %.2f).", n_dropped, jaccard_threshold)
    return gene_sets[gene_sets["set_id"].isin(kept_ids)]


def no_tested_pathways():
    return pd.DataFrame({"set_id": pd.Series(dtype="object"), "gene_set": pd.Series(dtype="object"),
                         "source_db": pd.Series(dtype="object"), "pval": pd.Series(dtype="float64"),
                         "fdr": pd.Series(dtype="float64")})


def run_gsea(signature, gene_sets, seed):
    """Return every tested pathway's statistics (unfiltered) and the gene sets as {set_id: genes}. When no
    pathway can be tested the statistics table is empty and a warning says why."""
    sets_by_id = gene_sets_as_dict(gene_sets)
    if len(signature) < BLITZGSEA_MIN_SIZE:
        logger.warning("Only %d seed gene(s): GSEA tests a pathway only if it holds at least %d of them, so "
                       "nothing can be tested. %s", len(signature), BLITZGSEA_MIN_SIZE, NO_GSEA_CONSEQUENCE)
        return no_tested_pathways(), sets_by_id

    library = {}
    for set_id, genes in sets_by_id.items():
        library[set_id] = sorted(genes)
    # Scores are one-sided (0-1), so the sign of the enrichment score means top vs bottom of the ranking,
    # not up- vs down-regulation.
    try:
        result = blitzgsea.gsea(signature, library, seed=seed, verbose=False, progress=False)
    except (ZeroDivisionError, ValueError) as error:
        # blitzgsea's null-distribution calibration breaks on small signatures: below ~7 genes, and whenever
        # the signature length equals one of its calibration set sizes (up to ~230 genes).
        logger.warning("blitzgsea could not calibrate GSEA on the %d-gene signature (%s: %s). %s",
                       len(signature), type(error).__name__, error, NO_GSEA_CONSEQUENCE)
        return no_tested_pathways(), sets_by_id
    if result.empty:
        logger.warning("No gene set holds %d or more of the %d seed genes, so GSEA tested nothing. %s",
                       BLITZGSEA_MIN_SIZE, len(signature), NO_GSEA_CONSEQUENCE)
        return no_tested_pathways(), sets_by_id

    result = result.reset_index().rename(columns={"Term": "set_id"})
    set_names = gene_sets.drop_duplicates("set_id").set_index("set_id")[["set_name", "source_db"]]
    result = result.join(set_names, on="set_id").rename(columns={"set_name": "gene_set"})
    return result, sets_by_id


def filter_results(result, pval_threshold, fdr_threshold):
    significant = (result["pval"] <= pval_threshold) & (result["fdr"] <= fdr_threshold)
    return result.loc[significant].sort_values(["fdr", "pval"])


def round_to_4_significant_digits(value):
    return float(f"{value:.4g}")


def run(args):
    manifest = load_manifest(args.manifest)
    gene_sets = GeneSetsSchema.validate(pd.read_parquet(find_artifact(manifest, args.manifest, "gene_sets.parquet")))
    ot_disease_subset = OTAssociationsSchema.validate(
        pd.read_parquet(find_artifact(manifest, args.manifest, "ot_disease_subset.parquet"))
    )

    holdout_gene_ids = set()
    if args.benchmark_holdout_file is not None:
        genes_df = GenesSchema.validate(pd.read_parquet(find_artifact(manifest, args.manifest, "genes.parquet")))
        holdout_gene_ids = read_holdout_genes(args.benchmark_holdout_file, genes_df)
        is_holdout = ot_disease_subset["gene_id"].isin(holdout_gene_ids)
        logger.info("Benchmark holdout: excluded %d gene(s) from the seed set.", ot_disease_subset.loc[is_holdout, "gene_id"].nunique())
        ot_disease_subset = ot_disease_subset[~is_holdout]

    signature = build_signature(ot_disease_subset)
    collapsed_gene_sets = collapse_near_duplicate_sets(gene_sets)
    result, _ = run_gsea(signature, collapsed_gene_sets, seed=manifest.seed)

    significant = filter_results(result, args.pval_threshold, args.fdr_threshold)
    logger.info("GSEA: %d seed genes, %d/%d gene sets tested, %d significant.",
                len(signature), len(result), gene_sets["set_id"].nunique(), len(significant))

    output = significant[list(DiseasePathwaysSchema.columns)].copy()
    output["pval"] = output["pval"].map(round_to_4_significant_digits)
    output["fdr"] = output["fdr"].map(round_to_4_significant_digits)
    output = DiseasePathwaysSchema.validate(output)
    output_path = args.manifest.parent / "disease_pathways.tsv"
    output.to_csv(output_path, sep="\t", index=False)

    seeds = SeedGenesSchema.validate(signature.rename(columns={"gene": "gene_id"}))
    seeds = seeds.sort_values(["score", "gene_id"], ascending=[False, True])
    seeds_path = args.manifest.parent / "disease_seed_genes.tsv"
    seeds.to_csv(seeds_path, sep="\t", index=False)
    # Always written, so a later run without a holdout replaces an earlier run's file.
    holdout = HoldoutGenesSchema.validate(pd.DataFrame({"gene_id": pd.Series(sorted(holdout_gene_ids), dtype="object")}))
    holdout.to_csv(args.manifest.parent / "benchmark_holdout_genes.tsv", sep="\t", index=False)
    logger.info("Stage 2 complete: %s, %s (%d seed genes, %d held out)", output_path, seeds_path, len(seeds),
                len(holdout))
    return output_path


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    run(parse_args(argv))


if __name__ == "__main__":
    main()
