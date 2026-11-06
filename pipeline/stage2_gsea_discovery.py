"""Stage 2 - disease-pathway discovery with GSEA.

Ranks the disease's genes by genetic evidence and runs preranked GSEA against the Reactome gene sets listed in
Stage 1's manifest, with the DisPathCatalogue v2 method (gsea_core.py): every gene of the gene sets is in the
ranking (genes without evidence score 0), enrichment scores on the mean-centred ranking, p-values from a gene-set
permutation null per set size, Benjamini-Hochberg FDR over all tested sets. A pathway is disease-relevant when
it is enriched at the top of the ranking (es > 0) with FDR below --fdr-threshold. The validation catalogue
(validation/build_graph_catalogue.py) selects its pathways the same way, from the same code.

Writes next to the manifest:
- disease_pathways.tsv: the disease-relevant pathways; Stage 3 builds the graph from them.
- gsea_results.tsv: every tested pathway with es, nes, pval, fdr, n_perm and size, both ES signs.
- disease_seed_genes.tsv: the genes with genetic evidence > 0 (minus any benchmark holdout); Stage 5 propagates
  from them.
- benchmark_holdout_genes.tsv: the holdout (empty without --benchmark-holdout-file), so Stage 4 can ignore its
  evidence.
"""

import argparse
import logging
from pathlib import Path

import pandas as pd

import gsea_core
from stage0_schemas import (
    DiseasePathwaysSchema,
    GeneSetsSchema,
    GenesSchema,
    HoldoutGenesSchema,
    OTAssociationsSchema,
    SeedGenesSchema,
    datatype_harmonic_sum,
    fail,
    find_artifact,
    load_manifest,
)
from stage1_ingest import build_symbol_index, resolve_gene

logger = logging.getLogger("gsea_discovery")

NO_GSEA_CONSEQUENCE = ("disease_pathways.tsv is written empty, so Stage 3 has no pathway to build the graph "
                       "from and stops.")
GSEA_RESULT_COLUMNS = ["set_id", "gene_set", "source_db", "es", "nes", "pval", "fdr", "n_perm", "geneset_size"]


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="stage2_gsea_discovery.py",
        description="Stage 2 - disease-pathway discovery. GSEA (DisPathCatalogue v2 method: permutation p-values, "
        "BH-FDR) of the disease's genetic-association ranking against Stage 1's Reactome gene sets. Reads "
        "gene_sets.parquet + ot_disease_subset.parquet via --manifest; writes disease_pathways.tsv, "
        "gsea_results.tsv, disease_seed_genes.tsv and benchmark_holdout_genes.tsv next to it.",
    )
    parser.add_argument("--manifest", required=True, type=Path, help="Path to Stage 1's manifest.json.")
    parser.add_argument("--fdr-threshold", type=float, default=0.05,
                        help="BH-FDR cutoff for a pathway enriched at the top of the ranking (es > 0). "
                        "Default: 0.05, as in the validation catalogue.")
    parser.add_argument("--min-genes", type=int, default=gsea_core.MIN_SIZE,
                        help=f"Run GSEA only if at least this many genes have genetic evidence and are in a gene "
                        f"set. Default: {gsea_core.MIN_SIZE}.")
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
    """Each gene's genetic_association score: its datasource scores aggregated as Open Targets aggregates a
    datatype (datatype_harmonic_sum). Genetic evidence only: Open Targets' overall score includes pathway
    evidence, and testing that against pathways would be circular. Only genes scoring > 0 are returned."""
    genetic = ot_disease_subset[ot_disease_subset["datatype_id"] == "genetic_association"]
    signature = genetic.groupby("gene_id")["score"].agg(datatype_harmonic_sum).reset_index()
    signature = signature[signature["score"] > 0]
    if signature.empty:
        fail("ot_disease_subset.parquet has no genetic_association score > 0 - cannot build a GSEA ranking signature.")
    return signature.rename(columns={"gene_id": "gene"})[["gene", "score"]].reset_index(drop=True)


def gene_sets_as_library(gene_sets):
    """{set_id: sorted gene IDs}, in set_id order."""
    return {set_id: sorted(genes) for set_id, genes in gene_sets.groupby("set_id")["gene_id"]}


def no_tested_pathways():
    return pd.DataFrame({column: pd.Series(dtype="float64") for column in GSEA_RESULT_COLUMNS}).astype(
        {"set_id": "object", "gene_set": "object", "source_db": "object", "n_perm": "int64", "geneset_size": "int64"})


def run_gsea(signature, gene_sets, seed, min_genes=gsea_core.MIN_SIZE):
    """Every tested pathway's statistics (unfiltered, both ES signs). Empty, with a warning saying why, when
    there is nothing to test."""
    library = gene_sets_as_library(gene_sets)
    background = set(gene_sets["gene_id"])
    n_in_background = int(signature["gene"].isin(background).sum())
    if n_in_background < min_genes:
        logger.warning("Only %d gene(s) with genetic evidence are in a gene set (fewer than %d), so GSEA does not "
                       "run. %s", n_in_background, min_genes, NO_GSEA_CONSEQUENCE)
        return no_tested_pathways()
    scores = dict(zip(signature["gene"], signature["score"]))
    result = gsea_core.run_gsea(scores, library, seed=seed, order_seed=seed)
    if result.empty:
        logger.warning("No gene set has %d to %d genes, so GSEA tested nothing. %s", gsea_core.MIN_SIZE,
                       gsea_core.MAX_SIZE, NO_GSEA_CONSEQUENCE)
        return no_tested_pathways()
    set_names = gene_sets.drop_duplicates("set_id").set_index("set_id")[["set_name", "source_db"]]
    result = result.join(set_names, on="set_id").rename(columns={"set_name": "gene_set"})
    return result[GSEA_RESULT_COLUMNS]


def filter_results(result, fdr_threshold):
    """The disease-relevant pathways: enriched at the top of the ranking (es > 0) with FDR below the cutoff.
    A negative ES means the set's genes sit among the weakest-evidence genes, not that it is depleted."""
    return gsea_core.significant_sets(result, fdr_threshold)


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
    result = run_gsea(signature, gene_sets, seed=manifest.seed, min_genes=args.min_genes)
    significant = filter_results(result, args.fdr_threshold)
    logger.info("GSEA: %d genes with genetic evidence (%d in a gene set), %d/%d gene sets tested, %d with es > 0 "
                "and FDR < %g.", len(signature), int(signature["gene"].isin(set(gene_sets["gene_id"])).sum()),
                len(result), gene_sets["set_id"].nunique(), len(significant), args.fdr_threshold)

    out_dir = args.manifest.parent
    result.sort_values(["fdr", "pval", "set_id"], na_position="last").to_csv(
        out_dir / "gsea_results.tsv", sep="\t", index=False, na_rep="NA")
    output = significant[list(DiseasePathwaysSchema.columns)].copy()
    output["pval"] = output["pval"].map(round_to_4_significant_digits)
    output["fdr"] = output["fdr"].map(round_to_4_significant_digits)
    output = DiseasePathwaysSchema.validate(output)
    output_path = out_dir / "disease_pathways.tsv"
    output.to_csv(output_path, sep="\t", index=False)

    seeds = SeedGenesSchema.validate(signature.rename(columns={"gene": "gene_id"}))
    seeds = seeds.sort_values(["score", "gene_id"], ascending=[False, True])
    seeds_path = out_dir / "disease_seed_genes.tsv"
    seeds.to_csv(seeds_path, sep="\t", index=False)
    # Always written, so a later run without a holdout replaces an earlier run's file.
    holdout = HoldoutGenesSchema.validate(pd.DataFrame({"gene_id": pd.Series(sorted(holdout_gene_ids), dtype="object")}))
    holdout.to_csv(out_dir / "benchmark_holdout_genes.tsv", sep="\t", index=False)
    logger.info("Stage 2 complete: %s, %s (%d seed genes, %d held out)", output_path, seeds_path, len(seeds),
                len(holdout))
    return output_path


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    run(parse_args(argv))


if __name__ == "__main__":
    main()
