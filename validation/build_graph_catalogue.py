"""Build a catalogue of disease graphs from the Open Targets disease-pathway GSEA catalogue (Reactome only).

The GSEA catalogue contains for >3000 diseases:
1. Pathways: the Reactome rows of genetic_disease_pathway_catalogue.parquet, with the Benjamini-Hochberg FDR
   recomputed over the disease's Reactome pathways only (both ES signs, non-null p-values), keeping es > 0 and
   FDR < --fdr. A disease with no such pathway is left out.
2. Graph: exactly as Stage 3 builds it (stage3_build_graph.build_graph): nodes are all genes on those pathways
   (gene lists from merged_gene_sets.gmt), edges the curated Reactome functional interactions between them,
   directed and signed, unweighted.
3. Seeds: disease-associated genes in genetic_gene_scores_input.parquet.
Gene symbols are mapped to Ensembl gene IDs.

Writes --out-dir/disease_graph_catalogue/, with one <diseaseId>_disease_graph/ subdirectory per disease with
graph files (graph_weight.npz, graph_sign.npz, graph_gene_index.json, graph_metadata.json),
disease_seed_genes.tsv and disease_pathways.tsv.

Stage 5 can score a disease with --graph-dir and graph_catalogue_diseases.parquet, one summary row per disease.
"""

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "pipeline"))
import stage3_build_graph as bg
from stage0_schemas import DiseasePathwaysSchema, SeedGenesSchema
from stage1_ingest import (
    CoverageTracker,
    build_genes_table,
    build_interactions_table,
    build_symbol_index,
    read_target_table,
    resolve_symbol,
)

logger = logging.getLogger("graph_catalogue")

CATALOGUE_COLUMNS = ["diseaseId", "pathway", "source", "es", "nes", "pval", "geneset_size"]
CATALOGUE_SUBDIR = "disease_graph_catalogue"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="build_graph_catalogue.py",
        description="Build one graph per disease from the Open Targets disease-pathway GSEA catalogue"
        "(Reactome pathways only), with each disease's genetic-association scores as seeds.",
    )
    parser.add_argument("--catalogue-dir", type=Path, required=True,
                        help="Directory with the GSEA catalogue files, e.g. validation/Disease_graph_catalogue.")
    parser.add_argument("--target-dir", type=Path, required=True,
                        help="Open Targets target parquet directory, for symbol -> Ensembl, "
                        "e.g. data/opentargets/26.06/target.")
    parser.add_argument("--fi-file", type=Path, required=True,
                        help="Reactome functional-interaction file, unzipped, "
                        "e.g. data/reactome/FIsInGene_04142025_with_annotations.txt.")
    parser.add_argument("--out-dir", type=Path, default=None,
                        help=f"Writes <out-dir>/{CATALOGUE_SUBDIR}/. Default: --catalogue-dir.")
    parser.add_argument("--fdr", type=float, default=0.05, help="FDR cutoff (BH over Reactome). Default: 0.05.")
    parser.add_argument("--diseases", default=None,
                        help="Only these diseases, comma-separated Open Targets IDs, e.g. EFO_0003767,MONDO_0004975.")
    parser.add_argument("--n-diseases", type=int, default=None,
                        help="Only the first N included diseases (sorted by ID), e.g. for a quick test run.")
    return parser.parse_args(argv)


def benjamini_hochberg(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    adjusted = pvals[order] * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    fdr = np.empty(n)
    fdr[order] = np.minimum(adjusted, 1.0)
    return fdr


def select_reactome_pathways(catalogue, fdr_cutoff):
    """Reactome rows with BH recomputed per disease over its Reactome pathways; keep es > 0 and FDR < cutoff."""
    reactome = catalogue[(catalogue["source"] == "Reactome") & catalogue["pval"].notna()].copy()
    reactome["fdr_reactome"] = reactome.groupby("diseaseId")["pval"].transform(benjamini_hochberg)
    selected = reactome[(reactome["es"] > 0) & (reactome["fdr_reactome"] < fdr_cutoff)]
    return reactome["diseaseId"].nunique(), selected.reset_index(drop=True)


def read_reactome_gene_sets(gmt_path, index):
    """Reactome lines of the catalogue's GMT as (set_id = pathway name as in the catalogue, gene_id)."""
    coverage = CoverageTracker(source="catalogue_gmt_reactome")
    rows = []
    for line in gmt_path.read_text().splitlines():
        fields = line.split("\t")
        if not fields[0].endswith("[Reactome]"):
            continue
        genes = set()
        for symbol in fields[2:]:
            if not symbol:
                continue
            gene_id, _ = resolve_symbol(symbol, index)
            coverage.record(symbol, gene_id)
            if gene_id is not None:
                genes.add(gene_id)
        for gene_id in sorted(genes):
            rows.append((fields[0], fields[1], gene_id))
    gene_sets = pd.DataFrame(rows, columns=["set_id", "reactome_id", "gene_id"])
    return gene_sets, coverage.to_entry()


def read_seeds(scores_path, index, disease_ids):
    """Per disease, the genetic-association score of each gene mapped to Ensembl (max if symbols collide)."""
    scores = pd.read_parquet(scores_path, columns=["diseaseId", "gene", "score"])
    scores = scores[scores["diseaseId"].isin(disease_ids) & (scores["score"] > 0)]
    symbol_to_gene = {}
    for symbol in scores["gene"].unique():
        symbol_to_gene[symbol] = resolve_symbol(symbol, index)[0]
    scores = scores.assign(gene_id=scores["gene"].map(symbol_to_gene))
    n_unmapped = scores["gene_id"].isna().groupby(scores["diseaseId"]).sum()
    seeds = scores.dropna(subset=["gene_id"]).groupby(["diseaseId", "gene_id"], as_index=False)["score"].max()
    return seeds, n_unmapped


def build_disease_graph(pathways, gene_sets, interactions_by_gene_a):
    """Stage 3's graph on the disease's pathways. Interactions are pre-filtered to the node set, which leaves
    Stage 3's result unchanged (it drops out-of-scope rows itself) but keeps each call fast."""
    disease_gene_sets = gene_sets[gene_sets["set_id"].isin(pathways)]
    nodes = set(disease_gene_sets["gene_id"])
    candidate_rows = interactions_by_gene_a.loc[interactions_by_gene_a.index.isin(nodes)]
    in_scope = candidate_rows[candidate_rows["gene_b"].isin(nodes)].reset_index()
    graph, _ = bg.build_graph(disease_gene_sets, pd.DataFrame({"set_id": sorted(pathways)}), in_scope)
    return graph


def write_disease_dir(disease_dir, graph, seeds_df, pathways_df, metadata):
    """Stage 3's four graph files plus the seeds and pathways, so Stage 5 runs with --graph-dir disease_dir."""
    disease_dir.mkdir(parents=True, exist_ok=True)
    weight, sign, gene_index = bg.graph_to_sparse(graph)
    sp.save_npz(disease_dir / "graph_weight.npz", weight)
    sp.save_npz(disease_dir / "graph_sign.npz", sign)
    (disease_dir / "graph_gene_index.json").write_text(json.dumps(gene_index, indent=2))
    (disease_dir / "graph_metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True))
    SeedGenesSchema.validate(seeds_df).to_csv(disease_dir / "disease_seed_genes.tsv", sep="\t", index=False)
    DiseasePathwaysSchema.validate(pathways_df).to_csv(disease_dir / "disease_pathways.tsv", sep="\t", index=False)


def run(args):
    catalogue_out = (args.out_dir or args.catalogue_dir) / CATALOGUE_SUBDIR
    catalogue_out.mkdir(parents=True, exist_ok=True)

    logger.info("Reading the GSEA catalogue")
    catalogue = pd.read_parquet(args.catalogue_dir / "genetic_disease_pathway_catalogue.parquet",
                                columns=CATALOGUE_COLUMNS, filters=[("source", "==", "Reactome")])
    if args.diseases:
        requested = {disease.strip() for disease in args.diseases.split(",") if disease.strip()}
        unknown = sorted(requested - set(catalogue["diseaseId"]))
        if unknown:
            logger.warning("Not in the GSEA catalogue (not tested, or not an Open Targets 26.03 ID): %s",
                           ", ".join(unknown))
        catalogue = catalogue[catalogue["diseaseId"].isin(requested)]
    n_tested, selected = select_reactome_pathways(catalogue, args.fdr)
    disease_ids = sorted(selected["diseaseId"].unique())
    logger.info("Diseases: %d tested by GSEA, %d with a Reactome pathway at es > 0 and FDR < %g, %d left out.",
                n_tested, len(disease_ids), args.fdr, n_tested - len(disease_ids))
    if args.n_diseases is not None:
        if args.n_diseases < 1:
            raise SystemExit(f"--n-diseases must be at least 1 (got {args.n_diseases}).")
        disease_ids = disease_ids[:args.n_diseases]
        selected = selected[selected["diseaseId"].isin(disease_ids)]
        logger.info("--n-diseases: building the first %d.", len(disease_ids))

    logger.info("Mapping gene symbols to Ensembl")
    genes_df = build_genes_table(read_target_table(sorted(args.target_dir.glob("*.parquet"))))
    index = build_symbol_index(genes_df)
    gene_sets, gmt_coverage = read_reactome_gene_sets(args.catalogue_dir / "merged_gene_sets.gmt", index)
    logger.info("Reactome gene sets: %d, %.2f%% of their symbols mapped to Ensembl.",
                gene_sets["set_id"].nunique(), gmt_coverage["percent_mapped"])
    reactome_ids = gene_sets.drop_duplicates("set_id").set_index("set_id")["reactome_id"]
    missing_sets = set(selected["pathway"]) - set(reactome_ids.index)
    if missing_sets:
        logger.warning("%d selected pathway(s) not in the GMT, ignored: %s", len(missing_sets), sorted(missing_sets)[:5])
        selected = selected[~selected["pathway"].isin(missing_sets)]

    fi_coverage = CoverageTracker(source="reactome_fi")
    interactions = build_interactions_table(args.fi_file, index, "reactome", fi_coverage, curated_only=True,
                                            min_score=None)
    logger.info("Curated functional interactions: %d rows, %.2f%% of symbols mapped.",
                len(interactions), fi_coverage.to_entry()["percent_mapped"])
    interactions_by_gene_a = interactions.set_index("gene_a")

    seeds, n_unmapped_seeds = read_seeds(args.catalogue_dir / "genetic_gene_scores_input.parquet", index, disease_ids)
    seeds_by_disease = dict(tuple(seeds.groupby("diseaseId")))
    provenance = {
        "gsea_catalogue": str(args.catalogue_dir),
        "target_dir": str(args.target_dir),
        "fi_file": str(args.fi_file),
        "fdr_cutoff": args.fdr,
    }

    summaries = []
    for n, (disease_id, disease_pathways) in enumerate(selected.groupby("diseaseId"), start=1):
        graph = build_disease_graph(set(disease_pathways["pathway"]), gene_sets, interactions_by_gene_a)
        seeds_df = seeds_by_disease.get(disease_id, pd.DataFrame(columns=["gene_id", "score"]))
        seeds_df = seeds_df[["gene_id", "score"]].sort_values(["score", "gene_id"], ascending=[False, True])
        pathways_df = pd.DataFrame({
            "set_id": disease_pathways["pathway"].map(reactome_ids).to_numpy(),
            "gene_set": disease_pathways["pathway"].str.replace(r" \[Reactome\]$", "", regex=True).to_numpy(),
            "source_db": "reactome",
            "pval": disease_pathways["pval"].to_numpy(),
            "fdr": disease_pathways["fdr_reactome"].to_numpy(),
        }).sort_values(["fdr", "pval"]).drop_duplicates("set_id")
        summary = {
            "diseaseId": disease_id,
            "n_pathways": len(pathways_df),
            "n_nodes": graph.number_of_nodes(),
            "n_edges": graph.number_of_edges(),
            "n_isolated_nodes": sum(1 for gene in graph.nodes if graph.degree(gene) == 0),
            "n_seeds": len(seeds_df),
            "n_seeds_in_graph": int(seeds_df["gene_id"].isin(set(graph.nodes)).sum()),
            "n_seed_symbols_unmapped": int(n_unmapped_seeds.get(disease_id, 0)),
        }
        write_disease_dir(catalogue_out / f"{disease_id}_disease_graph", graph, seeds_df, pathways_df,
                          {**summary, **provenance})
        summaries.append(summary)
        if n % 250 == 0:
            logger.info("  %d/%d diseases", n, len(disease_ids))

    diseases_df = pd.DataFrame(summaries)
    diseases_df.to_parquet(catalogue_out / "graph_catalogue_diseases.parquet", index=False)

    no_seed = int((diseases_df["n_seeds_in_graph"] == 0).sum())
    no_edge = int((diseases_df["n_edges"] == 0).sum())
    logger.info("Catalogue: %d diseases; median %d nodes, %d edges, %d seeds. "
                "%d disease(s) without a seed in the graph, %d without an edge.",
                len(diseases_df), diseases_df["n_nodes"].median(), diseases_df["n_edges"].median(),
                diseases_df["n_seeds_in_graph"].median(), no_seed, no_edge)
    logger.info("Written to %s", catalogue_out)
    return catalogue_out


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    run(parse_args(argv))


if __name__ == "__main__":
    main()
