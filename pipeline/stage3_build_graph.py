"""Stage 3 - pathway-based gene-gene graph.

Builds a directed, signed, unweighted graph whose nodes are all genes of the disease-relevant pathways
detected by GSEA (Stage 2's disease_pathways.tsv) and whose edges are the Reactome functional interactions
between those genes. Every edge has weight 1. Writes graph_weight.npz, graph_sign.npz (sparse matrices
sharing the row/column order in graph_gene_index.json) and graph_metadata.json.
"""

import argparse
import json
import logging
from datetime import datetime, timezone
from pathlib import Path

import networkx as nx
import pandas as pd
import scipy.sparse as sp

from stage0_schemas import (
    ENSEMBL_GENE_ID_PATTERN,
    DiseasePathwaysSchema,
    GeneSetsSchema,
    InteractionsSchema,
    SeedGenesSchema,
    fail,
    find_artifact,
    load_manifest,
)

logger = logging.getLogger("build_graph")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="stage3_build_graph.py",
        description="Stage 3 - pathway-based gene-gene graph. Builds a directed, signed, unweighted graph of the "
        "Reactome functional interactions between the genes of the disease-relevant pathways found by GSEA. "
        "Reads gene_sets/interactions via --manifest, disease_pathways.tsv and disease_seed_genes.tsv (only to "
        "check the seeds are in the graph); writes graph_weight.npz / graph_sign.npz / graph_gene_index.json / "
        "graph_metadata.json.",
    )
    parser.add_argument("--manifest", required=True, type=Path, help="Path to Stage 1's manifest.json.")
    parser.add_argument(
        "--disease-pathways", type=Path, default=None,
        help="Path to Stage 2's disease_pathways.tsv. Default: next to --manifest.",
    )
    parser.add_argument(
        "--seed-genes", type=Path, default=None,
        help="Path to Stage 2's disease_seed_genes.tsv. Default: next to --manifest.",
    )
    parser.add_argument("--out-dir", type=Path, default=None, help="Output directory. Default: --manifest's directory.")
    return parser.parse_args(argv)


def add_pathway_genes(graph, gene_sets_df, pathway_ids):
    """Add every gene of the pathways as a node (also genes without any interaction)."""
    selected = gene_sets_df[gene_sets_df["set_id"].isin(pathway_ids)]
    graph.add_nodes_from(selected["gene_id"].unique())


def add_interaction_edges(graph, interactions_df, pathway_genes):
    """Add interactions whose two genes are both already in the graph. Adding all interactions would turn
    the pathway graph into the whole interactome. Returns (n_added, n_out_of_scope, n_duplicates), counted
    over directed edges (an undirected row counts as two)."""
    # Stage 1 stores an undirected interaction as one row; it connects the genes in both directions.
    undirected = interactions_df[~interactions_df["directed"]]
    reversed_rows = undirected.rename(columns={"gene_a": "gene_b", "gene_b": "gene_a"})
    interactions_df = pd.concat([interactions_df, reversed_rows], ignore_index=True)
    # A DiGraph holds one edge per ordered pair: keep the highest-confidence row for each pair (on a tie the
    # original row wins over a reversed undirected copy).
    by_confidence = interactions_df.assign(rank=interactions_df["confidence"].fillna(-1.0))
    by_confidence = by_confidence.sort_values("rank", ascending=False, kind="stable")
    unique_pairs = by_confidence.drop_duplicates(subset=["gene_a", "gene_b"], keep="first")
    n_duplicates = len(interactions_df) - len(unique_pairs)

    in_scope = unique_pairs["gene_a"].isin(pathway_genes) & unique_pairs["gene_b"].isin(pathway_genes)
    n_out_of_scope = int((~in_scope).sum())

    for row in unique_pairs[in_scope].itertuples(index=False):
        confidence = None
        if not pd.isna(row.confidence):
            confidence = float(row.confidence)
        graph.add_edge(row.gene_a, row.gene_b, sign=int(row.sign), confidence=confidence,
                       evidence_type=row.evidence_type, source_db=row.source_db)
    return int(in_scope.sum()), n_out_of_scope, n_duplicates


def build_graph(gene_sets_df, disease_pathways_df, interactions_df):
    disease_pathway_ids = set(disease_pathways_df["set_id"])
    graph = nx.DiGraph()
    add_pathway_genes(graph, gene_sets_df, disease_pathway_ids)
    n_added, n_out_of_scope, n_duplicates = add_interaction_edges(graph, interactions_df, set(graph.nodes))

    stats = {
        "n_disease_pathways": len(disease_pathway_ids),
        "n_nodes": graph.number_of_nodes(),
        "n_edges": graph.number_of_edges(),
        "n_interaction_edges_added": n_added,
        "n_interaction_edges_out_of_scope": n_out_of_scope,
        "n_interaction_duplicate_pairs_collapsed": n_duplicates,
    }
    return graph, stats


def graph_to_sparse(graph):
    """Return (weight, sign, gene_index). The graph is unweighted: weight is 1 on every edge and 0 elsewhere.
    A stored sign of 0 means either no edge or an unsigned edge, so read sign only where weight > 0.
    Edge direction is the asymmetry between weight[i, j] and weight[j, i]."""
    gene_index = sorted(graph.nodes)
    for gene in gene_index:
        if not ENSEMBL_GENE_ID_PATTERN.match(gene):
            fail(f"Graph contains a non-Ensembl gene ID: {gene}")

    weight = nx.to_scipy_sparse_array(graph, nodelist=gene_index, weight=None, format="csr", dtype="float64")
    sign = nx.to_scipy_sparse_array(graph, nodelist=gene_index, weight="sign", format="csr", dtype="int8")
    return weight, sign, gene_index


def run(args):
    manifest = load_manifest(args.manifest)
    gene_sets_df = GeneSetsSchema.validate(pd.read_parquet(find_artifact(manifest, args.manifest, "gene_sets.parquet")))
    interactions_df = InteractionsSchema.validate(
        pd.read_parquet(find_artifact(manifest, args.manifest, "interactions.parquet"))
    )
    disease_pathways_path = args.disease_pathways or args.manifest.parent / "disease_pathways.tsv"
    if not disease_pathways_path.exists():
        fail(f"--disease-pathways {disease_pathways_path} does not exist (run stage2_gsea_discovery.py first?).")
    disease_pathways_df = DiseasePathwaysSchema.validate(pd.read_csv(disease_pathways_path, sep="\t"))
    seed_genes_path = args.seed_genes or args.manifest.parent / "disease_seed_genes.tsv"
    if not seed_genes_path.exists():
        fail(f"--seed-genes {seed_genes_path} does not exist (run stage2_gsea_discovery.py first?).")
    seed_gene_ids = set(SeedGenesSchema.validate(pd.read_csv(seed_genes_path, sep="\t"))["gene_id"])

    if disease_pathways_df.empty:
        fail(f"{disease_pathways_path} lists no disease-relevant pathway (GSEA found none), so there is no "
             f"graph to build. See Stage 2's log for why.")
    graph, stats = build_graph(gene_sets_df, disease_pathways_df, interactions_df)
    logger.info("Graph: %d nodes (genes on the %d disease-relevant pathways), %d edges (%d interaction edges "
                "added, %d out of scope).", stats["n_nodes"], stats["n_disease_pathways"], stats["n_edges"],
                stats["n_interaction_edges_added"], stats["n_interaction_edges_out_of_scope"])
    weight, sign, gene_index = graph_to_sparse(graph)
    # Seeds that are on no disease pathway are not graph nodes, so they cannot seed Stage 5's propagation.
    n_seeds_in_graph = len(seed_gene_ids & set(gene_index))
    logger.info("Seed genes in the graph: %d/%d.", n_seeds_in_graph, len(seed_gene_ids))
    if n_seeds_in_graph == 0:
        fail(f"None of the {len(seed_gene_ids)} seed gene(s) in {seed_genes_path} is on a disease-relevant "
             f"pathway, so there is nothing to propagate from.")

    reactome_sources = []
    for source in manifest.sources:
        if source.db == "reactome":
            reactome_sources.append(source.model_dump(mode="json"))
    metadata = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "seed": manifest.seed,
        "reactome_sources": reactome_sources,
        "disease": manifest.disease.model_dump(mode="json"),
        "n_seed_genes_in_graph": n_seeds_in_graph,
        "n_seed_genes_total": len(seed_gene_ids),
    }
    metadata.update(stats)

    out_dir = args.out_dir or args.manifest.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    sp.save_npz(out_dir / "graph_weight.npz", weight)
    sp.save_npz(out_dir / "graph_sign.npz", sign)
    (out_dir / "graph_gene_index.json").write_text(json.dumps(gene_index, indent=2))
    metadata_path = out_dir / "graph_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True))
    logger.info("Stage 3 complete: %s", metadata_path)
    return metadata_path


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    run(parse_args(argv))


if __name__ == "__main__":
    main()
