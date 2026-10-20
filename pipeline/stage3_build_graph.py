"""Stage 3 - pathway-based gene-gene graph.

Builds a directed, signed graph whose nodes are all genes of the disease-relevant pathways (Stage 2) plus
all genes of every pathway that contains a seed gene (Stage 2's disease_seed_genes.tsv, the disease's
genetic-association genes). Edges come from pathway co-membership and from Reactome functional
interactions between those genes. Writes graph_weight.npz, graph_sign.npz (sparse matrices sharing the
row/column order in graph_gene_index.json) and graph_metadata.json.
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
        description="Stage 3 - pathway-based gene-gene graph. Builds a directed, signed graph from the union "
        "of disease-relevant and seed-gene-containing pathways (co-membership edges + Reactome "
        "functional-interaction edges). Reads gene_sets/interactions via --manifest, disease_pathways.tsv and "
        "disease_seed_genes.tsv; writes graph_weight.npz / graph_sign.npz / graph_gene_index.json / "
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


def add_comembership_edges(graph, gene_sets_df, pathway_ids):
    """Connect every pair of genes in each pathway, in both directions, with sign 0 (unknown effect).
    Each pathway adds 1 / (size - 1) to a pair, so a few large pathways cannot dominate the graph.
    Returns the number of distinct gene pairs."""
    selected = gene_sets_df[gene_sets_df["set_id"].isin(pathway_ids)]
    graph.add_nodes_from(selected["gene_id"].unique())

    gene_pairs = set()
    for _, pathway in selected.groupby("set_id"):
        genes = sorted(pathway["gene_id"].unique())
        if len(genes) < 2:
            continue
        contribution = 1.0 / (len(genes) - 1)
        source_db = pathway["source_db"].iloc[0]
        for i in range(len(genes)):
            for j in range(i + 1, len(genes)):
                gene_pairs.add((genes[i], genes[j]))
                for source, target in [(genes[i], genes[j]), (genes[j], genes[i])]:
                    if graph.has_edge(source, target):
                        graph[source][target]["comembership_weight"] += contribution
                    else:
                        graph.add_edge(source, target, comembership_weight=contribution, sign=0,
                                       confidence=None, evidence_type=None, source_db=source_db)
    return len(gene_pairs)


def add_interaction_edges(graph, interactions_df, pathway_genes):
    """Add interactions whose two genes are both already in the graph. Adding all interactions would turn
    the pathway graph into the whole interactome. Returns (n_added, n_out_of_scope, n_duplicates)."""
    # A DiGraph holds one edge per ordered pair: keep the highest-confidence row for each pair.
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
        if not graph.has_edge(row.gene_a, row.gene_b):
            graph.add_edge(row.gene_a, row.gene_b, comembership_weight=0.0)
        # The interaction sign is more specific than co-membership's sign 0, so it overrides it.
        edge = graph[row.gene_a][row.gene_b]
        edge["sign"] = int(row.sign)
        edge["confidence"] = confidence
        edge["evidence_type"] = row.evidence_type
        edge["source_db"] = row.source_db
    return int(in_scope.sum()), n_out_of_scope, n_duplicates


def build_graph(gene_sets_df, disease_pathways_df, interactions_df, seed_gene_ids):
    disease_pathway_ids = set(disease_pathways_df["set_id"])
    seed_pathway_ids = set(gene_sets_df.loc[gene_sets_df["gene_id"].isin(seed_gene_ids), "set_id"])
    if not seed_pathway_ids:
        logger.warning("None of the %d seed genes is in a pathway of gene_sets.parquet.", len(seed_gene_ids))
    pathway_ids = disease_pathway_ids | seed_pathway_ids

    graph = nx.DiGraph()
    n_gene_pairs = add_comembership_edges(graph, gene_sets_df, pathway_ids)
    n_added, n_out_of_scope, n_duplicates = add_interaction_edges(graph, interactions_df, set(graph.nodes))

    stats = {
        "n_pathways_union": len(pathway_ids),
        "n_disease_pathways": len(disease_pathway_ids),
        "n_seed_pathways": len(seed_pathway_ids),
        "n_nodes": graph.number_of_nodes(),
        "n_edges": graph.number_of_edges(),
        "n_comembership_unordered_pairs": n_gene_pairs,
        "n_interaction_edges_added": n_added,
        "n_interaction_edges_out_of_scope": n_out_of_scope,
        "n_interaction_duplicate_pairs_collapsed": n_duplicates,
    }
    return graph, stats


def graph_to_sparse(graph):
    """Return (weight, sign, gene_index). weight = co-membership weight + interaction confidence.
    A stored sign of 0 means either no edge or an unsigned edge, so read sign only where weight > 0.
    Edge direction is the asymmetry between weight[i, j] and weight[j, i]."""
    gene_index = sorted(graph.nodes)
    for gene in gene_index:
        if not ENSEMBL_GENE_ID_PATTERN.match(gene):
            fail(f"Graph contains a non-Ensembl gene ID: {gene}")

    for _, _, edge in graph.edges(data=True):
        confidence = edge["confidence"]
        if confidence is None:
            confidence = 0.0
        edge["matrix_weight"] = edge["comembership_weight"] + confidence

    weight = nx.to_scipy_sparse_array(graph, nodelist=gene_index, weight="matrix_weight", format="csr", dtype="float64")
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

    graph, stats = build_graph(gene_sets_df, disease_pathways_df, interactions_df, seed_gene_ids)
    logger.info("Graph: %d nodes, %d edges (%d pathways: %d disease, %d with seed genes; "
                "%d interaction edges added, %d out of scope).",
                stats["n_nodes"], stats["n_edges"], stats["n_pathways_union"], stats["n_disease_pathways"],
                stats["n_seed_pathways"], stats["n_interaction_edges_added"], stats["n_interaction_edges_out_of_scope"])
    weight, sign, gene_index = graph_to_sparse(graph)
    # Seeds that are in no size-capped pathway are not graph nodes, so they cannot seed Stage 5's propagation.
    n_seeds_in_graph = len(seed_gene_ids & set(gene_index))
    logger.info("Seed genes in the graph: %d/%d.", n_seeds_in_graph, len(seed_gene_ids))

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
