"""Stage 4 - optional user weights, and genetic-evidence gene scores.

Weights Stage 3's unweighted graph with the user's own node and/or edge weights, both optional:
- --node-weights: a TSV with gene_id, weight. Each edge gets the mean (or, with --node-weight-mode product,
  the product) of its two genes' weights; genes not in the file weigh --default-node-weight.
- --edge-weights: a TSV with gene_a, gene_b, weight. A row weights the graph edge(s) between the two genes, in
  whichever directions the graph has; edges not in the file weigh --default-edge-weight.
With both, an edge's weight is the product of the two. Rows for genes or edges that are not in the graph are
skipped (logged): Stage 4 never adds or removes an edge. Without either file every edge keeps weight 1.
Writes the same four graph files as Stage 3 (graph_weight.npz is replaced), so Stage 5 reads them the same way.
The weights are always applied to Stage 3's structure (every edge reset to 1 first), so re-running Stage 4
with other files, or with none, does not compound.

Also scores each graph gene with its Open Targets genetic_association evidence for the disease (no other
datatype), for Stage 6's annotation: added to graph_metadata.json and written to genetic_evidence_scores.tsv.
These scores are not used as weights. Genes in Stage 2's benchmark_holdout_genes.tsv score 0.0: their own Open Targets
evidence is what a benchmark is meant to recover, so it must not show as their genetic_evidence_score.
"""

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from pandera.errors import SchemaError

from stage0_schemas import (
    GeneticEvidenceScoresSchema,
    HoldoutGenesSchema,
    OTAssociationsSchema,
    UserEdgeWeightsSchema,
    UserNodeWeightsSchema,
    fail,
    find_artifact,
    load_graph,
    load_manifest,
)

logger = logging.getLogger("weights")

# The only Open Targets datatype counted as evidence. affected_pathway would be circular (Open Targets files its
# Reactome evidence under it, and the graph is built from Reactome), and known_drug would favour genes that are
# already drug targets for the disease.
EVIDENCE_DATATYPES = ["genetic_association"]


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="stage4_weights.py",
        description="Stage 4 - optional user weights. Weights the Stage 3 graph with the user's own node "
        "and/or edge weights (--node-weights / --edge-weights, both optional; without them every edge keeps "
        "weight 1) and rewrites the graph files. Also scores each graph gene from its Open Targets "
        "genetic_association evidence for Stage 6 (not used as a weight); extends graph_metadata.json and "
        "writes genetic_evidence_scores.tsv.",
    )
    parser.add_argument("--manifest", required=True, type=Path, help="Path to Stage 1's manifest.json.")
    parser.add_argument("--graph-dir", type=Path, default=None, help="Stage 3 graph directory. Default: next to --manifest.")
    parser.add_argument("--out-dir", type=Path, default=None, help="Output directory. Default: --graph-dir.")
    parser.add_argument("--node-weights", type=Path, default=None,
                        help="Optional TSV with columns gene_id (Ensembl), weight (> 0).")
    parser.add_argument("--edge-weights", type=Path, default=None,
                        help="Optional TSV with columns gene_a, gene_b (Ensembl), weight (> 0).")
    parser.add_argument("--node-weight-mode", choices=["mean", "product"], default="mean",
                        help="How an edge combines its two genes' node weights. Default: mean.")
    parser.add_argument("--default-node-weight", type=float, default=1.0,
                        help="Weight of a gene missing from --node-weights. Default: 1.")
    parser.add_argument("--default-edge-weight", type=float, default=1.0,
                        help="Weight of an edge missing from --edge-weights. Default: 1.")
    parser.add_argument(
        "--holdout-genes", type=Path, default=None,
        help="Stage 2's benchmark_holdout_genes.tsv; these genes' evidence is ignored. Default: next to --manifest "
        "(skipped with a warning if absent).",
    )
    return parser.parse_args(argv)


def compute_genetic_evidence_scores(ot_disease_subset, datatypes, gene_index):
    """One score per graph gene: its strongest admitted evidence row, or 0.0 if it has none."""
    admitted = ot_disease_subset[ot_disease_subset["datatype_id"].isin(datatypes)]
    if admitted.empty:
        logger.warning("No ot_disease_subset rows of datatype %s - every gene scores 0.0.", datatypes)
    best_score_per_gene = admitted.groupby("gene_id")["score"].max()

    scores = pd.Series(0.0, index=gene_index, name="genetic_evidence_score")
    for gene_id, score in best_score_per_gene.items():
        if gene_id in scores.index:
            scores[gene_id] = score
    return scores


def read_holdout_genes(holdout_path, explicit):
    if not holdout_path.exists():
        if explicit:
            fail(f"--holdout-genes {holdout_path} does not exist.")
        logger.warning("No %s (Stage 2 from before holdout tracking?) - no gene's evidence is excluded.", holdout_path)
        return set()
    return set(HoldoutGenesSchema.validate(pd.read_csv(holdout_path, sep="\t"))["gene_id"])


def read_weights_file(path, schema, option):
    if not path.exists():
        fail(f"{option} {path} does not exist.")
    try:
        return schema.validate(pd.read_csv(path, sep="\t"))
    except SchemaError as error:
        fail(f"{option} {path} is not a valid weights file: {error}")


def node_weight_factors(node_weights_df, gene_index, edges, mode, default):
    """Per edge, the combination of its two genes' node weights. Returns (factors, n_used, n_not_in_graph)."""
    position = {gene: i for i, gene in enumerate(gene_index)}
    node_weight = np.full(len(gene_index), default)
    n_used = 0
    for gene_id, weight in zip(node_weights_df["gene_id"], node_weights_df["weight"]):
        if gene_id in position:
            node_weight[position[gene_id]] = weight
            n_used += 1
    source, target = node_weight[edges.row], node_weight[edges.col]
    if mode == "mean":
        factors = (source + target) / 2.0
    elif mode == "product":
        factors = source * target
    else:
        raise ValueError(f"unknown node-weight-mode {mode!r}")
    return factors, n_used, len(node_weights_df) - n_used


def edge_weight_factors(edge_weights_df, gene_index, edges, default):
    """Per edge, the user's weight for its gene pair (either direction). Returns (factors, n_used, n_not_in_graph)."""
    pairs = pd.DataFrame({"a": np.minimum(edges.row, edges.col), "b": np.maximum(edges.row, edges.col)})
    position = {gene: i for i, gene in enumerate(gene_index)}
    user = edge_weights_df.assign(i=edge_weights_df["gene_a"].map(position), j=edge_weights_df["gene_b"].map(position))
    in_graph_genes = user["i"].notna() & user["j"].notna()
    user = user[in_graph_genes].astype({"i": int, "j": int})
    user = user.assign(a=np.minimum(user["i"], user["j"]), b=np.maximum(user["i"], user["j"]))
    duplicated = user.duplicated(subset=["a", "b"], keep=False)
    if duplicated.any():
        rows = user.loc[duplicated, ["gene_a", "gene_b"]].head(5).values.tolist()
        fail(f"--edge-weights lists the same gene pair more than once (in either order), e.g. {rows}.")
    matched = pairs.merge(user[["a", "b", "weight"]], on=["a", "b"], how="left")["weight"]
    n_used = int(user.set_index(["a", "b"]).index.isin(pairs.set_index(["a", "b"]).index).sum())
    return matched.fillna(default).to_numpy(), n_used, len(edge_weights_df) - n_used


def apply_user_weights(weight, gene_index, node_weights_df, edge_weights_df, node_weight_mode,
                       default_node_weight, default_edge_weight):
    """Weight every existing edge from the user's files, starting from 1 (Stage 3's unweighted graph).
    No edge is added or removed. Returns (new weight matrix, stats)."""
    edges = weight.tocoo()
    new_weights = np.ones(edges.nnz)
    stats = {"node_weight_mode": None, "n_node_weights_used": 0, "n_node_weights_not_in_graph": 0,
             "n_edge_weights_used": 0, "n_edge_weights_not_in_graph": 0}
    if node_weights_df is not None:
        factors, n_used, n_skipped = node_weight_factors(node_weights_df, gene_index, edges, node_weight_mode,
                                                         default_node_weight)
        new_weights *= factors
        stats.update(node_weight_mode=node_weight_mode, n_node_weights_used=n_used,
                     n_node_weights_not_in_graph=n_skipped)
    if edge_weights_df is not None:
        factors, n_used, n_skipped = edge_weight_factors(edge_weights_df, gene_index, edges, default_edge_weight)
        new_weights *= factors
        stats.update(n_edge_weights_used=n_used, n_edge_weights_not_in_graph=n_skipped)
    return sp.csr_matrix((new_weights, (edges.row, edges.col)), shape=weight.shape), stats


def run(args):
    for option, value in [("--default-node-weight", args.default_node_weight),
                          ("--default-edge-weight", args.default_edge_weight)]:
        if not (np.isfinite(value) and value > 0):
            fail(f"{option} must be > 0 (got {value}).")
    manifest = load_manifest(args.manifest)
    ot_disease_subset = OTAssociationsSchema.validate(
        pd.read_parquet(find_artifact(manifest, args.manifest, "ot_disease_subset.parquet"))
    )
    graph_dir = args.graph_dir or args.manifest.parent
    weight, sign, gene_index, metadata = load_graph(graph_dir)
    holdout_gene_ids = read_holdout_genes(args.holdout_genes or args.manifest.parent / "benchmark_holdout_genes.tsv",
                                          args.holdout_genes is not None)
    if holdout_gene_ids:
        logger.info("Benchmark holdout: ignoring the evidence of %d gene(s).", len(holdout_gene_ids))
        ot_disease_subset = ot_disease_subset[~ot_disease_subset["gene_id"].isin(holdout_gene_ids)]

    scores = compute_genetic_evidence_scores(ot_disease_subset, EVIDENCE_DATATYPES, gene_index)
    node_scores = scores.to_numpy()
    n_genes_with_evidence = int((node_scores > 0).sum())
    logger.info("Genetic evidence: %d/%d genes scored.", n_genes_with_evidence, len(gene_index))

    node_weights_df = edge_weights_df = None
    if args.node_weights is not None:
        node_weights_df = read_weights_file(args.node_weights, UserNodeWeightsSchema, "--node-weights")
    if args.edge_weights is not None:
        edge_weights_df = read_weights_file(args.edge_weights, UserEdgeWeightsSchema, "--edge-weights")
    edge_weight, weight_stats = apply_user_weights(weight, gene_index, node_weights_df, edge_weights_df,
                                                   args.node_weight_mode, args.default_node_weight,
                                                   args.default_edge_weight)
    if node_weights_df is not None:
        logger.info("Node weights: %d genes weighted, %d not in the graph (skipped); edges combine them by %s.",
                    weight_stats["n_node_weights_used"], weight_stats["n_node_weights_not_in_graph"],
                    args.node_weight_mode)
    if edge_weights_df is not None:
        logger.info("Edge weights: %d gene pairs weighted, %d not graph edges (skipped).",
                    weight_stats["n_edge_weights_used"], weight_stats["n_edge_weights_not_in_graph"])
    if node_weights_df is None and edge_weights_df is None:
        logger.info("No --node-weights / --edge-weights: the graph stays unweighted.")

    out_dir = args.out_dir or graph_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    sp.save_npz(out_dir / "graph_weight.npz", edge_weight)
    sp.save_npz(out_dir / "graph_sign.npz", sign)
    (out_dir / "graph_gene_index.json").write_text(json.dumps(gene_index, indent=2))

    # Stage 6 reads the gene scores from here, not from genetic_evidence_scores.tsv.
    metadata["genetic_evidence_score"] = {}
    for gene_id, score in zip(gene_index, node_scores):
        metadata["genetic_evidence_score"][gene_id] = float(score)
    metadata["stage4_genetic_evidence"] = {
        "datatypes": EVIDENCE_DATATYPES,
        "n_genes_with_evidence": n_genes_with_evidence,
        "n_genes_total": len(gene_index),
        "holdout_genes": sorted(holdout_gene_ids),
    }
    weights_used = node_weights_df is not None or edge_weights_df is not None
    metadata["stage4_user_weights"] = {
        "weighted": weights_used,
        "node_weights_file": str(args.node_weights) if args.node_weights is not None else None,
        "edge_weights_file": str(args.edge_weights) if args.edge_weights is not None else None,
        "default_node_weight": args.default_node_weight,
        "default_edge_weight": args.default_edge_weight,
        **weight_stats,
    }
    (out_dir / "graph_metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True))

    scores_df = pd.DataFrame({"gene_id": gene_index, "genetic_evidence_score": node_scores})
    scores_df = scores_df.sort_values("genetic_evidence_score", ascending=False).reset_index(drop=True)
    scores_df = GeneticEvidenceScoresSchema.validate(scores_df)
    rounded = []
    for score in scores_df["genetic_evidence_score"]:
        rounded.append(float(f"{score:.4g}"))
    scores_df["genetic_evidence_score"] = rounded
    scores_df.to_csv(out_dir / "genetic_evidence_scores.tsv", sep="\t", index=False)

    logger.info("Stage 4 complete: %s", out_dir)
    return out_dir


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    run(parse_args(argv))


if __name__ == "__main__":
    main()
