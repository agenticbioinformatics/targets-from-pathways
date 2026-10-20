"""Stage 4 - genetic-evidence weights.

Scores each graph gene with its Open Targets evidence for the disease (only the datatypes in --datatypes)
and re-weights the existing graph edges from the two endpoint scores. Writes the same four graph files
as Stage 3 (so Stage 5 reads either the same way) plus gene_weights.tsv for inspection.
graph_weight.npz is replaced; use a different --out-dir to keep Stage 3's structural weights.
Genes in Stage 2's benchmark_holdout_genes.tsv score 0.0: their own Open Targets evidence is what a benchmark
is meant to recover, so it must not weight their edges or Stage 6's composite.
"""

import argparse
import json
import logging
from pathlib import Path

import pandas as pd
import scipy.sparse as sp

from stage0_schemas import (
    GeneWeightsSchema,
    HoldoutGenesSchema,
    OTAssociationsSchema,
    fail,
    find_artifact,
    load_graph,
    load_manifest,
)

logger = logging.getLogger("genetic_evidence_weights")

# An inclusion list. affected_pathway must stay out: Open Targets files its Reactome evidence under it,
# and weighting a Reactome-built graph with Reactome evidence would be circular.
DEFAULT_DATATYPES = ["genetic_association", "known_drug"]


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="stage4_genetic_evidence_weights.py",
        description="Stage 4 - genetic-evidence weighting. Scores each graph gene from non-pathway-derived "
        "Open Targets datatypes (--datatypes) and rewrites graph_weight.npz with edge weights derived from "
        "those scores; extends graph_metadata.json and writes gene_weights.tsv.",
    )
    parser.add_argument("--manifest", required=True, type=Path, help="Path to Stage 1's manifest.json.")
    parser.add_argument("--graph-dir", type=Path, default=None, help="Stage 3 graph directory. Default: next to --manifest.")
    parser.add_argument("--out-dir", type=Path, default=None, help="Output directory. Default: --graph-dir.")
    parser.add_argument(
        "--datatypes", default=",".join(DEFAULT_DATATYPES),
        help=f"Comma-separated OT datatype_id values counted as evidence. Default: {','.join(DEFAULT_DATATYPES)}.",
    )
    parser.add_argument(
        "--holdout-genes", type=Path, default=None,
        help="Stage 2's benchmark_holdout_genes.tsv; these genes' evidence is ignored. Default: next to --manifest "
        "(skipped with a warning if absent).",
    )
    parser.add_argument("--edge-weight-mode", choices=["avg", "product"], default="avg",
                        help="How to combine the two endpoint gene scores into an edge weight.")
    return parser.parse_args(argv)


def compute_genetic_evidence_scores(ot_disease_subset, datatypes, gene_index):
    """One score per graph gene: its strongest admitted evidence row, or 0.0 if it has none."""
    admitted = ot_disease_subset[ot_disease_subset["datatype_id"].isin(datatypes)]
    if admitted.empty:
        logger.warning("No ot_disease_subset rows match --datatypes %s - every gene scores 0.0.", datatypes)
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


def compute_edge_weights(existing_weight, node_scores, mode):
    """Weight each existing edge from its two endpoint scores. No edge is added or removed."""
    edges = existing_weight.tocoo()
    source_scores = node_scores[edges.row]
    target_scores = node_scores[edges.col]
    if mode == "avg":
        new_weights = (source_scores + target_scores) / 2.0
    elif mode == "product":
        new_weights = source_scores * target_scores
    else:
        raise ValueError(f"unknown edge-weight-mode {mode!r}")
    return sp.csr_matrix((new_weights, (edges.row, edges.col)), shape=existing_weight.shape)


def run(args):
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

    datatypes = []
    for datatype in args.datatypes.split(","):
        if datatype.strip():
            datatypes.append(datatype.strip())
    if not datatypes:
        fail("--datatypes resolved to an empty list.")

    scores = compute_genetic_evidence_scores(ot_disease_subset, datatypes, gene_index)
    node_scores = scores.to_numpy()
    edge_weight = compute_edge_weights(weight, node_scores, args.edge_weight_mode)
    n_genes_with_evidence = int((node_scores > 0).sum())
    logger.info("Genetic evidence: %d/%d genes scored, %d edges weighted (%s).",
                n_genes_with_evidence, len(gene_index), edge_weight.nnz, args.edge_weight_mode)

    out_dir = args.out_dir or graph_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    sp.save_npz(out_dir / "graph_weight.npz", edge_weight)
    sp.save_npz(out_dir / "graph_sign.npz", sign)
    (out_dir / "graph_gene_index.json").write_text(json.dumps(gene_index, indent=2))

    # Stage 5 reads the gene scores from here, not from gene_weights.tsv.
    metadata["genetic_evidence_score"] = {}
    for gene_id, score in zip(gene_index, node_scores):
        metadata["genetic_evidence_score"][gene_id] = float(score)
    metadata["stage4_genetic_evidence"] = {
        "datatypes": datatypes,
        "edge_weight_mode": args.edge_weight_mode,
        "n_genes_with_evidence": n_genes_with_evidence,
        "n_genes_total": len(gene_index),
        "n_edges_weighted": int(edge_weight.nnz),
        "holdout_genes": sorted(holdout_gene_ids),
    }
    (out_dir / "graph_metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True))

    weights_df = pd.DataFrame({"gene_id": gene_index, "genetic_evidence_score": node_scores})
    weights_df = weights_df.sort_values("genetic_evidence_score", ascending=False).reset_index(drop=True)
    weights_df = GeneWeightsSchema.validate(weights_df)
    rounded = []
    for score in weights_df["genetic_evidence_score"]:
        rounded.append(float(f"{score:.4g}"))
    weights_df["genetic_evidence_score"] = rounded
    weights_df.to_csv(out_dir / "gene_weights.tsv", sep="\t", index=False)

    logger.info("Stage 4 complete: %s", out_dir)
    return out_dir


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    run(parse_args(argv))


if __name__ == "__main__":
    main()
