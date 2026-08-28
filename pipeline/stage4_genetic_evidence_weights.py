"""Stage 4 — Genetic-evidence node/edge weighting.

Takes Stage 3's pathway-based gene-gene graph (its sparse serialization —
``graph_weight.npz``, ``graph_sign.npz``, ``graph_gene_index.json``,
``graph_metadata.json``, read via ``--graph-dir``) and Open Targets
association data (``ot_disease_subset.parquet``, read via ``--manifest``,
same as ``stage2_gsea_discovery.py``), and maps non-pathway-derived genetic
evidence onto it as node and edge weights — never a raw source file.

**Datatype restriction is load-bearing, not theoretical.** ``--datatypes``
(default ``genetic_association,known_drug``) is an *inclusion* list, so
everything else — ``affected_pathway`` and any literature/pathway-derived
datatype in particular — is excluded simply by omission. This matters
concretely, not just in principle: Open Targets 26.06's ``reactome``
evidence datasource is itself mapped to the ``affected_pathway`` datatype
(see ``stage1_ingest.py``'s ``_DATASOURCE_TO_DATATYPE``), so an unrestricted
(or aggregated/overall) association score would fold Reactome-derived
evidence into the very node/edge weights that then feed Stage 5's scoring
over a graph whose topology *is* Reactome pathway structure — validating
Reactome against itself. This is the same orthogonality discipline
``stage2_gsea_discovery.py`` applies to Stage 2's ranking signature (see its
module docstring); Stage 4 is the second of the two places it matters.

Where a gene has more than one admitted-datatype row (multiple datasources,
e.g. ``eva`` and ``clinical_precedence``, or both ``genetic_association``
and ``known_drug``), the max score wins — the same "strongest single piece
of evidence, not Open Targets' own cross-datasource aggregate" convention
``stage2_gsea_discovery.py.build_signature`` already uses, applied here across
datatypes as well as datasources.

Every node in Stage 3's graph gets a score: 0.0, not null, for a gene with
no matching evidence row, so downstream weighting never has to special-case
a missing value.

**Edge weights** are derived *only* from the two endpoint node scores
(``--edge-weight-mode``: ``avg`` — the default — or ``product``), computed
at exactly the ``(gene_a, gene_b)`` positions Stage 3's ``weight`` matrix
already has a real entry for — genetic evidence never invents a new edge,
only re-weights an existing one.

**Output is Stage 3's exact serialization format, not a parallel one** —
the same four filenames, same shapes, in the same directory layout — so
``stage5_score_candidates.py`` (Stage 5) can point ``--graph`` at either a Stage 3
directory or a Stage 4 directory and load it through the *same* code path,
with no branching to detect which stage produced it:

- ``graph_weight.npz`` — **replaced**. Stage 3 stores structural
  (comembership + interaction-confidence) weight here; once this script has
  run, the same file holds the genetic-evidence-derived edge weight instead
  (see the attribute-name comment below), at exactly the same ``(row,
  col)`` sparsity pattern as Stage 3's original — no edge is added or
  removed, only re-weighted. This means a real edge between two genes with
  no genetic evidence at all is now stored as an explicit ``0.0`` rather
  than a nonzero structural weight; that is still a present edge (the
  sparsity pattern hasn't changed), so anything checking edge existence by
  presence in the sparse structure — not by "is the stored value nonzero" —
  is unaffected. If you need Stage 3's original structural weight preserved
  alongside this one, pass a different ``--out-dir`` than ``--graph-dir``.
- ``graph_sign.npz`` / ``graph_gene_index.json`` — copied through
  unchanged.
- ``graph_metadata.json`` — Stage 3's metadata, extended (not replaced)
  with the ``genetic_evidence_score`` mapping and a ``stage4_genetic_evidence``
  stats block (``--datatypes``, ``--edge-weight-mode``, coverage counts).
- ``gene_weights.tsv`` — a separate, human-facing artifact (schema-validated
  against ``schemas.GeneWeightsSchema``, sorted by score descending) for
  eyeballing. Not part of the "same format as Stage 3" contract above —
  Stage 5 never reads this file, only ``graph_metadata.json``'s
  ``genetic_evidence_score`` mapping.
"""

# ==========================================================================
# Node/edge weight attribute names (the slots this script writes into
# Stage 3's serialization format — see the module docstring's "Output"
# section above for the full read/write contract):
#
#   NODE weight -> "genetic_evidence_score", a {gene_id: score} mapping
#                  under that exact key in graph_metadata.json.
#   EDGE weight -> graph_weight.npz's stored values — Stage 3's "weight"
#                  slot, reused rather than a separate attribute name, so
#                  Stage 5 always reads edge weight from one place
#                  regardless of which stage last wrote it.
# ==========================================================================

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from stage0_schemas import GeneWeightsSchema, Manifest, OTAssociationsSchema

logger = logging.getLogger("genetic_evidence_weights")

# --datatypes' default is an *inclusion* list, not an exclusion list —
# everything not named here (affected_pathway, literature, somatic_mutation,
# rna_expression, animal_model, ...) is excluded simply by omission.
# affected_pathway is the one that matters most to keep out: OT 26.06 maps
# its own "reactome" evidence datasource to affected_pathway (see
# stage1_ingest.py's _DATASOURCE_TO_DATATYPE), so including it here would let
# Reactome-derived evidence weight a graph whose topology already *is*
# Reactome pathway structure — the same orthogonality Stage 2's GSEA
# signature (stage2_gsea_discovery.py) has to preserve, for the same reason.
DEFAULT_DATATYPES = ["genetic_association", "known_drug"]

EDGE_WEIGHT_MODES = ("avg", "product")


# ==========================================================================
# CLI
# ==========================================================================


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="stage4_genetic_evidence_weights.py",
        description="Stage 4: genetic-evidence node/edge weighting over Stage 3's pathway-based "
        "gene-gene graph.",
    )
    p.add_argument("--manifest", required=True, type=Path, help="Path to Stage 1's manifest.json.")
    p.add_argument(
        "--graph-dir",
        type=Path,
        default=None,
        help="Directory containing Stage 3's graph_weight.npz/graph_sign.npz/graph_gene_index.json/"
        "graph_metadata.json. Defaults to alongside --manifest, where stage3_build_graph.py writes them.",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Where to write the updated graph artifacts and gene_weights.tsv. Defaults to --graph-dir.",
    )
    p.add_argument(
        "--datatypes",
        default=",".join(DEFAULT_DATATYPES),
        help="Comma-separated OT datatype_id values admitted as genetic evidence. Default: "
        f"{','.join(DEFAULT_DATATYPES)} — deliberately excludes affected_pathway and any "
        "literature/pathway-derived datatype; see module docstring.",
    )
    p.add_argument(
        "--edge-weight-mode",
        choices=EDGE_WEIGHT_MODES,
        default="avg",
        help="How to combine the two endpoint node scores into an edge weight.",
    )
    return p


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def _fail(msg: str) -> None:
    logger.error(msg)
    sys.exit(1)


# ==========================================================================
# Manifest-driven artifact loading (Open Targets association data)
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


def load_disease_subset(manifest: Manifest, manifest_path: Path) -> pd.DataFrame:
    path = _resolve_artifact_path(manifest, manifest_path, "ot_disease_subset.parquet")
    return OTAssociationsSchema.validate(pd.read_parquet(path))


# ==========================================================================
# Stage 3 graph loading
# ==========================================================================


def load_stage3_graph(graph_dir: Path) -> tuple[sp.csr_matrix, sp.csr_matrix, list[str], dict]:
    weight_path = graph_dir / "graph_weight.npz"
    sign_path = graph_dir / "graph_sign.npz"
    index_path = graph_dir / "graph_gene_index.json"
    metadata_path = graph_dir / "graph_metadata.json"
    for path in (weight_path, sign_path, index_path, metadata_path):
        if not path.exists():
            _fail(f"Stage 3 artifact {path} does not exist (run stage3_build_graph.py first?).")

    weight = sp.load_npz(weight_path).tocsr()
    sign = sp.load_npz(sign_path).tocsr()
    gene_index = json.loads(index_path.read_text())
    metadata = json.loads(metadata_path.read_text())

    n = len(gene_index)
    if weight.shape != (n, n) or sign.shape != (n, n):
        _fail(
            f"Stage 3 artifacts in {graph_dir} are misaligned: weight.shape={weight.shape}, "
            f"sign.shape={sign.shape}, but gene index has {n} entries."
        )
    return weight, sign, gene_index, metadata


# ==========================================================================
# Genetic evidence scores
# ==========================================================================


def compute_genetic_evidence_scores(
    ot_disease_subset: pd.DataFrame, datatypes: list[str], gene_index: list[str]
) -> pd.Series:
    """One score per gene in ``gene_index`` — 0.0 for a gene with no
    admitted-datatype evidence row. The max score wins where a gene has more
    than one admitted row (multiple datasources and/or datatypes); see
    module docstring for why this mirrors stage2_gsea_discovery.py's convention."""
    admitted = ot_disease_subset[ot_disease_subset["datatype_id"].isin(datatypes)]
    if admitted.empty:
        logger.warning(
            "No ot_disease_subset.parquet rows matched --datatypes %s — every gene's "
            "genetic_evidence_score will be 0.0.",
            datatypes,
        )
    per_gene_max = admitted.groupby("gene_id")["score"].max()

    scores = pd.Series(0.0, index=gene_index, name="genetic_evidence_score")
    scores.update(per_gene_max)  # only touches genes already in gene_index; extras are ignored

    return scores


# ==========================================================================
# Edge weights (derived from endpoint node scores only)
# ==========================================================================


def compute_edge_weights(existing_weight: sp.csr_matrix, node_scores: np.ndarray, mode: str) -> sp.csr_matrix:
    """Derive an edge-weight matrix from the two endpoint node scores, at
    exactly the (u, v) positions Stage 3's weight matrix already has a real
    entry for — genetic evidence never adds a new edge, only annotates
    existing ones."""
    coo = existing_weight.tocoo()
    u_scores = node_scores[coo.row]
    v_scores = node_scores[coo.col]
    if mode == "avg":
        data = (u_scores + v_scores) / 2.0
    elif mode == "product":
        data = u_scores * v_scores
    else:
        raise ValueError(f"unknown edge-weight-mode {mode!r}")
    return sp.csr_matrix((data, (coo.row, coo.col)), shape=existing_weight.shape)


# ==========================================================================
# Orchestration
# ==========================================================================


def run(args: argparse.Namespace) -> Path:
    manifest = load_manifest(args.manifest)
    ot_disease_subset = load_disease_subset(manifest, args.manifest)

    graph_dir = args.graph_dir or args.manifest.parent
    weight, sign, gene_index, metadata = load_stage3_graph(graph_dir)

    datatypes = [d.strip() for d in args.datatypes.split(",") if d.strip()]
    if not datatypes:
        _fail("--datatypes resolved to an empty list.")

    scores = compute_genetic_evidence_scores(ot_disease_subset, datatypes, gene_index)
    node_scores = scores.to_numpy()

    edge_weight = compute_edge_weights(weight, node_scores, args.edge_weight_mode)
    logger.info(
        "Genetic evidence: %d/%d genes scored, %d edges weighted (%s).",
        int((node_scores > 0).sum()), len(gene_index), edge_weight.nnz, args.edge_weight_mode,
    )

    out_dir = args.out_dir or graph_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # graph_weight.npz is REPLACED with the genetic-evidence-derived edge
    # weight (same sparsity pattern as Stage 3's structural weight) — see
    # the attribute-name comment at the top of this file and the module
    # docstring's "Output" section for why this reuses Stage 3's "weight"
    # slot instead of a separate file.
    sp.save_npz(out_dir / "graph_weight.npz", edge_weight)
    sp.save_npz(out_dir / "graph_sign.npz", sign)
    (out_dir / "graph_gene_index.json").write_text(json.dumps(gene_index, indent=2))

    metadata = {
        **metadata,
        "genetic_evidence_score": {gene_id: float(score) for gene_id, score in zip(gene_index, node_scores)},
        "stage4_genetic_evidence": {
            "datatypes": datatypes,
            "edge_weight_mode": args.edge_weight_mode,
            "n_genes_with_evidence": int((node_scores > 0).sum()),
            "n_genes_total": len(gene_index),
            "n_edges_weighted": int(edge_weight.nnz),
        },
    }
    (out_dir / "graph_metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True))

    weights_df = pd.DataFrame({"gene_id": gene_index, "genetic_evidence_score": node_scores})
    weights_df = weights_df.sort_values("genetic_evidence_score", ascending=False).reset_index(drop=True)
    weights_df = GeneWeightsSchema.validate(weights_df)
    tsv_path = out_dir / "gene_weights.tsv"
    # Rounded only here, for eyeballing — graph_node_genetic_evidence_score.npy
    # keeps full precision for anything downstream that does arithmetic with it.
    to_write = weights_df.copy()
    to_write["genetic_evidence_score"] = to_write["genetic_evidence_score"].map(lambda v: float(f"{v:.4g}"))
    to_write.to_csv(tsv_path, sep="\t", index=False)

    logger.info("Stage 4 complete: %s", out_dir)
    return out_dir


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    args = parse_args(argv)
    run(args)


if __name__ == "__main__":
    main()
