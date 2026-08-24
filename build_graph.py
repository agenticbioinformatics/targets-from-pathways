"""Stage 3 — Pathway-based gene-gene graph construction.

Builds the **pathway-based gene-gene graph**: a directed, signed graph whose
node set is exactly the genes appearing in (Stage 2's disease-relevant
pathways) ∪ (the target's own pathways) — never the whole genome, never the
whole interactome. That scoping is what makes it "pathway-based" rather than
a generic PPI/functional-interaction network with a pathway-derived weight
bolted on; see "Interactions are scoped to the pathway-gene pool" below for
why this matters concretely.

Inputs, both read only through the artifacts named — never a raw source
file, and never by re-deriving anything Stage 1/2 already computed:

- ``--manifest``: Stage 1's manifest.json, for ``gene_sets.parquet``,
  ``interactions.parquet`` (already parsed, filtered, and canonicalized —
  this module does not touch the Reactome FI file), the resolved target
  gene, and (for graph metadata) the Reactome version(s) and random seed.
- ``--disease-pathways``: Stage 2's ``disease_pathways.tsv``, validated
  against ``schemas.DiseasePathwaysSchema``. Stage 2 doesn't register this
  file as a manifest output artifact (it isn't part of Stage 1's contract),
  so it defaults to sitting alongside ``manifest.json`` — where
  ``gsea_discovery.py`` writes it — with an explicit override for when
  it doesn't.

Two edges sources, both added to one ``networkx.DiGraph``:

1. **Co-membership**: for every pathway in the union above, every pair of
   its genes gets a bidirectional pair of edges (sign=0, "unknown effect")
   weighted ``1 / (|pathway| - 1)``, *summed* over every qualifying pathway
   that co-contains the pair. Without the down-weighting, one 200-gene
   pathway alone contributes ~20,000 edges and a handful of large Reactome
   sets would dominate the entire graph's topology — the hub-domination
   failure mode behind v1's results. The realized edge count is logged
   against Stage 1's ``scale_report.json`` projection; the projection
   necessarily runs far higher, since it sums over *every* retained gene
   set while this stage only ever uses the small disease+target union.
2. **Interactions**: directional, signed rows read straight from
   ``interactions.parquet``. Which curation those rows come from is
   README.md's open decision ("what counts as pathway topology?") and this
   module does not resolve it — it reads whatever ``source_db``-tagged,
   already-normalized rows Stage 1 wrote, unchanged by which option gets
   picked later.

**Interactions are scoped to the pathway-gene pool.** A row is only added
if *both* ``gene_a`` and ``gene_b`` already belong to the disease+target
pathway union — never genome-wide. ``interactions.parquet`` is a
disease-agnostic, curated Reactome-FI-derived table; adding every row
unconditionally would silently turn a "pathway-based" graph back into the
whole interactome (hundreds of thousands of edges unrelated to this
disease or target), exactly the scale/hub-domination problem the
co-membership weighting exists to avoid, and would contradict README.md's
own description of this graph as one "spanning the target's and the
disease's pathway neighborhoods." Out-of-scope rows are counted and
logged, not silently dropped.

**Duplicate interaction rows.** ``interactions.parquet`` can legitimately
carry more than one row for the same ordered ``(gene_a, gene_b)`` pair —
distinct Reactome FI reaction annotations can produce different
sign/confidence values for the same gene pair, and Stage 1 only
deduplicates exact full-row duplicates. Since a plain ``DiGraph`` (as
specified — not a ``MultiDiGraph``) can hold only one edge per ordered
pair, the highest-``confidence`` row wins; how many pairs this affected is
logged, never silent.

When both edge sources touch the same ordered pair, the interaction row's
``sign`` overrides co-membership's default ``sign=0`` (it is more specific
evidence); the co-membership weight is preserved either way.

**Serialization**: two ``scipy.sparse`` CSR matrices sharing one gene index
(an ordered list, row/col ``i`` <-> ``gene_index[i]``), not GraphML —
GraphML on a graph with hundreds of thousands of edges is slow to write,
slow to re-read, and large on disk, and nothing downstream needs anything
GraphML gives beyond what the matrices plus index already carry.
``networkx`` is used only for construction, in-process inspection, and
small test fixtures; it is never the serialized artifact.

- ``weight`` (float64): ``comembership_weight + confidence`` (confidence
  treated as 0 where null) — a single combined closeness score for
  anything that walks the graph by weight, Stage 5's diffusion in
  particular. A true stored zero means "no edge."
- ``sign`` (int8): the edge's sign where one exists, ``-1``/``0``/``1``.
  A stored ``0`` is ambiguous by construction — it means either "no edge"
  or "an edge exists with unknown/co-membership-only effect" — so
  ``sign`` is only meaningful at a position where ``weight > 0``; check
  ``weight`` first. Edge *direction* is not stored as a separate matrix at
  all: it's already exactly what row/col asymmetry between ``weight[i,j]``
  and ``weight[j,i]`` encodes in a directed adjacency matrix, so a third
  matrix would only duplicate that information.

Graph metadata (seed, Reactome source versions, target/disease identity,
edge-count comparison against the Stage 1 projection) is written alongside
as ``graph_metadata.json``.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from datetime import datetime, timezone
from pathlib import Path

import networkx as nx
import pandas as pd
import scipy.sparse as sp

from schemas import (
    DiseasePathwaysSchema,
    ENSEMBL_GENE_ID_PATTERN,
    GeneSetsSchema,
    InteractionsSchema,
    Manifest,
)

logger = logging.getLogger("build_graph")


# ==========================================================================
# CLI
# ==========================================================================


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="build_graph.py",
        description="Stage 3: pathway-based gene-gene graph construction over a Stage 1 manifest "
        "and Stage 2's disease-relevant pathway list.",
    )
    p.add_argument("--manifest", required=True, type=Path, help="Path to Stage 1's manifest.json.")
    p.add_argument(
        "--disease-pathways",
        type=Path,
        default=None,
        help="Path to Stage 2's disease_pathways.tsv. Defaults to alongside --manifest, where "
        "gsea_discovery.py writes it.",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Where to write graph artifacts. Defaults to --manifest's directory.",
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


def load_interactions(manifest: Manifest, manifest_path: Path) -> pd.DataFrame:
    path = _resolve_artifact_path(manifest, manifest_path, "interactions.parquet")
    return InteractionsSchema.validate(pd.read_parquet(path))


def load_disease_pathways(path: Path) -> pd.DataFrame:
    if not path.exists():
        _fail(f"--disease-pathways {path} does not exist (run gsea_discovery.py first?).")
    return DiseasePathwaysSchema.validate(pd.read_csv(path, sep="\t"))


# ==========================================================================
# Pathway union: (disease-relevant) ∪ (target-containing)
# ==========================================================================


def compute_target_pathway_ids(gene_sets_df: pd.DataFrame, target_gene_id: str) -> set[str]:
    ids = set(gene_sets_df.loc[gene_sets_df["gene_id"] == target_gene_id, "set_id"])
    if not ids:
        logger.warning(
            "Target gene %s is not a member of any pathway in gene_sets.parquet — "
            "the target-containing side of the union is empty (README.md open decision #3: "
            "no additional cap is applied here beyond Stage 1's own gene_sets.parquet size cap).",
            target_gene_id,
        )
    return ids


def compute_pathway_union(
    gene_sets_df: pd.DataFrame, disease_pathways_df: pd.DataFrame, target_gene_id: str
) -> tuple[set[str], set[str], set[str]]:
    disease_ids = set(disease_pathways_df["set_id"])
    target_ids = compute_target_pathway_ids(gene_sets_df, target_gene_id)
    union_ids = disease_ids | target_ids
    logger.info(
        "Pathway union: %d disease-relevant + %d target-containing = %d unique pathways "
        "(%d overlap both).",
        len(disease_ids), len(target_ids), len(union_ids), len(disease_ids & target_ids),
    )
    return union_ids, disease_ids, target_ids


# ==========================================================================
# Co-membership edges
# ==========================================================================


def add_comembership_edges(graph: nx.DiGraph, gene_sets_df: pd.DataFrame, pathway_union_ids: set[str]) -> int:
    """Add every pathway-union pathway's co-membership edges, weighted
    ``1 / (|pathway| - 1)`` per pathway and summed where multiple qualifying
    pathways co-contain the same pair. Returns the number of unordered gene
    pairs that received at least one contribution (each becomes two directed
    edges, added to ``graph`` symmetrically with sign=0)."""
    subset = gene_sets_df[gene_sets_df["set_id"].isin(pathway_union_ids)]
    graph.add_nodes_from(subset["gene_id"].unique())

    unordered_pairs_seen: set[tuple[str, str]] = set()
    n_pathways_used = 0
    for set_id, group in subset.groupby("set_id"):
        genes = sorted(group["gene_id"].unique())
        n = len(genes)
        if n < 2:
            continue
        n_pathways_used += 1
        contribution = 1.0 / (n - 1)
        source_db = group["source_db"].iloc[0]
        for i in range(n):
            for j in range(i + 1, n):
                a, b = genes[i], genes[j]
                unordered_pairs_seen.add((a, b))
                for u, v in ((a, b), (b, a)):
                    if graph.has_edge(u, v):
                        graph[u][v]["comembership_weight"] += contribution
                    else:
                        graph.add_edge(
                            u, v,
                            comembership_weight=contribution,
                            sign=0,
                            confidence=None,
                            evidence_type=None,
                            source_db=source_db,
                        )
    logger.info(
        "Co-membership: %d pathways contributed edges (of %d in the union; the rest had <2 "
        "genes after collapsing/mapping), inducing %d unordered gene pairs (%d directed edges).",
        n_pathways_used, len(pathway_union_ids), len(unordered_pairs_seen), 2 * len(unordered_pairs_seen),
    )
    return len(unordered_pairs_seen)


# ==========================================================================
# Interaction edges
# ==========================================================================


def dedupe_interactions_by_confidence(interactions_df: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """Collapse to one row per ordered (gene_a, gene_b): highest confidence
    wins, ties broken by original row order. See module docstring
    "Duplicate interaction rows" for why this is needed at all."""
    if interactions_df.empty:
        return interactions_df, 0
    n_rows_before = len(interactions_df)
    ordered = interactions_df.assign(_rank=interactions_df["confidence"].fillna(-1.0))
    ordered = ordered.sort_values("_rank", ascending=False, kind="stable")
    deduped = ordered.drop_duplicates(subset=["gene_a", "gene_b"], keep="first").drop(columns="_rank")
    n_dropped = n_rows_before - len(deduped)
    return deduped.reset_index(drop=True), n_dropped


def add_interaction_edges(
    graph: nx.DiGraph, interactions_df: pd.DataFrame, pathway_gene_pool: set[str]
) -> tuple[int, int, int]:
    """Add interactions.parquet rows as directed, signed edges, restricted to
    pairs where both genes are already in the pathway union (see module
    docstring "Interactions are scoped to the pathway-gene pool"). Returns
    (n_added_or_updated, n_out_of_scope_dropped, n_duplicate_pairs_collapsed)."""
    deduped, n_duplicate_pairs = dedupe_interactions_by_confidence(interactions_df)
    if deduped.empty:
        return 0, 0, n_duplicate_pairs

    in_scope = deduped["gene_a"].isin(pathway_gene_pool) & deduped["gene_b"].isin(pathway_gene_pool)
    n_out_of_scope = int((~in_scope).sum())
    scoped = deduped.loc[in_scope]

    for row in scoped.itertuples(index=False):
        confidence = None if pd.isna(row.confidence) else float(row.confidence)
        if graph.has_edge(row.gene_a, row.gene_b):
            graph[row.gene_a][row.gene_b].update(
                sign=int(row.sign), confidence=confidence, evidence_type=row.evidence_type, source_db=row.source_db
            )
        else:
            graph.add_edge(
                row.gene_a, row.gene_b,
                comembership_weight=0.0,
                sign=int(row.sign),
                confidence=confidence,
                evidence_type=row.evidence_type,
                source_db=row.source_db,
            )
    logger.info(
        "Interactions: %d rows in scope (both genes in the pathway union) added/updated as edges; "
        "%d row(s) out of scope dropped (at least one gene outside the union — see module docstring); "
        "%d duplicate (gene_a, gene_b) pair(s) collapsed to their highest-confidence row.",
        len(scoped), n_out_of_scope, n_duplicate_pairs,
    )
    return len(scoped), n_out_of_scope, n_duplicate_pairs


# ==========================================================================
# Orchestration: build the pathway-based gene-gene graph
# ==========================================================================


def build_pathway_based_gene_gene_graph(
    gene_sets_df: pd.DataFrame, disease_pathways_df: pd.DataFrame, interactions_df: pd.DataFrame, target_gene_id: str
) -> tuple[nx.DiGraph, dict]:
    pathway_union_ids, disease_ids, target_ids = compute_pathway_union(
        gene_sets_df, disease_pathways_df, target_gene_id
    )

    pathway_based_gene_gene_graph = nx.DiGraph()
    n_unordered_pairs = add_comembership_edges(pathway_based_gene_gene_graph, gene_sets_df, pathway_union_ids)
    pathway_gene_pool = set(pathway_based_gene_gene_graph.nodes)
    n_interaction_edges, n_out_of_scope, n_duplicate_pairs = add_interaction_edges(
        pathway_based_gene_gene_graph, interactions_df, pathway_gene_pool
    )

    stats = {
        "n_pathways_union": len(pathway_union_ids),
        "n_disease_pathways": len(disease_ids),
        "n_target_pathways": len(target_ids),
        "n_nodes": pathway_based_gene_gene_graph.number_of_nodes(),
        "n_edges": pathway_based_gene_gene_graph.number_of_edges(),
        "n_comembership_unordered_pairs": n_unordered_pairs,
        "n_interaction_edges_added": n_interaction_edges,
        "n_interaction_edges_out_of_scope": n_out_of_scope,
        "n_interaction_duplicate_pairs_collapsed": n_duplicate_pairs,
    }
    return pathway_based_gene_gene_graph, stats


def log_edge_count_vs_projection(stats: dict, manifest: Manifest) -> None:
    projected = manifest.scale_report.projected_comembership_edge_count
    logger.info(
        "Co-membership edge count: %d unordered pairs realized here vs. %d projected in Stage 1's "
        "scale_report.json. The projection sums over *every* retained gene set (%d sets); this "
        "stage only uses the %d-pathway disease+target union, so a much smaller realized count is "
        "expected, not a red flag by itself.",
        stats["n_comembership_unordered_pairs"], projected,
        manifest.scale_report.sets_retained, stats["n_pathways_union"],
    )


# ==========================================================================
# Sparse serialization
# ==========================================================================


def validate_gene_index(gene_index: list[str]) -> None:
    if len(gene_index) != len(set(gene_index)):
        _fail("Gene index contains duplicate gene IDs — this would silently misalign matrix rows.")
    bad = [g for g in gene_index if not ENSEMBL_GENE_ID_PATTERN.match(g)]
    if bad:
        _fail(f"Gene index contains non-Ensembl gene ID(s): {bad[:20]}")


def graph_to_sparse(graph: nx.DiGraph) -> tuple[sp.csr_matrix, sp.csr_matrix, list[str]]:
    """(weight, sign, gene_index) — see module docstring "Serialization" for
    the exact contract, in particular that ``sign`` is only meaningful where
    the corresponding ``weight`` entry is > 0."""
    gene_index = sorted(graph.nodes)
    validate_gene_index(gene_index)

    for _u, _v, attrs in graph.edges(data=True):
        confidence = attrs["confidence"] if attrs["confidence"] is not None else 0.0
        attrs["_matrix_weight"] = attrs["comembership_weight"] + confidence

    weight = nx.to_scipy_sparse_array(graph, nodelist=gene_index, weight="_matrix_weight", format="csr", dtype="float64")
    sign = nx.to_scipy_sparse_array(graph, nodelist=gene_index, weight="sign", format="csr", dtype="int8")
    return weight, sign, gene_index


# ==========================================================================
# Metadata
# ==========================================================================


def build_graph_metadata(manifest: Manifest, stats: dict) -> dict:
    reactome_sources = [s.model_dump(mode="json") for s in manifest.sources if s.db == "reactome"]
    return {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "seed": manifest.seed,
        "reactome_sources": reactome_sources,
        "resolved_target": manifest.resolved_target.model_dump(mode="json"),
        "disease": manifest.disease.model_dump(mode="json"),
        **stats,
    }


# ==========================================================================
# Orchestration
# ==========================================================================


def run(args: argparse.Namespace) -> Path:
    manifest = load_manifest(args.manifest)
    gene_sets_df = load_gene_sets(manifest, args.manifest)
    interactions_df = load_interactions(manifest, args.manifest)
    disease_pathways_path = args.disease_pathways or (args.manifest.parent / "disease_pathways.tsv")
    disease_pathways_df = load_disease_pathways(disease_pathways_path)
    target_gene_id = manifest.resolved_target.gene_id

    pathway_based_gene_gene_graph, stats = build_pathway_based_gene_gene_graph(
        gene_sets_df, disease_pathways_df, interactions_df, target_gene_id
    )
    log_edge_count_vs_projection(stats, manifest)
    logger.info(
        "Pathway-based gene-gene graph: %d nodes, %d directed edges.",
        stats["n_nodes"], stats["n_edges"],
    )

    weight, sign, gene_index = graph_to_sparse(pathway_based_gene_gene_graph)
    metadata = build_graph_metadata(manifest, stats)

    out_dir = args.out_dir or args.manifest.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    sp.save_npz(out_dir / "graph_weight.npz", weight)
    sp.save_npz(out_dir / "graph_sign.npz", sign)
    (out_dir / "graph_gene_index.json").write_text(json.dumps(gene_index, indent=2))
    metadata_path = out_dir / "graph_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True))
    logger.info("Stage 3 complete: %s", metadata_path)
    return metadata_path


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    args = parse_args(argv)
    run(args)


if __name__ == "__main__":
    main()
