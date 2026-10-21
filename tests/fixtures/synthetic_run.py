"""Test-only builder for a tiny, hand-designed Stage 1 output directory, used by the tests to exercise
Stages 2-7 end to end without the real Open Targets / Reactome data (see tests/conftest.py's
synthetic_stage1 fixture).

write_synthetic_stage1(out_dir) writes {gene_sets,ot_disease_subset,interactions,genes}.parquet,
ot_target_subset.parquet and manifest.json, validated against stage0_schemas.py exactly like
stage1_ingest.py's real output, so downstream stages exercise the identical code path they use on a real
Stage 1 run.

Gene universe (all fake ENSG ids, 30 genes total). "Seeds" are the genes with
genetic_association evidence, which Stage 2 writes to disease_seed_genes.tsv:
- G1..G11 (ENSG00000000001-11): seeds. G2..G11 are the 10 top
  disease-associated genes by genetic_association score; G1 is mid-range (0.60).
- G12, G13: two weaker seeds, added only to R-HSA-100 to make it a
  near-duplicate *superset* of R-HSA-101.
- G17..G21: weak seeds forming an unrelated pathway.
- N22..N24 (ENSG00000000022-24): non-seed genes on the disease pathways GSEA
  finds - the candidates Stage 5 should rank. N22 has a known_drug row only,
  which Stage 4 ignores (it counts genetic_association evidence only).
- N25..N28: non-seed genes on a pathway with one seed, too few for GSEA - they
  stay out of the Stage 3 graph.
- N29..N33: non-seed genes on a pathway with no seed gene at all - they must
  stay out of the Stage 3 graph too.

Gene sets (mimicking a Reactome ancestor/descendant hierarchy pair):
- R-HSA-100 "Ancestor Broad Signaling" = {G2..G13, N22, N23} (14 genes).
- R-HSA-101 "Child Broad Signaling"    = {G2..G11, N22, N23} (12 genes);
  Jaccard(100, 101) = 12/14 = 0.857 > 0.7, so Stage 2's collapsing drops the
  ancestor (100) before GSEA, so 100 is not in the Stage 3 graph (but its
  genes shared with 101 are).
- R-HSA-200 "Signaling Pathway A" = {G1..G6, N24} (7 genes) - distinct from
  101 (Jaccard = 5/14 = 0.36).
- R-HSA-300 "Unrelated Weak Pathway" = {G17..G21} (5 genes), all weak seeds -
  a *coherent* low-ranked cluster, which GSEA correctly reports as significant
  too (depletion is still enrichment).
- R-HSA-400 "Scattered Noise Pathway" = {G1, G4, G9, G13, G20} spanning the
  top, middle and bottom of the ranking - no coherent enrichment, so it fails
  the significance filter, so it is not in the graph.
- R-HSA-500 "Seed-Adjacent Pathway" = {G1, N25..N28} - only one seed, too few
  for GSEA to test, so it is not in the graph.
- R-HSA-600 "Seedless Pathway" = {N29..N33} - no seed, not significant, so
  it is not in the graph.

One gene (G11) additionally carries a high-scoring *affected_pathway*-datatype
row (not genetic_association) - a deliberate decoy: if stage2_gsea_discovery.py
mistakenly used the overall score instead of genetic_association alone, G11
would rank far higher than its genetic_association score of 0.71 implies.

Also writes genes.parquet (approved symbol "GENE##" per gene) to exercise
--benchmark-holdout-file, which resolves holdout entries via this table. G17's
synonym "OLDG17" is deliberately not its approved symbol, to exercise the
synonym-resolution path.

Also writes interactions.parquet, for Stage 3 (the graph's only edges, all
of weight 1):
- G1 -> G2 (sign=+1, confidence=0.9): both seeds, in the graph.
- G4 -> G9 (sign=-1, confidence=0.6): also in the graph.
- G3 -> N22, G5 -> N22, G7 - N22 (undirected): three seeds reach N22.
- G6 -> N24 (sign=+1): one seed reaches N24.
- N24 -> N23 (sign=+1): N23 is reached only through N24, two steps from a
  seed, so it should rank below N22 and N24.
- G18 - G19 (undirected): within R-HSA-300, the low-ranked cluster.
- G12 -> G1 (sign=+1, confidence=0.99): G12 is only in R-HSA-100, which Stage
  2 collapses away, so this row is dropped as out of scope.
- G2 -> N25 (sign=+1, confidence=0.8): N25 is only in R-HSA-500, which GSEA
  cannot test, so this row is dropped too.
- N29 -> N22 (sign=+1, confidence=0.95): N29 is only in the seedless pathway,
  so this row is dropped as well.

Also writes ot_target_subset.parquet (id / tractability / safetyLiabilities
columns, the Open Targets `target` schema subset Stage 6 reads) and records
it under manifest.sources[opentargets]:
- G1: clinically tractable ("Approved Drug"), no safety liabilities.
- G2: clinically tractable, TWO safety liabilities.
- N22: clinically tractable, one liability -> clinical / has_liabilities.
- N24: discovery-tractable only -> discovery / unknown.
- N25: clinically tractable, no liabilities -> clinical / unknown.
- every other gene is ABSENT from the table -> Stage 6 marks it
  tractability=unknown, safety=unknown and does NOT penalise it.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from stage0_schemas import (
    CoverageEntry,
    GeneSetsSchema,
    GeneSetSizeBucket,
    GenesSchema,
    InteractionsSchema,
    Manifest,
    OTAssociationsSchema,
    OutputArtifact,
    ResolvedDisease,
    ScaleReport,
    Source,
    SourceFile,
)

DISEASE_ID = "MONDO_0000001"


def _gene(n: int) -> str:
    return f"ENSG{n:011d}"


G = {i: _gene(i) for i in [*range(1, 14), *range(17, 22)]}  # seeds G1..G13, G17..G21
N = {i: _gene(i) for i in range(22, 34)}  # non-seed genes N22..N33


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def build_gene_sets() -> pd.DataFrame:
    rows = []

    def add_set(set_id: str, set_name: str, genes: list[str], hierarchy_level: int, parent_id: str | None):
        for gene_id in genes:
            rows.append(
                {
                    "set_id": set_id,
                    "set_name": set_name,
                    "source_db": "reactome",
                    "source_version": "97",
                    "gene_id": gene_id,
                    "hierarchy_level": hierarchy_level,
                    "parent_id": parent_id,
                }
            )

    ancestor_genes = [*[G[i] for i in range(2, 14)], N[22], N[23]]  # G2..G13, N22, N23
    child_genes = [*[G[i] for i in range(2, 12)], N[22], N[23]]  # G2..G11, N22, N23
    pathway_a_genes = [*[G[i] for i in range(1, 7)], N[24]]  # G1..G6, N24
    unrelated_genes = [G[i] for i in range(17, 22)]  # G17..G21
    scattered_genes = [G[1], G[4], G[9], G[13], G[20]]
    seed_adjacent_genes = [G[1], *[N[i] for i in range(25, 29)]]  # G1, N25..N28
    seedless_genes = [N[i] for i in range(29, 34)]  # N29..N33

    add_set("R-HSA-100", "Ancestor Broad Signaling", ancestor_genes, 1, None)
    add_set("R-HSA-101", "Child Broad Signaling", child_genes, 2, "R-HSA-100")
    add_set("R-HSA-200", "Signaling Pathway A", pathway_a_genes, 1, None)
    add_set("R-HSA-300", "Unrelated Weak Pathway", unrelated_genes, 1, None)
    add_set("R-HSA-400", "Scattered Noise Pathway", scattered_genes, 1, None)
    add_set("R-HSA-500", "Seed-Adjacent Pathway", seed_adjacent_genes, 1, None)
    add_set("R-HSA-600", "Seedless Pathway", seedless_genes, 1, None)

    return pd.DataFrame(rows)


def build_genes() -> pd.DataFrame:
    rows = []
    for i, gene_id in {**G, **N}.items():
        synonyms = ["OLDG17"] if i == 17 else []
        rows.append(
            {"gene_id": gene_id, "symbol": f"GENE{i:02d}", "synonyms": synonyms, "biotype": "protein_coding"}
        )
    return pd.DataFrame(rows)


def build_ot_disease_subset() -> pd.DataFrame:
    # (gene, genetic_association score)
    genetic_scores = {
        G[1]: 0.60,
        G[2]: 0.98, G[3]: 0.95, G[4]: 0.92, G[5]: 0.89, G[6]: 0.86,
        G[7]: 0.83, G[8]: 0.80, G[9]: 0.77, G[10]: 0.74, G[11]: 0.71,
        G[12]: 0.50, G[13]: 0.45,
        G[17]: 0.05, G[18]: 0.04, G[19]: 0.03, G[20]: 0.02, G[21]: 0.06,
    }
    rows = [
        {
            "gene_id": gene_id,
            "disease_id": DISEASE_ID,
            "datatype_id": "genetic_association",
            "datasource_id": "eva",
            "score": score,
        }
        for gene_id, score in genetic_scores.items()
    ]
    # N22 is not a seed (no genetic_association row); its known_drug evidence must not count in Stage 4.
    rows.append(
        {"gene_id": N[22], "disease_id": DISEASE_ID, "datatype_id": "known_drug", "datasource_id": "chembl",
         "score": 0.50}
    )
    # Decoy: G11's aggregated/affected_pathway-datatype evidence is much
    # stronger than its genetic_association score. stage2_gsea_discovery.py must
    # ignore this row when building the ranking signature.
    rows.append(
        {
            "gene_id": G[11],
            "disease_id": DISEASE_ID,
            "datatype_id": "affected_pathway",
            "datasource_id": "reactome",
            "score": 0.99,
        }
    )
    # Literature evidence, as Open Targets has for many well-studied genes: stronger for some strong seeds,
    # weak for a few weak ones, and present for two candidates. The pipeline never reads it (Stages 2 and 4
    # use genetic_association only); validation/compare_scores.py ranks genes by it as the Open Targets baseline.
    literature_scores = {
        G[2]: 0.60, G[3]: 0.20, G[5]: 0.45, G[8]: 0.10, G[9]: 0.30,
        G[17]: 0.15, G[20]: 0.05,
        N[23]: 0.25, N[27]: 0.40,
    }
    for gene_id, score in literature_scores.items():
        rows.append({"gene_id": gene_id, "disease_id": DISEASE_ID, "datatype_id": "literature",
                     "datasource_id": "europepmc", "score": score})
    return pd.DataFrame(rows)


def build_interactions() -> pd.DataFrame:
    rows = [
        {"gene_a": G[1], "gene_b": G[2], "directed": True, "sign": 1, "source_db": "reactome",
         "source_version": "04142025", "evidence_type": "curated", "confidence": 0.9},
        {"gene_a": G[4], "gene_b": G[9], "directed": True, "sign": -1, "source_db": "reactome",
         "source_version": "04142025", "evidence_type": "curated", "confidence": 0.6},
        # Seeds reach N22 and N24 directly; N23 only through N24.
        {"gene_a": G[3], "gene_b": N[22], "directed": True, "sign": 1, "source_db": "reactome",
         "source_version": "04142025", "evidence_type": "curated", "confidence": 0.9},
        {"gene_a": G[5], "gene_b": N[22], "directed": True, "sign": -1, "source_db": "reactome",
         "source_version": "04142025", "evidence_type": "curated", "confidence": 0.7},
        {"gene_a": G[7], "gene_b": N[22], "directed": False, "sign": 0, "source_db": "reactome",
         "source_version": "04142025", "evidence_type": "predicted", "confidence": 0.5},
        {"gene_a": G[6], "gene_b": N[24], "directed": True, "sign": 1, "source_db": "reactome",
         "source_version": "04142025", "evidence_type": "curated", "confidence": 0.8},
        {"gene_a": N[24], "gene_b": N[23], "directed": True, "sign": 1, "source_db": "reactome",
         "source_version": "04142025", "evidence_type": "curated", "confidence": 0.8},
        {"gene_a": G[18], "gene_b": G[19], "directed": False, "sign": 0, "source_db": "reactome",
         "source_version": "04142025", "evidence_type": "predicted", "confidence": 0.6},
        # Out of scope for Stage 3: G12 is only in R-HSA-100, which Stage 2 collapses away.
        {"gene_a": G[12], "gene_b": G[1], "directed": True, "sign": 1, "source_db": "reactome",
         "source_version": "04142025", "evidence_type": "curated", "confidence": 0.99},
        # Out of scope for Stage 3: N25 is only in R-HSA-500, which GSEA cannot test.
        {"gene_a": G[2], "gene_b": N[25], "directed": True, "sign": 1, "source_db": "reactome",
         "source_version": "04142025", "evidence_type": "curated", "confidence": 0.8},
        # Out of scope for Stage 3: N29 is only in the seedless pathway, so this row must be dropped.
        {"gene_a": N[29], "gene_b": N[22], "directed": True, "sign": 1, "source_db": "reactome",
         "source_version": "04142025", "evidence_type": "curated", "confidence": 0.95},
    ]
    return pd.DataFrame(rows)


def build_target_annotations() -> pd.DataFrame:
    """Open Targets `target` schema subset — only the columns Stage 6 reads."""

    def tract(*flags: tuple[str, str, bool]) -> list[dict]:
        return [{"modality": m, "id": i, "value": v} for (m, i, v) in flags]

    def liabs(*events: str) -> list[dict]:
        return [{"event": e, "eventId": f"EFO_{n:07d}"} for n, e in enumerate(events, start=1)]

    rows = [
        {"id": G[1],
         "tractability": tract(("SM", "Approved Drug", True), ("SM", "Structure with Ligand", True)),
         "safetyLiabilities": liabs()},
        {"id": G[2],
         "tractability": tract(("SM", "Phase 1 Clinical", True), ("AB", "GO CC high conf", True)),
         "safetyLiabilities": liabs("cardiovascular toxicity", "hepatotoxicity")},
        {"id": N[22],
         "tractability": tract(("SM", "Phase 2 Clinical", True)),
         "safetyLiabilities": liabs("nephrotoxicity")},
        {"id": N[24],
         "tractability": tract(("SM", "Structure with Ligand", True), ("SM", "High-Quality Pocket", False)),
         "safetyLiabilities": liabs()},
        {"id": N[25],
         "tractability": tract(("SM", "Approved Drug", True)),
         "safetyLiabilities": liabs()},
    ]
    return pd.DataFrame(rows, columns=["id", "tractability", "safetyLiabilities"])


def write_synthetic_stage1(out_dir: Path) -> Path:
    """Write the synthetic Stage 1 output into out_dir and return the manifest path."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    gene_sets_df = GeneSetsSchema.validate(build_gene_sets())
    gene_sets_path = out_dir / "gene_sets.parquet"
    gene_sets_df.to_parquet(gene_sets_path, index=False)

    ot_subset_df = OTAssociationsSchema.validate(build_ot_disease_subset())
    ot_subset_path = out_dir / "ot_disease_subset.parquet"
    ot_subset_df.to_parquet(ot_subset_path, index=False)

    genes_df = GenesSchema.validate(build_genes())
    genes_path = out_dir / "genes.parquet"
    genes_df.to_parquet(genes_path, index=False)

    interactions_df = InteractionsSchema.validate(build_interactions())
    interactions_path = out_dir / "interactions.parquet"
    interactions_df.to_parquet(interactions_path, index=False)

    # Not a Stage 1 output artifact — this stands in for the Open Targets
    # `target` parquet Stage 1 downloads, so it goes under manifest.sources.
    target_subset_path = out_dir / "ot_target_subset.parquet"
    build_target_annotations().to_parquet(target_subset_path, index=False)

    def artifact(path: Path) -> OutputArtifact:
        # A bare file name: stages look for each artifact next to the manifest first.
        return OutputArtifact(path=path.name, sha256=_sha256_bytes(path.read_bytes()))

    n_associated = int(ot_subset_df["gene_id"].nunique())

    manifest = Manifest(
        run_id="synthetic000000001",
        git_commit="unknown",
        created_at=datetime.now(timezone.utc),
        seed=0,
        disease=ResolvedDisease(efo_id=DISEASE_ID, name="synthetic disease", n_associated_genes=n_associated),
        sources=[
            Source(
                db="reactome",
                version="97",
                files=[SourceFile(path="fake_reactome.gmt.zip", sha256="0" * 64, bytes=0)],
            ),
            Source(
                db="opentargets",
                version="26.06",
                files=[
                    SourceFile(path="fake_association.parquet", sha256="1" * 64, bytes=0),
                    SourceFile(
                        # Stage 6 looks for this basename next to the manifest.
                        path="ot_target_subset.parquet",
                        sha256=_sha256_bytes(target_subset_path.read_bytes()),
                        bytes=target_subset_path.stat().st_size,
                    ),
                ],
            ),
        ],
        cli_parameters={"disease": DISEASE_ID, "min_set_size": 1, "max_set_size": 200},
        output_artifacts=[
            artifact(gene_sets_path), artifact(ot_subset_path), artifact(genes_path), artifact(interactions_path),
        ],
        coverage_report={
            "reactome_gmt": CoverageEntry(
                source="reactome_gmt", raw_gene_count=30, mapped_count=30, percent_mapped=100.0, dropped_count=0
            ),
        },
        scale_report=ScaleReport(
            gene_set_size_distribution_before_cap=[GeneSetSizeBucket(set_size=6, count=7)],
            gene_set_size_distribution_after_cap=[GeneSetSizeBucket(set_size=6, count=7)],
            sets_retained=7,
            interaction_counts_by_sign={},
            interaction_counts_by_evidence_type={},
        ),
    )
    manifest_path = out_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest.model_dump(mode="json"), indent=2, sort_keys=True))
    return manifest_path
