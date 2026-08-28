"""One-off generator for a tiny, hand-designed Stage 1 output directory used
to exercise stage2_gsea_discovery.py, stage3_build_graph.py, and stage4_genetic_evidence_weights.py
end to end (see README.md's Stage 2/3/4 run instructions). Not run by any
test suite — re-run manually if the example needs to change.

Writes example_data/stage1_run/{gene_sets,ot_disease_subset,interactions}.parquet
and manifest.json, validated against stage0_schemas.py exactly like stage1_ingest.py's real
output, so downstream stages exercise the identical code path they use on a
real Stage 1 run.

Gene universe (all fake ENSG ids, 17 genes total):
- TARGET (ENSG00000000001): the --target gene. Has a mid-range
  genetic_association score and sits in one pathway (R-HSA-200).
- G2..G11 (ENSG00000000002-11): the 10 top disease-associated genes by
  genetic_association score.
- G12, G13: two more, weaker disease-associated genes, added only to
  R-HSA-100 to make it a near-duplicate *superset* of R-HSA-101.
- G17..G21: weakly disease-associated genes forming an unrelated pathway.

Gene sets (mimicking a Reactome ancestor/descendant hierarchy pair):
- R-HSA-100 "Ancestor Broad Signaling" = {G2..G13} (12 genes) — the ancestor.
- R-HSA-101 "Child Broad Signaling"    = {G2..G11} (10 genes) — its child;
  Jaccard(100, 101) = 10/12 = 0.833 > 0.7, so collapsing must drop the
  ancestor (100) and keep the smaller, more specific child (101).
- R-HSA-200 "Target Signaling Pathway" = {TARGET, G2..G6} (6 genes) —
  distinct from 101 (Jaccard = 5/11 = 0.45), contains the target gene.
- R-HSA-300 "Unrelated Weak Pathway"   = {G17..G21} (5 genes), all weakly
  disease-associated — a *coherent* low-ranked cluster, which GSEA correctly
  reports as significant too (depletion is still enrichment); it is not a
  negative control by itself.
- R-HSA-400 "Scattered Noise Pathway"  = {TARGET, G4, G9, G13, G20} (5 genes,
  the minimum blitzgsea will test) spanning the top, middle, and bottom of
  the ranking rather than clustering anywhere — the actual negative control:
  no coherent enrichment signal, so it should fail the significance filter
  and demonstrate --pval-threshold/--fdr-threshold actually excluding a
  tested pathway.

One gene (G11) additionally carries a high-scoring *affected_pathway*-datatype
row (not genetic_association) — a deliberate decoy: if stage2_gsea_discovery.py
mistakenly used the aggregated/overall score instead of genetic_association
alone, G11 would rank far higher than its genetic_association score of 0.71
implies.

Also writes genes.parquet (approved symbol "GENE##" per gene, TARGET's
symbol "TARGET1") to exercise --benchmark-holdout-file, which resolves
holdout entries via this table. G17's synonym "OLDG17" is deliberately not
its approved symbol, to exercise the synonym-resolution path a real
literature-curated benchmark file might need.

Also writes interactions.parquet, for stage3_build_graph.py (Stage 3) and
stage4_genetic_evidence_weights.py (Stage 4):
- TARGET -> G2 (sign=+1, confidence=0.9): both genes are in the pathway
  union (R-HSA-200), so this overrides that pair's co-membership sign.
- G4 -> G9 (sign=-1, confidence=0.6): also both in the union (R-HSA-101).
- G12 -> TARGET (sign=+1, confidence=0.99): G12 is a member of R-HSA-100
  only, which is excluded from Stage 3's pathway union entirely (it's the
  redundant ancestor Stage 2's collapsing drops) — this row must be
  dropped as out of scope, not silently included.
"""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / "pipeline"))
from stage0_schemas import (  # noqa: E402
    CoverageEntry,
    GeneSetsSchema,
    GeneSetSizeBucket,
    GenesSchema,
    InteractionsSchema,
    Manifest,
    OTAssociationsSchema,
    OutputArtifact,
    ResolvedDisease,
    ResolvedTarget,
    ScaleReport,
    Source,
    SourceFile,
)

OUT_DIR = Path(__file__).parent / "stage1_run"

TARGET = "ENSG00000000001"
DISEASE_ID = "MONDO_0000001"


def _gene(n: int) -> str:
    return f"ENSG{n:011d}"


G = {i: _gene(i) for i in range(2, 22)}  # G2..G21


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

    ancestor_genes = [G[i] for i in range(2, 14)]  # G2..G13
    child_genes = [G[i] for i in range(2, 12)]  # G2..G11
    target_genes = [TARGET, *[G[i] for i in range(2, 7)]]  # TARGET, G2..G6
    unrelated_genes = [G[i] for i in range(17, 22)]  # G17..G21
    scattered_genes = [TARGET, G[4], G[9], G[13], G[20]]

    add_set("R-HSA-100", "Ancestor Broad Signaling", ancestor_genes, 1, None)
    add_set("R-HSA-101", "Child Broad Signaling", child_genes, 2, "R-HSA-100")
    add_set("R-HSA-200", "Target Signaling Pathway", target_genes, 1, None)
    add_set("R-HSA-300", "Unrelated Weak Pathway", unrelated_genes, 1, None)
    add_set("R-HSA-400", "Scattered Noise Pathway", scattered_genes, 1, None)

    return pd.DataFrame(rows)


def build_genes() -> pd.DataFrame:
    rows = [{"gene_id": TARGET, "symbol": "TARGET1", "synonyms": [], "biotype": "protein_coding"}]
    for i, gene_id in G.items():
        synonyms = ["OLDG17"] if i == 17 else []
        rows.append(
            {"gene_id": gene_id, "symbol": f"GENE{i:02d}", "synonyms": synonyms, "biotype": "protein_coding"}
        )
    return pd.DataFrame(rows)


def build_ot_disease_subset() -> pd.DataFrame:
    # (gene, genetic_association score)
    genetic_scores = {
        TARGET: 0.60,
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
    return pd.DataFrame(rows)


def build_interactions() -> pd.DataFrame:
    rows = [
        {"gene_a": TARGET, "gene_b": G[2], "directed": True, "sign": 1, "source_db": "reactome",
         "source_version": "04142025", "evidence_type": "curated", "confidence": 0.9},
        {"gene_a": G[4], "gene_b": G[9], "directed": True, "sign": -1, "source_db": "reactome",
         "source_version": "04142025", "evidence_type": "curated", "confidence": 0.6},
        # Out of scope for Stage 3: G12 belongs only to R-HSA-100, the
        # ancestor pathway excluded from the union — this row must be
        # dropped, not silently included.
        {"gene_a": G[12], "gene_b": TARGET, "directed": True, "sign": 1, "source_db": "reactome",
         "source_version": "04142025", "evidence_type": "curated", "confidence": 0.99},
    ]
    return pd.DataFrame(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    gene_sets_df = GeneSetsSchema.validate(build_gene_sets())
    gene_sets_path = OUT_DIR / "gene_sets.parquet"
    gene_sets_df.to_parquet(gene_sets_path, index=False)

    ot_subset_df = OTAssociationsSchema.validate(build_ot_disease_subset())
    ot_subset_path = OUT_DIR / "ot_disease_subset.parquet"
    ot_subset_df.to_parquet(ot_subset_path, index=False)

    genes_df = GenesSchema.validate(build_genes())
    genes_path = OUT_DIR / "genes.parquet"
    genes_df.to_parquet(genes_path, index=False)

    interactions_df = InteractionsSchema.validate(build_interactions())
    interactions_path = OUT_DIR / "interactions.parquet"
    interactions_df.to_parquet(interactions_path, index=False)

    def artifact(path: Path) -> OutputArtifact:
        return OutputArtifact(path=str(path), sha256=_sha256_bytes(path.read_bytes()))

    n_associated = int(ot_subset_df["gene_id"].nunique())

    manifest = Manifest(
        run_id="example0000000001",
        git_commit="unknown",
        created_at=datetime.now(timezone.utc),
        seed=0,
        resolved_target=ResolvedTarget(input="TARGET", gene_id=TARGET, symbol="TARGET"),
        disease=ResolvedDisease(efo_id=DISEASE_ID, name="example disease", n_associated_genes=n_associated),
        sources=[
            Source(
                db="reactome",
                version="97",
                files=[SourceFile(path="example_data/fake_reactome.gmt.zip", sha256="0" * 64, bytes=0)],
            ),
            Source(
                db="opentargets",
                version="26.06",
                files=[SourceFile(path="example_data/fake_association.parquet", sha256="1" * 64, bytes=0)],
            ),
        ],
        cli_parameters={"target": "TARGET", "disease": DISEASE_ID, "min_set_size": 1, "max_set_size": 200},
        output_artifacts=[
            artifact(gene_sets_path), artifact(ot_subset_path), artifact(genes_path), artifact(interactions_path),
        ],
        coverage_report={
            "reactome_gmt": CoverageEntry(
                source="reactome_gmt", raw_gene_count=21, mapped_count=21, percent_mapped=100.0, dropped_count=0
            ),
        },
        scale_report=ScaleReport(
            gene_set_size_distribution_before_cap=[GeneSetSizeBucket(set_size=6, count=4)],
            gene_set_size_distribution_after_cap=[GeneSetSizeBucket(set_size=6, count=4)],
            sets_retained=4,
            interaction_counts_by_sign={},
            interaction_counts_by_evidence_type={},
            projected_comembership_edge_count=0,
        ),
    )
    manifest_path = OUT_DIR / "manifest.json"
    manifest_path.write_text(json.dumps(manifest.model_dump(mode="json"), indent=2, sort_keys=True))
    print(f"Wrote {manifest_path}")


if __name__ == "__main__":
    main()
