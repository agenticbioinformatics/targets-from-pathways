"""End-to-end schema-conformance test across the whole pipeline.

The highest-value test in the suite: it runs the real stage scripts in
order and re-reads **every inter-stage artifact**, validating each against
its ``stage0_schemas`` contract. A renamed or dropped column in any stage's
output, or a gene id that stopped being an Ensembl id, or ``source_db``
attribution getting lost, fails here — not on Day 3 morning.

It has two links, joined at Stage 1's output contract:

- **link 1** — ``stage1_ingest.run_ingest`` on the tiny committed
  ``tests/fixtures/`` (the *real* Stage 1 code path, acquisition layer
  patched to the fixtures): ``genes`` / ``gene_sets`` / ``interactions`` /
  ``ot_disease_subset`` / ``manifest.json`` all validate; every gene id is
  an Ensembl id; ``source_db == "reactome"`` on every normalized row.

- **link 2** — Stages 2→6 via ``run_pipeline`` on the shared 17-gene
  synthetic Stage 1 output the Stage 2-6 unit tests are built around
  (``example_data/example_run``, copied into ``tmp_path``). blitzgsea
  (Stage 2) cannot calibrate a GSEA on the 3-gene signature the
  ``tests/fixtures/`` inputs produce, so the chain continues from that
  larger synthetic Stage 1 output — which is emitted and schema-validated
  exactly like ``stage1_ingest``'s. Every remaining hand-off
  (``disease_pathways`` → graph artifacts → ``gene_weights`` →
  ``candidate_scores`` → ``candidates_annotated``) validates; Ensembl ids
  and ``source_db`` survive to the final table; and the final table carries
  a ``composite_score`` in [0, 1] for **every** graph gene (candidates ==
  graph nodes − target, none missing).
"""

from __future__ import annotations

import json
import shutil

import pandas as pd
import pytest
import scipy.sparse as sp

import run_pipeline
import stage1_ingest
from conftest import make_args
from stage0_schemas import (
    ENSEMBL_GENE_ID_PATTERN,
    AnnotatedCandidatesSchema,
    CandidateScoresSchema,
    DiseasePathwaysSchema,
    GeneSetsSchema,
    GeneWeightsSchema,
    GenesSchema,
    InteractionsSchema,
    Manifest,
    OTAssociationsSchema,
    assert_foreign_key,
)

REPO = run_pipeline.PIPELINE_DIR.parent
EXAMPLE_RUN = REPO / "example_data" / "example_run"
STAGE1_FILES = [
    "manifest.json", "genes.parquet", "gene_sets.parquet", "interactions.parquet",
    "ot_disease_subset.parquet", "ot_target_subset.parquet",
]


def _all_ensembl(ids) -> bool:
    return all(ENSEMBL_GENE_ID_PATTERN.match(str(g)) for g in ids)


@pytest.mark.skipif(
    not (EXAMPLE_RUN / "manifest.json").exists(),
    reason="example_data/example_run not populated (run example_data/build_gsea_example.py)",
)
def test_pipeline_schema_conformance_end_to_end(tmp_path, patched_acquisition):
    # ================================================================
    # link 1 — the real Stage 1 code path (tests/fixtures/)
    # ================================================================
    s1 = tmp_path / "s1"
    stage1_ingest.run_ingest(make_args(out_dir=s1, data_dir=tmp_path / "s1_data"))

    genes = pd.read_parquet(s1 / "genes.parquet")
    gene_sets = pd.read_parquet(s1 / "gene_sets.parquet")
    interactions = pd.read_parquet(s1 / "interactions.parquet")
    ot = pd.read_parquet(s1 / "ot_disease_subset.parquet")
    GenesSchema.validate(genes)
    GeneSetsSchema.validate(gene_sets)
    InteractionsSchema.validate(interactions)
    OTAssociationsSchema.validate(ot)
    manifest = Manifest.model_validate(json.loads((s1 / "manifest.json").read_text()))

    # foreign keys: every downstream gene id resolves against genes.parquet
    assert_foreign_key(gene_sets, "gene_id", genes)
    assert_foreign_key(interactions, "gene_a", genes)
    assert_foreign_key(interactions, "gene_b", genes)
    assert_foreign_key(ot, "gene_id", genes)

    assert _all_ensembl(genes["gene_id"])
    assert set(gene_sets["source_db"]) == {"reactome"}
    assert set(interactions["source_db"]) == {"reactome"}
    assert manifest.resolved_target.gene_id in set(genes["gene_id"])

    # ================================================================
    # link 2 — Stages 2→6 from the 17-gene synthetic Stage 1 output
    # ================================================================
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    for name in STAGE1_FILES:
        shutil.copy(EXAMPLE_RUN / name, run_dir / name)

    # the Stage 1 output contract Stage 2 reads
    genes = pd.read_parquet(run_dir / "genes.parquet")
    GenesSchema.validate(genes)
    GeneSetsSchema.validate(pd.read_parquet(run_dir / "gene_sets.parquet"))
    InteractionsSchema.validate(pd.read_parquet(run_dir / "interactions.parquet"))
    OTAssociationsSchema.validate(pd.read_parquet(run_dir / "ot_disease_subset.parquet"))
    manifest = Manifest.model_validate(json.loads((run_dir / "manifest.json").read_text()))
    target_id = manifest.resolved_target.gene_id
    gene_universe = set(genes["gene_id"])
    assert _all_ensembl(gene_universe)

    run_pipeline.run_pipeline(from_manifest=run_dir / "manifest.json")

    # --- Stage 2 hand-off ---
    disease_pathways = pd.read_csv(run_dir / "disease_pathways.tsv", sep="\t")
    DiseasePathwaysSchema.validate(disease_pathways)
    assert set(disease_pathways["source_db"]) <= {"reactome"}  # source_db attribution preserved

    # --- Stage 3 hand-off (sparse serialization, not a pandera schema) ---
    graph_index = json.loads((run_dir / "graph_gene_index.json").read_text())
    assert _all_ensembl(graph_index)
    assert set(graph_index) <= gene_universe  # graph nodes are canonical genes, no leakage
    weight = sp.load_npz(run_dir / "graph_weight.npz")
    sign = sp.load_npz(run_dir / "graph_sign.npz")
    n = len(graph_index)
    assert weight.shape == sign.shape == (n, n)
    graph_meta = json.loads((run_dir / "graph_metadata.json").read_text())
    assert graph_meta["resolved_target"]["gene_id"] == target_id

    # --- Stage 4 hand-off ---
    gene_weights = pd.read_csv(run_dir / "gene_weights.tsv", sep="\t")
    GeneWeightsSchema.validate(gene_weights)
    assert _all_ensembl(gene_weights["gene_id"])
    assert "stage4_genetic_evidence" in json.loads((run_dir / "graph_metadata.json").read_text())

    # --- Stage 5 hand-off ---
    candidate_scores = pd.read_csv(run_dir / "candidate_scores.tsv", sep="\t")
    CandidateScoresSchema.validate(candidate_scores)
    assert _all_ensembl(candidate_scores["gene"])
    assert set(candidate_scores["gene"]) <= gene_universe

    # --- Stage 6 hand-off — the final output ---
    annotated = pd.read_csv(run_dir / "candidates_annotated.tsv", sep="\t")
    AnnotatedCandidatesSchema.validate(annotated)
    assert _all_ensembl(annotated["gene"])
    assert set(annotated["gene"]) <= gene_universe

    # a composite_score for every input gene: candidates == graph nodes − target
    assert set(annotated["gene"]) == set(graph_index) - {target_id}
    assert annotated["composite_score"].notna().all()
    assert annotated["composite_score"].between(0.0, 1.0).all()
