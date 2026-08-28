"""Inter-stage data contracts for the targets-from-pathways v2 pipeline.

Every stage after Stage 1 reads its inputs from, and writes its outputs as,
one of the artifacts defined here — never a raw source file, and never a
column layout described only in a docstring or comment. Concretely:

- Tabular artifacts (``genes.parquet``, ``gene_sets.parquet``,
  ``interactions.parquet``, ``ot_disease_subset.parquet``, and Stage 2's
  ``disease_pathways.tsv``) are pandas DataFrames validated against the
  ``pandera.DataFrameSchema`` objects below via ``SCHEMA.validate(df)``, on
  read *and* on write. ``disease_pathways.tsv`` is TSV rather than parquet
  (it's meant to be eyeballed), but the contract discipline is the same:
  Stage 3 and Stage 7 read it as ``DiseasePathwaysSchema``, never by
  re-deriving its column layout from ``stage2_gsea_discovery.py``'s source.
- ``manifest.json`` is validated against the ``Manifest`` pydantic model via
  ``Manifest.model_validate(json.load(f))`` on read and
  ``Manifest.model_dump(mode="json")`` on write.

``source_db`` is retained as a column on ``GeneSetsSchema`` and
``InteractionsSchema`` even though this hackathon build only ever emits
``"reactome"`` (see ``SUPPORTED_SOURCE_DBS``). This is deliberate: it is the
seam a second curation (e.g. a future WikiPathways or OmniPath adapter) will
need to plug into downstream stages without a refactor. Do not add KEGG,
WikiPathways, or OmniPath support in this version — extend
``SUPPORTED_SOURCE_DBS`` and the corresponding ingest adapter only when that
work is actually scoped.

Gene identity is Ensembl gene ID throughout (``genes.gene_id`` is the single
primary key every other table's ``gene_id``/``gene_a``/``gene_b`` foreign
keys resolve against); HGNC symbols are carried only for display. This is
enforced at the column level via ``ENSEMBL_GENE_ID_PATTERN`` — v1 leaked a
Reactome physical-entity ID (``GRB2-1``) into a gene list by skipping exactly
this check.
"""

from __future__ import annotations

import re
from datetime import datetime
from typing import Any

import numpy as np
import pandas as pd
from pandera.pandas import Check, Column, DataFrameSchema
from pydantic import BaseModel, Field, field_validator

__all__ = [
    "ENSEMBL_GENE_ID_PATTERN",
    "SHA256_PATTERN",
    "SUPPORTED_SOURCE_DBS",
    "SUPPORTED_PROVENANCE_DBS",
    "EVIDENCE_TYPES",
    "SIGN_VALUES",
    "GenesSchema",
    "GeneSetsSchema",
    "InteractionsSchema",
    "OTAssociationsSchema",
    "DiseasePathwaysSchema",
    "GeneWeightsSchema",
    "CandidateScoresSchema",
    "AnnotatedCandidatesSchema",
    "TRACTABILITY_BUCKETS",
    "SAFETY_FLAGS",
    "assert_foreign_key",
    "ResolvedTarget",
    "ResolvedDisease",
    "SourceFile",
    "Source",
    "OutputArtifact",
    "CoverageEntry",
    "GeneSetSizeBucket",
    "ScaleReport",
    "Manifest",
]

# --------------------------------------------------------------------------
# Shared vocabularies
# --------------------------------------------------------------------------

# Ensembl human gene IDs: "ENSG" + 11 digits. Rejects HGNC symbols and
# Reactome physical-entity identifiers (e.g. "GRB2-1"), which is the exact
# class of bug this pattern exists to catch (see module docstring).
ENSEMBL_GENE_ID_PATTERN = re.compile(r"^ENSG\d{11}$")

SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")

# Extend this set — and only this set — when a second curation is added.
SUPPORTED_SOURCE_DBS = frozenset({"reactome"})

# Manifest.sources.db is a broader vocabulary than SUPPORTED_SOURCE_DBS: the
# latter gates the *pathway-curation* identity carried on GeneSetsSchema/
# InteractionsSchema.source_db (reactome, and only reactome, for this
# hackathon build). Stage 1's manifest additionally records acquisition
# provenance for Open Targets (target/disease/association parquet), which is
# a genetic-evidence source, not a pathway curation — so it is admitted here
# without being added to SUPPORTED_SOURCE_DBS or to either DataFrameSchema.
SUPPORTED_PROVENANCE_DBS = SUPPORTED_SOURCE_DBS | {"opentargets"}

EVIDENCE_TYPES = frozenset({"curated", "predicted"})

SIGN_VALUES = frozenset({-1, 0, 1})


def _is_list_of_str(series: pd.Series) -> pd.Series:
    # A DataFrame built in-process holds Python `list`; one that has round-
    # tripped through parquet (pyarrow's LIST type -> pandas) comes back as
    # `numpy.ndarray` of the same `str` elements, not `list` — verified
    # against a real genes.parquet written by stage1_ingest.py and re-read with a
    # bare `pd.read_parquet`. Both container types are accepted here so
    # "validate on read" (module docstring) doesn't spuriously fail for
    # every stage after Stage 1; a bare string (the real malformed case this
    # check exists to catch) is still neither and is still rejected.
    return series.apply(
        lambda v: isinstance(v, (list, np.ndarray)) and all(isinstance(x, str) for x in v)
    )


# --------------------------------------------------------------------------
# Tabular artifacts (pandera)
# --------------------------------------------------------------------------

GenesSchema = DataFrameSchema(
    {
        "gene_id": Column(
            str,
            checks=Check.str_matches(ENSEMBL_GENE_ID_PATTERN.pattern),
            unique=True,
            nullable=False,
        ),
        "symbol": Column(str, nullable=False),
        "synonyms": Column(
            object,
            checks=Check(_is_list_of_str, error="synonyms must be a list[str]"),
            nullable=False,
        ),
        "biotype": Column(str, nullable=False),
    },
    strict=True,
)

GeneSetsSchema = DataFrameSchema(
    {
        "set_id": Column(str, nullable=False),
        "set_name": Column(str, nullable=False),
        "source_db": Column(str, checks=Check.isin(SUPPORTED_SOURCE_DBS), nullable=False),
        "source_version": Column(str, nullable=False),
        "gene_id": Column(
            str,
            checks=Check.str_matches(ENSEMBL_GENE_ID_PATTERN.pattern),
            nullable=False,
        ),
        # coerce=True here only: a NaN-bearing int column round-trips through
        # pandas/parquet as float64, so this column alone needs coercion to
        # the nullable Int64 dtype. Every other column keeps strict dtype
        # checking — no schema-wide coerce, since that would silently accept
        # e.g. a non-bool "directed" column (see InteractionsSchema).
        "hierarchy_level": Column(pd.Int64Dtype(), nullable=True, coerce=True),
        "parent_id": Column(str, nullable=True),
    },
    unique=["set_id", "gene_id"],  # one row per (set, gene) pair
    strict=True,
)

InteractionsSchema = DataFrameSchema(
    {
        "gene_a": Column(
            str,
            checks=Check.str_matches(ENSEMBL_GENE_ID_PATTERN.pattern),
            nullable=False,
        ),
        "gene_b": Column(
            str,
            checks=Check.str_matches(ENSEMBL_GENE_ID_PATTERN.pattern),
            nullable=False,
        ),
        "directed": Column(bool, nullable=False),
        "sign": Column(int, checks=Check.isin(SIGN_VALUES), nullable=False),
        "source_db": Column(str, checks=Check.isin(SUPPORTED_SOURCE_DBS), nullable=False),
        "source_version": Column(str, nullable=False),
        "evidence_type": Column(str, checks=Check.isin(EVIDENCE_TYPES), nullable=False),
        # coerce=True: same int64-vs-float64 round-trip rationale as
        # hierarchy_level above. Known residual gap (same class as the
        # "directed" trap, narrower blast radius): a genuinely bool-dtype
        # column here would silently cast to 0.0/1.0 and pass in_range(0, 1)
        # rather than raising, since bool->float64 is a lossless numpy cast.
        # Accepted for now — Reactome FI confidence is always numeric in
        # practice — but pinned by
        # test_interactions_confidence_bool_column_coerces_silently so a
        # future change in this trade-off is a deliberate decision, not a
        # surprise.
        "confidence": Column(float, checks=Check.in_range(0.0, 1.0), nullable=True, coerce=True),
    },
    strict=True,
)

OTAssociationsSchema = DataFrameSchema(
    {
        "gene_id": Column(
            str,
            checks=Check.str_matches(ENSEMBL_GENE_ID_PATTERN.pattern),
            nullable=False,
        ),
        "disease_id": Column(str, nullable=False),
        "datatype_id": Column(str, nullable=False),
        "datasource_id": Column(str, nullable=False),
        # coerce=True: same int64-vs-float64 round-trip rationale as
        # GeneSetsSchema.hierarchy_level and InteractionsSchema.confidence
        # above, and the same accepted bool->float64 residual gap.
        "score": Column(float, checks=Check.in_range(0.0, 1.0), nullable=False, coerce=True),
    },
    strict=True,
)


DiseasePathwaysSchema = DataFrameSchema(
    {
        # Reactome stable ID (or equivalent for a future source_db) — the
        # join key back to gene_sets.parquet for Stage 3's graph construction.
        "set_id": Column(str, nullable=False, unique=True),
        # Human-readable pathway name (gene_sets.set_name), for display.
        "gene_set": Column(str, nullable=False),
        "source_db": Column(str, checks=Check.isin(SUPPORTED_SOURCE_DBS), nullable=False),
        # Nominal GSEA p-value, pre-correction.
        "pval": Column(float, checks=Check.in_range(0.0, 1.0), nullable=False),
        # Benjamini-Hochberg FDR, corrected across every pathway actually
        # tested in this run (see stage2_gsea_discovery.py's module docstring) —
        # not comparable across runs with a different tested-pathway set.
        "fdr": Column(float, checks=Check.in_range(0.0, 1.0), nullable=False),
        # Whether this pathway's (collapsed) gene membership includes the
        # manifest's resolved --target gene.
        "contains_target": Column(bool, nullable=False),
    },
    strict=True,
)


GeneWeightsSchema = DataFrameSchema(
    {
        "gene_id": Column(
            str,
            checks=Check.str_matches(ENSEMBL_GENE_ID_PATTERN.pattern),
            unique=True,
            nullable=False,
        ),
        # Non-pathway-datatype Open Targets evidence (see
        # stage4_genetic_evidence_weights.py's module docstring) — 0.0, not null,
        # for a gene with no matching evidence, so every graph node gets a
        # row here.
        "genetic_evidence_score": Column(float, checks=Check.in_range(0.0, 1.0), nullable=False),
    },
    strict=True,
)


CandidateScoresSchema = DataFrameSchema(
    {
        # Ensembl gene ID of a candidate alternative target — every non-target
        # node of the Stage 3/4 graph gets exactly one row. Same identifier
        # namespace as GenesSchema.gene_id / OTAssociationsSchema.gene_id
        # (all Open-Targets-canonicalised in Stage 1), so Stage 6 joins this
        # column straight onto Open Targets annotation with no ID mapping;
        # the pattern check here is what makes a namespace slip fail loudly.
        "gene": Column(
            str,
            checks=Check.str_matches(ENSEMBL_GENE_ID_PATTERN.pattern),
            unique=True,
            nullable=False,
        ),
        # Direction- and weight-aware topology proximity to the target
        # (higher = closer); 0.0 for a node the target cannot reach. Always
        # present — topology is Stage 5's always-on method.
        "topology_score": Column(
            float, checks=Check.greater_than_or_equal_to(0.0), nullable=False, coerce=True
        ),
        # Personalised random-walk-with-restart stationary probability,
        # restart on the target. Present only when Stage 5 was run with
        # `rwr` in --method (hence required=False); a probability in [0, 1].
        "rwr_score": Column(
            float, checks=Check.in_range(0.0, 1.0), nullable=False, coerce=True, required=False
        ),
    },
    strict=True,
)


# Stage 6 vocabularies. Tractability is an ordinal proxy derived from Open
# Targets' per-modality flags; safety has deliberately NO "safe" value — OT
# supplies known *liabilities* or silence, never positive proof of safety,
# so a gene with no liabilities annotation is "unknown", not "safe".
TRACTABILITY_BUCKETS = frozenset({"clinical", "discovery", "unknown"})
SAFETY_FLAGS = frozenset({"has_liabilities", "unknown"})


AnnotatedCandidatesSchema = DataFrameSchema(
    {
        "gene": Column(
            str,
            checks=Check.str_matches(ENSEMBL_GENE_ID_PATTERN.pattern),
            unique=True,
            nullable=False,
        ),
        # --- carried through from Stage 5 (schemas.CandidateScoresSchema) ---
        "topology_score": Column(
            float, checks=Check.greater_than_or_equal_to(0.0), nullable=False, coerce=True
        ),
        "rwr_score": Column(
            float, checks=Check.in_range(0.0, 1.0), nullable=False, coerce=True, required=False
        ),
        # --- joined from Stage 4's graph_metadata.json (absent if Stage 4
        #     was not run) ---
        "genetic_evidence_score": Column(
            float, checks=Check.in_range(0.0, 1.0), nullable=False, coerce=True, required=False
        ),
        # --- Open Targets target-level annotation (this stage) ---
        "tractability": Column(str, checks=Check.isin(TRACTABILITY_BUCKETS), nullable=False),
        "safety": Column(str, checks=Check.isin(SAFETY_FLAGS), nullable=False),
        # Count of Open Targets safetyLiabilities entries; 0 is "none recorded"
        # (which maps to safety == "unknown", never "safe").
        "n_safety_liabilities": Column(
            int, checks=Check.greater_than_or_equal_to(0), nullable=False, coerce=True
        ),
        # Weighted average of the normalised components, in [0, 1].
        "composite_score": Column(float, checks=Check.in_range(0.0, 1.0), nullable=False, coerce=True),
        # Human-readable "k=v|k=v" of the normalised component values that
        # fed composite_score — the evidence trace Stage 8's report needs.
        "composite_breakdown": Column(str, nullable=False),
    },
    strict=True,
)


def assert_foreign_key(
    child_df: pd.DataFrame, child_col: str, genes_df: pd.DataFrame, *, gene_col: str = "gene_id"
) -> None:
    """Raise ValueError if ``child_df[child_col]`` contains an ID absent from ``genes_df[gene_col]``.

    Pandera validates one DataFrame at a time and cannot express a foreign
    key across two tables, so this small helper exists for stages to check
    ``gene_sets.gene_id``, ``interactions.gene_a``/``gene_b``, and
    ``ot_disease_subset.gene_id`` against ``genes.gene_id`` after validating
    each table independently.
    """
    known = set(genes_df[gene_col])
    missing = sorted(set(child_df[child_col]) - known)
    if missing:
        raise ValueError(
            f"{child_col} contains {len(missing)} gene_id(s) not present in genes.gene_id: "
            f"{missing[:20]}"
        )


# --------------------------------------------------------------------------
# Manifest (pydantic)
# --------------------------------------------------------------------------


class ResolvedTarget(BaseModel):
    input: str
    gene_id: str
    symbol: str

    @field_validator("gene_id")
    @classmethod
    def _valid_ensembl_id(cls, v: str) -> str:
        if not ENSEMBL_GENE_ID_PATTERN.match(v):
            raise ValueError(f"gene_id {v!r} is not a valid Ensembl gene ID")
        return v


class ResolvedDisease(BaseModel):
    efo_id: str
    name: str
    n_associated_genes: int = Field(ge=0)


class SourceFile(BaseModel):
    path: str
    sha256: str
    bytes: int = Field(ge=0)

    @field_validator("sha256")
    @classmethod
    def _valid_sha256(cls, v: str) -> str:
        if not SHA256_PATTERN.match(v):
            raise ValueError(f"sha256 {v!r} is not 64 lowercase hex characters")
        return v


class Source(BaseModel):
    db: str
    version: str
    files: list[SourceFile]

    @field_validator("db")
    @classmethod
    def _supported_db(cls, v: str) -> str:
        if v not in SUPPORTED_PROVENANCE_DBS:
            raise ValueError(f"db {v!r} is not in SUPPORTED_PROVENANCE_DBS {sorted(SUPPORTED_PROVENANCE_DBS)}")
        return v


class OutputArtifact(BaseModel):
    path: str
    sha256: str

    @field_validator("sha256")
    @classmethod
    def _valid_sha256(cls, v: str) -> str:
        if not SHA256_PATTERN.match(v):
            raise ValueError(f"sha256 {v!r} is not 64 lowercase hex characters")
        return v


class CoverageEntry(BaseModel):
    """One entry of coverage_report.json's per-source mapping summary."""

    source: str
    raw_gene_count: int = Field(ge=0)
    mapped_count: int = Field(ge=0)
    percent_mapped: float = Field(ge=0.0, le=100.0)
    dropped_count: int = Field(ge=0)
    example_dropped_ids: list[str] = Field(default_factory=list, max_length=20)

    @field_validator("mapped_count")
    @classmethod
    def _mapped_le_raw(cls, v: int, info) -> int:
        raw = info.data.get("raw_gene_count")
        if raw is not None and v > raw:
            raise ValueError(f"mapped_count ({v}) cannot exceed raw_gene_count ({raw})")
        return v

    @field_validator("percent_mapped")
    @classmethod
    def _percent_matches_counts(cls, v: float, info) -> float:
        # Catches the class of bug where percent_mapped is computed
        # separately from mapped_count/raw_gene_count and silently drifts
        # out of sync — which matters here because the 90% coverage floor
        # gating the whole pipeline run (see README.md) reads percent_mapped
        # directly. Tolerance of 0.51 accommodates integer-rounded reporting.
        raw = info.data.get("raw_gene_count")
        mapped = info.data.get("mapped_count")
        if raw is not None and mapped is not None:
            expected = (100.0 * mapped / raw) if raw > 0 else 0.0
            if abs(v - expected) > 0.51:
                raise ValueError(
                    f"percent_mapped ({v}) is inconsistent with "
                    f"mapped_count/raw_gene_count ({mapped}/{raw} = {expected:.2f})"
                )
        return v

    @field_validator("dropped_count")
    @classmethod
    def _dropped_matches_counts(cls, v: int, info) -> int:
        # Unlike percent_mapped, these are integer counts, so the identity
        # is exact — no rounding tolerance needed.
        raw = info.data.get("raw_gene_count")
        mapped = info.data.get("mapped_count")
        if raw is not None and mapped is not None and v != raw - mapped:
            raise ValueError(
                f"dropped_count ({v}) must equal raw_gene_count - mapped_count "
                f"({raw} - {mapped} = {raw - mapped})"
            )
        return v


class GeneSetSizeBucket(BaseModel):
    set_size: int = Field(ge=1)
    count: int = Field(ge=0)


class ScaleReport(BaseModel):
    """scale_report.json: gene-set size distribution and the projected Stage 3 edge count."""

    gene_set_size_distribution_before_cap: list[GeneSetSizeBucket]
    gene_set_size_distribution_after_cap: list[GeneSetSizeBucket]
    sets_retained: int = Field(ge=0)
    interaction_counts_by_sign: dict[str, int]
    interaction_counts_by_evidence_type: dict[str, int]
    projected_comembership_edge_count: int = Field(ge=0)


class Manifest(BaseModel):
    run_id: str
    git_commit: str
    created_at: datetime
    seed: int
    resolved_target: ResolvedTarget
    disease: ResolvedDisease
    sources: list[Source]
    cli_parameters: dict[str, Any]
    output_artifacts: list[OutputArtifact]
    coverage_report: dict[str, CoverageEntry]
    scale_report: ScaleReport
