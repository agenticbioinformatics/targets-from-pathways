"""Stage 0 - data contracts shared by all stages, plus small helpers every stage needs.

Tables are validated with pandera on read and on write; manifest.json is validated with pydantic.
"""

import hashlib
import json
import logging
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd
import scipy.sparse as sp
from pandera.pandas import Check, Column, DataFrameSchema
from pydantic import BaseModel, Field, field_validator, model_validator

# Gene identity is the Ensembl gene ID everywhere. The pattern rejects HGNC symbols and Reactome
# physical-entity IDs such as "GRB2-1", which once leaked into gene lists.
ENSEMBL_GENE_ID_PATTERN = re.compile(r"^ENSG\d{11}$")
SHA256_PATTERN = r"^[0-9a-f]{64}$"

# Pathway databases allowed in gene_sets/interactions. Kept as a column so a second database
# can be added later without changing downstream stages.
SUPPORTED_SOURCE_DBS = {"reactome"}
# Databases allowed in manifest.sources: Open Targets is a data source but not a pathway database.
SUPPORTED_PROVENANCE_DBS = {"reactome", "opentargets"}
EVIDENCE_TYPES = {"curated", "predicted"}
SIGN_VALUES = {-1, 0, 1}
TRACTABILITY_BUCKETS = {"clinical", "discovery", "unknown"}
# There is no "safe" value: Open Targets lists known liabilities or nothing, never proof of safety.
SAFETY_FLAGS = {"has_liabilities", "unknown"}


def is_list_of_strings(series):
    # After a parquet round-trip, list columns come back as numpy arrays instead of Python lists.
    results = []
    for value in series:
        ok = isinstance(value, (list, np.ndarray))
        if ok:
            for item in value:
                if not isinstance(item, str):
                    ok = False
        results.append(ok)
    return pd.Series(results, index=series.index)


def ensembl_column(unique=False):
    return Column(str, checks=Check.str_matches(ENSEMBL_GENE_ID_PATTERN.pattern), unique=unique, nullable=False)


# coerce=True is used only on numeric columns: a parquet round-trip can turn int into float.
# Other columns are not coerced, so e.g. a non-bool "directed" column still fails.
GenesSchema = DataFrameSchema(
    {
        "gene_id": ensembl_column(unique=True),
        "symbol": Column(str, nullable=False),
        "synonyms": Column(object, checks=Check(is_list_of_strings, error="synonyms must be a list[str]"), nullable=False),
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
        "gene_id": ensembl_column(),
        "hierarchy_level": Column(pd.Int64Dtype(), nullable=True, coerce=True),
        "parent_id": Column(str, nullable=True),
    },
    unique=["set_id", "gene_id"],
    strict=True,
)

InteractionsSchema = DataFrameSchema(
    {
        "gene_a": ensembl_column(),
        "gene_b": ensembl_column(),
        "directed": Column(bool, nullable=False),
        "sign": Column(int, checks=Check.isin(SIGN_VALUES), nullable=False),
        "source_db": Column(str, checks=Check.isin(SUPPORTED_SOURCE_DBS), nullable=False),
        "source_version": Column(str, nullable=False),
        "evidence_type": Column(str, checks=Check.isin(EVIDENCE_TYPES), nullable=False),
        # Known gap, pinned by a test: a bool column would be coerced to 0.0/1.0 and pass.
        "confidence": Column(float, checks=Check.in_range(0.0, 1.0), nullable=True, coerce=True),
    },
    strict=True,
)

OTAssociationsSchema = DataFrameSchema(
    {
        "gene_id": ensembl_column(),
        "disease_id": Column(str, nullable=False),
        "datatype_id": Column(str, nullable=False),
        "datasource_id": Column(str, nullable=False),
        "score": Column(float, checks=Check.in_range(0.0, 1.0), nullable=False, coerce=True),
    },
    strict=True,
)

DiseasePathwaysSchema = DataFrameSchema(
    {
        "set_id": Column(str, nullable=False, unique=True),
        "gene_set": Column(str, nullable=False),
        "source_db": Column(str, checks=Check.isin(SUPPORTED_SOURCE_DBS), nullable=False),
        # Coerced so that an empty table (no pathway could be tested) read back from TSV still validates.
        "pval": Column(float, checks=Check.in_range(0.0, 1.0), nullable=False, coerce=True),
        # FDR is corrected over the pathways tested in this run only, so it is not comparable across runs.
        "fdr": Column(float, checks=Check.in_range(0.0, 1.0), nullable=False, coerce=True),
    },
    strict=True,
)

# Stage 2's seed genes: the disease's genes with a genetic_association score > 0 (minus any benchmark holdout),
# with their strongest genetic_association score. GSEA ranks them (Stage 2) and Stage 5 propagates from them.
# score is optional so that a hand-made seed list (gene_id only) can be passed to Stage 5.
SeedGenesSchema = DataFrameSchema(
    {
        "gene_id": ensembl_column(unique=True),
        "score": Column(float, checks=Check.in_range(0.0, 1.0), nullable=False, coerce=True, required=False),
    },
    strict=True,
)

# Stage 2's benchmark holdout, resolved to Ensembl IDs (empty when there is none). Stage 4 must not use these
# genes' own evidence, or a held-out gene would be ranked by the very association it is meant to recover.
HoldoutGenesSchema = DataFrameSchema(
    {"gene_id": ensembl_column(unique=True)},
    strict=True,
)

GeneticEvidenceScoresSchema = DataFrameSchema(
    {
        "gene_id": ensembl_column(unique=True),
        "genetic_evidence_score": Column(float, checks=Check.in_range(0.0, 1.0), nullable=False),
    },
    strict=True,
)

# Stage 4's optional user-supplied weights (--node-weights / --edge-weights). A weight must be > 0: a zero
# would silently cut the edge out of propagation, and Stage 3's graph has no way to say "no edge" but leaving it out.
def positive_weight_column():
    return Column(float, checks=[Check.greater_than(0.0), Check(np.isfinite, element_wise=True)],
                  nullable=False, coerce=True)


UserNodeWeightsSchema = DataFrameSchema(
    {"gene_id": ensembl_column(unique=True), "weight": positive_weight_column()},
    strict=True,
)
UserEdgeWeightsSchema = DataFrameSchema(
    {"gene_a": ensembl_column(), "gene_b": ensembl_column(), "weight": positive_weight_column()},
    strict=True,
)

# Stage 5 writes a score column per requested --method; at least one method is always requested.
PROPAGATION_SCORE_COLUMNS = ["rwr_score", "bfwalk_score"]


def has_a_propagation_score(df):
    for column in PROPAGATION_SCORE_COLUMNS:
        if column in df.columns:
            return True
    return False


HAS_PROPAGATION_SCORE = Check(has_a_propagation_score, error=f"needs at least one of {PROPAGATION_SCORE_COLUMNS}")


def ranking_column(columns):
    """The score candidates are ranked by: the last of PROPAGATION_SCORE_COLUMNS present (bfwalk over rwr)."""
    present = [column for column in PROPAGATION_SCORE_COLUMNS if column in columns]
    if not present:
        raise ValueError(f"no propagation score column; expected one of {PROPAGATION_SCORE_COLUMNS}")
    return present[-1]

CandidateScoresSchema = DataFrameSchema(
    {
        "gene": ensembl_column(unique=True),
        "rwr_score": Column(float, checks=Check.in_range(0.0, 1.0), nullable=False, coerce=True, required=False),
        "bfwalk_score": Column(float, checks=Check.greater_than_or_equal_to(0.0), nullable=False, coerce=True, required=False),
    },
    checks=HAS_PROPAGATION_SCORE,
    strict=True,
)

AnnotatedCandidatesSchema = DataFrameSchema(
    {
        "gene": ensembl_column(unique=True),
        "rwr_score": Column(float, checks=Check.in_range(0.0, 1.0), nullable=False, coerce=True, required=False),
        "bfwalk_score": Column(float, checks=Check.greater_than_or_equal_to(0.0), nullable=False, coerce=True, required=False),
        "genetic_evidence_score": Column(
            float, checks=Check.in_range(0.0, 1.0), nullable=False, coerce=True, required=False
        ),
        "tractability": Column(str, checks=Check.isin(TRACTABILITY_BUCKETS), nullable=False),
        "safety": Column(str, checks=Check.isin(SAFETY_FLAGS), nullable=False),
        "n_safety_liabilities": Column(int, checks=Check.greater_than_or_equal_to(0), nullable=False, coerce=True),
    },
    checks=HAS_PROPAGATION_SCORE,
    strict=True,
)


def assert_foreign_key(child_df, child_column, genes_df, gene_column="gene_id"):
    """Raise ValueError if child_df[child_column] has gene IDs missing from genes_df (pandera checks one table only)."""
    known_genes = set(genes_df[gene_column])
    missing = sorted(set(child_df[child_column]) - known_genes)
    if missing:
        raise ValueError(f"{child_column} contains {len(missing)} gene_id(s) not present in genes.gene_id: {missing[:20]}")


class ResolvedDisease(BaseModel):
    efo_id: str
    name: str
    n_associated_genes: int = Field(ge=0)


class SourceFile(BaseModel):
    path: str
    sha256: str = Field(pattern=SHA256_PATTERN)
    bytes: int = Field(ge=0)


class Source(BaseModel):
    db: str
    version: str
    files: List[SourceFile]

    @field_validator("db")
    @classmethod
    def check_db(cls, value):
        if value not in SUPPORTED_PROVENANCE_DBS:
            raise ValueError(f"db {value!r} is not in SUPPORTED_PROVENANCE_DBS {sorted(SUPPORTED_PROVENANCE_DBS)}")
        return value


class OutputArtifact(BaseModel):
    path: str
    sha256: str = Field(pattern=SHA256_PATTERN)


class CoverageEntry(BaseModel):
    source: str
    raw_gene_count: int = Field(ge=0)
    mapped_count: int = Field(ge=0)
    percent_mapped: float = Field(ge=0.0, le=100.0)
    dropped_count: int = Field(ge=0)
    example_dropped_ids: List[str] = Field(default_factory=list, max_length=20)

    # The 90% coverage floor reads percent_mapped, so it must agree with the counts.
    @model_validator(mode="after")
    def check_counts_agree(self):
        raw = self.raw_gene_count
        mapped = self.mapped_count
        if mapped > raw:
            raise ValueError(f"mapped_count ({mapped}) cannot exceed raw_gene_count ({raw})")
        expected_percent = 0.0
        if raw > 0:
            expected_percent = 100.0 * mapped / raw
        # 0.51 tolerance allows for rounding.
        if abs(self.percent_mapped - expected_percent) > 0.51:
            raise ValueError(
                f"percent_mapped ({self.percent_mapped}) is inconsistent with "
                f"mapped_count/raw_gene_count ({mapped}/{raw} = {expected_percent:.2f})"
            )
        if self.dropped_count != raw - mapped:
            raise ValueError(
                f"dropped_count ({self.dropped_count}) must equal raw_gene_count - mapped_count ({raw} - {mapped})"
            )
        return self


class GeneSetSizeBucket(BaseModel):
    set_size: int = Field(ge=1)
    count: int = Field(ge=0)


class ScaleReport(BaseModel):
    gene_set_size_distribution_before_cap: List[GeneSetSizeBucket]
    gene_set_size_distribution_after_cap: List[GeneSetSizeBucket]
    sets_retained: int = Field(ge=0)
    interaction_counts_by_sign: Dict[str, int]
    interaction_counts_by_evidence_type: Dict[str, int]


class Manifest(BaseModel):
    run_id: str
    git_commit: str
    created_at: datetime
    seed: int
    disease: ResolvedDisease
    sources: List[Source]
    cli_parameters: Dict[str, Any]
    output_artifacts: List[OutputArtifact]
    coverage_report: Dict[str, CoverageEntry]
    scale_report: ScaleReport


def fail(message):
    logging.getLogger("pipeline").error(message)
    sys.exit(1)


def load_manifest(manifest_path):
    return Manifest.model_validate(json.loads(Path(manifest_path).read_text()))


def sha256_of_file(path):
    hasher = hashlib.sha256()
    with open(path, "rb") as file:
        while True:
            chunk = file.read(1024 * 1024)
            if not chunk:
                break
            hasher.update(chunk)
    return hasher.hexdigest()


def find_artifact(manifest, manifest_path, filename):
    """Return the path of a Stage 1 output file listed in the manifest. The file next to the manifest wins:
    the listed path is relative to whatever directory Stage 1 ran in, and in a copied or moved run directory
    it would point at another run's file. The file's sha256 must match the manifest."""
    for artifact in manifest.output_artifacts:
        listed_path = Path(artifact.path)
        if listed_path.name != filename:
            continue
        next_to_manifest = Path(manifest_path).parent / filename
        if next_to_manifest.exists():
            path = next_to_manifest
        elif listed_path.exists():
            path = listed_path
        else:
            fail(f"Manifest artifact {filename!r} not found next to {manifest_path} (also checked {listed_path}).")
        if sha256_of_file(path) != artifact.sha256:
            fail(f"{path} does not match the sha256 recorded in {manifest_path} - it was changed after Stage 1 "
                 f"or belongs to another run.")
        return path
    fail(f"Manifest {manifest_path} lists no output artifact named {filename!r}.")


def load_graph(graph_dir):
    """Load the graph files written by Stage 3 (or Stage 4, which uses the same format)."""
    graph_dir = Path(graph_dir)
    weight_path = graph_dir / "graph_weight.npz"
    sign_path = graph_dir / "graph_sign.npz"
    index_path = graph_dir / "graph_gene_index.json"
    metadata_path = graph_dir / "graph_metadata.json"
    for path in [weight_path, sign_path, index_path, metadata_path]:
        if not path.exists():
            fail(f"Graph file {path} does not exist (run stage3_build_graph.py first?).")

    weight = sp.load_npz(weight_path).tocsr()
    sign = sp.load_npz(sign_path).tocsr()
    gene_index = json.loads(index_path.read_text())
    metadata = json.loads(metadata_path.read_text())

    n_genes = len(gene_index)
    if weight.shape != (n_genes, n_genes) or sign.shape != (n_genes, n_genes):
        fail(
            f"Graph files in {graph_dir} are misaligned: weight.shape={weight.shape}, "
            f"sign.shape={sign.shape}, but the gene index has {n_genes} entries."
        )
    return weight, sign, gene_index, metadata
