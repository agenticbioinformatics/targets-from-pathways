"""Tests that schemas.py rejects malformed data as loudly as it accepts valid data.

Per Stage 0 of PROMPTS.md: each schema must reject a deliberately malformed
frame (wrong dtype, missing column, an out-of-vocabulary value, and a
gene_id that looks like an HGNC symbol rather than an Ensembl ID) with a
clear error naming the offending column.
"""

from copy import deepcopy

import pandas as pd
import pandera.errors
import pytest
from pydantic import ValidationError

from schemas import (
    GeneSetsSchema,
    GenesSchema,
    InteractionsSchema,
    Manifest,
    OTAssociationsSchema,
    assert_foreign_key,
)

# --------------------------------------------------------------------------
# Valid fixtures
# --------------------------------------------------------------------------


@pytest.fixture
def valid_genes() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_id": ["ENSG00000141510", "ENSG00000073756"],
            "symbol": ["TP53", "PTGS2"],
            "synonyms": [["P53"], ["COX2", "PGHS2"]],
            "biotype": ["protein_coding", "protein_coding"],
        }
    )


@pytest.fixture
def valid_gene_sets() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "set_id": ["R-HSA-1", "R-HSA-1"],
            "set_name": ["Pathway A", "Pathway A"],
            "source_db": ["reactome", "reactome"],
            "source_version": ["90", "90"],
            "gene_id": ["ENSG00000141510", "ENSG00000073756"],
            "hierarchy_level": [1, None],
            "parent_id": [None, "R-HSA-0"],
        }
    )


@pytest.fixture
def valid_interactions() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_a": ["ENSG00000141510"],
            "gene_b": ["ENSG00000073756"],
            "directed": [True],
            "sign": [1],
            "source_db": ["reactome"],
            "source_version": ["90"],
            "evidence_type": ["curated"],
            "confidence": [0.9],
        }
    )


@pytest.fixture
def valid_ot_associations() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_id": ["ENSG00000141510"],
            "disease_id": ["EFO_0005755"],
            "datatype_id": ["genetic_association"],
            "datasource_id": ["gwas_catalog"],
            "score": [0.5],
        }
    )


@pytest.fixture
def valid_manifest_dict() -> dict:
    return {
        "run_id": "run-1",
        "git_commit": "a" * 40,
        "created_at": "2026-07-30T00:00:00Z",
        "seed": 42,
        "resolved_target": {"input": "PTGS2", "gene_id": "ENSG00000073756", "symbol": "PTGS2"},
        "disease": {"efo_id": "EFO_0005755", "name": "rheumatic disease", "n_associated_genes": 100},
        "sources": [
            {
                "db": "reactome",
                "version": "90",
                "files": [{"path": "ReactomePathways.gmt", "sha256": "a" * 64, "bytes": 1024}],
            }
        ],
        "cli_parameters": {"seed": 42, "min_set_size": 10},
        "output_artifacts": [{"path": "genes.parquet", "sha256": "b" * 64}],
        "coverage_report": {
            "reactome": {
                "source": "reactome",
                "raw_gene_count": 100,
                "mapped_count": 90,
                "percent_mapped": 90.0,
                "dropped_count": 10,
                "example_dropped_ids": ["GRB2-1"],
            }
        },
        "scale_report": {
            "gene_set_size_distribution_before_cap": [{"set_size": 10, "count": 5}],
            "gene_set_size_distribution_after_cap": [{"set_size": 10, "count": 5}],
            "sets_retained": 5,
            "interaction_counts_by_sign": {"1": 10, "0": 5, "-1": 3},
            "interaction_counts_by_evidence_type": {"curated": 15, "predicted": 3},
            "projected_comembership_edge_count": 45,
        },
    }


# --------------------------------------------------------------------------
# Valid frames pass
# --------------------------------------------------------------------------


def test_valid_genes_passes(valid_genes):
    GenesSchema.validate(valid_genes)


def test_valid_gene_sets_passes(valid_gene_sets):
    GeneSetsSchema.validate(valid_gene_sets)


def test_valid_interactions_passes(valid_interactions):
    InteractionsSchema.validate(valid_interactions)


def test_valid_ot_associations_passes(valid_ot_associations):
    OTAssociationsSchema.validate(valid_ot_associations)


def test_valid_manifest_passes(valid_manifest_dict):
    Manifest.model_validate(valid_manifest_dict)


# --------------------------------------------------------------------------
# GenesSchema rejects malformed input
# --------------------------------------------------------------------------


def test_genes_rejects_hgnc_symbol_as_gene_id(valid_genes):
    """The exact v1 bug class: a symbol-shaped ID leaking in as gene_id."""
    bad = valid_genes.copy()
    bad.loc[0, "gene_id"] = "GRB2-1"
    with pytest.raises(pandera.errors.SchemaError, match="gene_id"):
        GenesSchema.validate(bad)


def test_genes_rejects_duplicate_gene_id(valid_genes):
    bad = valid_genes.copy()
    bad.loc[1, "gene_id"] = bad.loc[0, "gene_id"]
    with pytest.raises(pandera.errors.SchemaError, match="gene_id"):
        GenesSchema.validate(bad)


def test_genes_rejects_missing_column(valid_genes):
    bad = valid_genes.drop(columns=["biotype"])
    with pytest.raises(pandera.errors.SchemaError, match="biotype"):
        GenesSchema.validate(bad)


def test_genes_accepts_synonyms_as_ndarray(valid_genes):
    """A parquet round-trip turns list[str] cells into numpy.ndarray of str
    (verified against a real genes.parquet), not `list` — GenesSchema must
    accept that container too, or "validate on read" fails for every stage
    after Stage 1 for a column Stage 1 itself writes validly."""
    import numpy as np

    ok = valid_genes.copy()
    ok["synonyms"] = ok["synonyms"].apply(lambda xs: np.array(xs, dtype=object))
    GenesSchema.validate(ok)


def test_genes_rejects_scalar_synonyms(valid_genes):
    """synonyms must be list[str], not a bare string."""
    bad = valid_genes.copy()
    bad["synonyms"] = ["P53", "COX2"]
    with pytest.raises(pandera.errors.SchemaError, match="synonyms"):
        GenesSchema.validate(bad)


def test_genes_rejects_unexpected_extra_column(valid_genes):
    """strict=True must reject a column not declared on the schema."""
    bad = valid_genes.copy()
    bad["debug_note"] = ["leftover", "leftover"]
    with pytest.raises(pandera.errors.SchemaError, match="debug_note"):
        GenesSchema.validate(bad)


def test_genes_rejects_null_symbol(valid_genes):
    bad = valid_genes.copy()
    bad.loc[0, "symbol"] = None
    with pytest.raises(pandera.errors.SchemaError, match="symbol"):
        GenesSchema.validate(bad)


def test_genes_accepts_empty_dataframe(valid_genes):
    """An empty (but correctly typed) frame is valid input, not an error."""
    empty = valid_genes.iloc[0:0]
    assert len(empty) == 0
    GenesSchema.validate(empty)


def test_genes_rejects_wrong_dtype_biotype(valid_genes):
    bad = valid_genes.copy()
    bad["biotype"] = [1, 1]
    with pytest.raises(pandera.errors.SchemaError, match="biotype"):
        GenesSchema.validate(bad)


# --------------------------------------------------------------------------
# GeneSetsSchema rejects malformed input
# --------------------------------------------------------------------------


def test_gene_sets_rejects_unsupported_source_db(valid_gene_sets):
    bad = valid_gene_sets.copy()
    bad.loc[0, "source_db"] = "kegg"
    with pytest.raises(pandera.errors.SchemaError, match="source_db"):
        GeneSetsSchema.validate(bad)


def test_gene_sets_rejects_duplicate_set_gene_pair(valid_gene_sets):
    bad = valid_gene_sets.copy()
    bad.loc[1, "gene_id"] = bad.loc[0, "gene_id"]  # same (set_id, gene_id) twice
    with pytest.raises(pandera.errors.SchemaError):
        GeneSetsSchema.validate(bad)


def test_gene_sets_rejects_non_ensembl_gene_id(valid_gene_sets):
    bad = valid_gene_sets.copy()
    bad.loc[0, "gene_id"] = "TP53"
    with pytest.raises(pandera.errors.SchemaError, match="gene_id"):
        GeneSetsSchema.validate(bad)


def test_gene_sets_rejects_unexpected_extra_column(valid_gene_sets):
    bad = valid_gene_sets.copy()
    bad["debug_note"] = ["leftover", "leftover"]
    with pytest.raises(pandera.errors.SchemaError, match="debug_note"):
        GeneSetsSchema.validate(bad)


def test_gene_sets_hierarchy_level_coerces_parquet_float_roundtrip(valid_gene_sets):
    """The exact scenario coerce=True exists for: a NaN-bearing int column
    that has round-tripped through parquet and now arrives as float64."""
    ok = valid_gene_sets.copy()
    ok["hierarchy_level"] = pd.Series([1.0, None], dtype="float64")
    validated = GeneSetsSchema.validate(ok)
    assert str(validated["hierarchy_level"].dtype) == "Int64"
    assert validated["hierarchy_level"].isna().tolist() == [False, True]
    assert validated["hierarchy_level"].iloc[0] == 1


def test_gene_sets_rejects_missing_column(valid_gene_sets):
    bad = valid_gene_sets.drop(columns=["source_version"])
    with pytest.raises(pandera.errors.SchemaError, match="source_version"):
        GeneSetsSchema.validate(bad)


def test_gene_sets_rejects_wrong_dtype_set_name(valid_gene_sets):
    bad = valid_gene_sets.copy()
    bad["set_name"] = [1, 2]
    with pytest.raises(pandera.errors.SchemaError, match="set_name"):
        GeneSetsSchema.validate(bad)


def test_gene_sets_rejects_non_numeric_hierarchy_level(valid_gene_sets):
    """coerce=True must still raise loudly on genuinely bad data, not just
    accept anything (the exact failure mode that made coerce=True risky at
    the whole-schema level for InteractionsSchema.directed)."""
    bad = valid_gene_sets.copy()
    bad["hierarchy_level"] = ["not-a-number", None]
    with pytest.raises(pandera.errors.SchemaError, match="hierarchy_level"):
        GeneSetsSchema.validate(bad)


# --------------------------------------------------------------------------
# InteractionsSchema rejects malformed input
# --------------------------------------------------------------------------


def test_interactions_rejects_out_of_vocabulary_sign(valid_interactions):
    bad = valid_interactions.copy()
    bad.loc[0, "sign"] = 2
    with pytest.raises(pandera.errors.SchemaError, match="sign"):
        InteractionsSchema.validate(bad)


def test_interactions_rejects_wrong_dtype_directed(valid_interactions):
    bad = valid_interactions.copy()
    bad["directed"] = ["yes"]
    with pytest.raises(pandera.errors.SchemaError, match="directed"):
        InteractionsSchema.validate(bad)


def test_interactions_rejects_unknown_evidence_type(valid_interactions):
    bad = valid_interactions.copy()
    bad.loc[0, "evidence_type"] = "inferred"
    with pytest.raises(pandera.errors.SchemaError, match="evidence_type"):
        InteractionsSchema.validate(bad)


def test_interactions_rejects_confidence_out_of_range(valid_interactions):
    bad = valid_interactions.copy()
    bad.loc[0, "confidence"] = 1.5
    with pytest.raises(pandera.errors.SchemaError, match="confidence"):
        InteractionsSchema.validate(bad)


def test_interactions_rejects_non_ensembl_gene_a(valid_interactions):
    """gene_a/gene_b carry the same FK-format check as GenesSchema.gene_id —
    this is the same v1 "GRB2-1" bug class, and it was untested here."""
    bad = valid_interactions.copy()
    bad.loc[0, "gene_a"] = "GRB2-1"
    with pytest.raises(pandera.errors.SchemaError, match="gene_a"):
        InteractionsSchema.validate(bad)


def test_interactions_rejects_missing_column(valid_interactions):
    bad = valid_interactions.drop(columns=["evidence_type"])
    with pytest.raises(pandera.errors.SchemaError, match="evidence_type"):
        InteractionsSchema.validate(bad)


def test_interactions_rejects_unsupported_source_db(valid_interactions):
    bad = valid_interactions.copy()
    bad.loc[0, "source_db"] = "kegg"
    with pytest.raises(pandera.errors.SchemaError, match="source_db"):
        InteractionsSchema.validate(bad)


def test_interactions_rejects_unexpected_extra_column(valid_interactions):
    bad = valid_interactions.copy()
    bad["debug_note"] = ["leftover"]
    with pytest.raises(pandera.errors.SchemaError, match="debug_note"):
        InteractionsSchema.validate(bad)


def test_interactions_confidence_coerces_int_roundtrip(valid_interactions):
    """A confidence column that is all 0/1 can arrive as int64 (no NaNs to
    force float64 during construction); coerce=True must accept it."""
    ok = valid_interactions.copy()
    ok["confidence"] = pd.Series([1], dtype="int64")
    validated = InteractionsSchema.validate(ok)
    assert validated["confidence"].dtype == "float64"
    assert validated["confidence"].iloc[0] == 1.0


def test_interactions_confidence_bool_column_coerces_silently():
    """Documents a known, accepted residual gap: coerce=True on a numeric
    column will silently cast an (incorrect) bool-dtype column to 0.0/1.0,
    which then satisfies in_range(0, 1) without raising. This is the same
    class of bug as the "directed" schema-wide-coerce trap, but deliberately
    left in place here (narrow blast radius, see the code comment on
    InteractionsSchema.confidence) rather than fixed, since numeric-column
    coercion is required for genuine int64/float64 parquet round-tripping."""
    df = pd.DataFrame(
        {
            "gene_a": ["ENSG00000141510"],
            "gene_b": ["ENSG00000073756"],
            "directed": [True],
            "sign": [1],
            "source_db": ["reactome"],
            "source_version": ["90"],
            "evidence_type": ["curated"],
            "confidence": pd.Series([True], dtype="bool"),
        }
    )
    validated = InteractionsSchema.validate(df)
    assert validated["confidence"].iloc[0] == 1.0


# --------------------------------------------------------------------------
# OTAssociationsSchema rejects malformed input
# --------------------------------------------------------------------------


def test_ot_associations_rejects_score_out_of_range(valid_ot_associations):
    bad = valid_ot_associations.copy()
    bad.loc[0, "score"] = 1.2
    with pytest.raises(pandera.errors.SchemaError, match="score"):
        OTAssociationsSchema.validate(bad)


def test_ot_associations_rejects_non_ensembl_gene_id(valid_ot_associations):
    """Same FK-format check, same v1 "GRB2-1" bug class, untested here."""
    bad = valid_ot_associations.copy()
    bad.loc[0, "gene_id"] = "GRB2-1"
    with pytest.raises(pandera.errors.SchemaError, match="gene_id"):
        OTAssociationsSchema.validate(bad)


def test_ot_associations_rejects_unexpected_extra_column(valid_ot_associations):
    bad = valid_ot_associations.copy()
    bad["debug_note"] = ["leftover"]
    with pytest.raises(pandera.errors.SchemaError, match="debug_note"):
        OTAssociationsSchema.validate(bad)


def test_ot_associations_rejects_missing_column(valid_ot_associations):
    bad = valid_ot_associations.drop(columns=["datasource_id"])
    with pytest.raises(pandera.errors.SchemaError, match="datasource_id"):
        OTAssociationsSchema.validate(bad)


def test_ot_associations_rejects_wrong_dtype_datatype_id(valid_ot_associations):
    bad = valid_ot_associations.copy()
    bad["datatype_id"] = [1]
    with pytest.raises(pandera.errors.SchemaError, match="datatype_id"):
        OTAssociationsSchema.validate(bad)


# --------------------------------------------------------------------------
# assert_foreign_key
# --------------------------------------------------------------------------


def test_assert_foreign_key_passes_when_all_ids_known(valid_genes, valid_gene_sets):
    # valid_gene_sets.gene_id is drawn entirely from valid_genes.gene_id.
    assert_foreign_key(valid_gene_sets, "gene_id", valid_genes)


def test_assert_foreign_key_raises_on_missing_gene_id(valid_genes, valid_gene_sets):
    bad = valid_gene_sets.copy()
    bad.loc[0, "gene_id"] = "ENSG00000000000"  # not in valid_genes
    with pytest.raises(ValueError, match="ENSG00000000000"):
        assert_foreign_key(bad, "gene_id", valid_genes)


def test_assert_foreign_key_checks_both_interaction_endpoints(valid_genes, valid_interactions):
    bad = valid_interactions.copy()
    bad.loc[0, "gene_b"] = "ENSG00000000000"  # not in valid_genes
    with pytest.raises(ValueError, match="gene_b"):
        assert_foreign_key(bad, "gene_b", valid_genes)


# --------------------------------------------------------------------------
# Manifest rejects malformed input
# --------------------------------------------------------------------------


def test_manifest_rejects_missing_field(valid_manifest_dict):
    bad = deepcopy(valid_manifest_dict)
    del bad["seed"]
    with pytest.raises(ValidationError, match="seed"):
        Manifest.model_validate(bad)


def test_manifest_rejects_non_ensembl_target_gene_id(valid_manifest_dict):
    bad = deepcopy(valid_manifest_dict)
    bad["resolved_target"]["gene_id"] = "PTGS2"
    with pytest.raises(ValidationError, match="gene_id"):
        Manifest.model_validate(bad)


def test_manifest_rejects_short_sha256(valid_manifest_dict):
    bad = deepcopy(valid_manifest_dict)
    bad["output_artifacts"][0]["sha256"] = "deadbeef"
    with pytest.raises(ValidationError, match="sha256"):
        Manifest.model_validate(bad)


def test_manifest_rejects_unsupported_source_db(valid_manifest_dict):
    bad = deepcopy(valid_manifest_dict)
    bad["sources"][0]["db"] = "kegg"
    with pytest.raises(ValidationError, match="db"):
        Manifest.model_validate(bad)


def test_manifest_accepts_opentargets_source_db(valid_manifest_dict):
    """opentargets is a valid Source.db (acquisition provenance) even though
    it is not in SUPPORTED_SOURCE_DBS (pathway-curation source_db values)."""
    ok = deepcopy(valid_manifest_dict)
    ok["sources"].append(
        {
            "db": "opentargets",
            "version": "26.06",
            "files": [{"path": "target/part-00000.parquet", "sha256": "c" * 64, "bytes": 2048}],
        }
    )
    Manifest.model_validate(ok)


def test_manifest_rejects_mapped_count_exceeding_raw(valid_manifest_dict):
    bad = deepcopy(valid_manifest_dict)
    bad["coverage_report"]["reactome"]["mapped_count"] = 1000
    with pytest.raises(ValidationError, match="mapped_count"):
        Manifest.model_validate(bad)


def test_manifest_rejects_percent_mapped_inconsistent_with_counts(valid_manifest_dict):
    """percent_mapped gates the pipeline's 90%-coverage-floor check (see
    README.md); it must not be able to silently drift from the counts it is
    supposed to summarize (e.g. a stale value left over from a different
    computation)."""
    bad = deepcopy(valid_manifest_dict)
    # raw=100, mapped=90 implies percent_mapped == 90.0, not 50.0.
    bad["coverage_report"]["reactome"]["percent_mapped"] = 50.0
    with pytest.raises(ValidationError, match="percent_mapped"):
        Manifest.model_validate(bad)


def test_manifest_accepts_percent_mapped_within_rounding_tolerance(valid_manifest_dict):
    """Reported percent_mapped is allowed to be a rounded/truncated value,
    not necessarily bit-exact to the raw division."""
    ok = deepcopy(valid_manifest_dict)
    ok["coverage_report"]["reactome"]["raw_gene_count"] = 93
    ok["coverage_report"]["reactome"]["mapped_count"] = 87
    ok["coverage_report"]["reactome"]["percent_mapped"] = 93.5  # true value 93.548...
    ok["coverage_report"]["reactome"]["dropped_count"] = 6  # 93 - 87
    Manifest.model_validate(ok)


def test_manifest_rejects_dropped_count_inconsistent_with_counts(valid_manifest_dict):
    bad = deepcopy(valid_manifest_dict)
    # raw=100, mapped=90 implies dropped_count == 10, not 3.
    bad["coverage_report"]["reactome"]["dropped_count"] = 3
    with pytest.raises(ValidationError, match="dropped_count"):
        Manifest.model_validate(bad)
