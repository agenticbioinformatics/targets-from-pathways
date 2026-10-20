"""Stage 1 - ingest and canonicalize.

Downloads (or reads from --data-dir) the version-pinned Open Targets and Reactome files, maps every
gene identifier to an Ensembl gene ID and writes genes / gene_sets / interactions / ot_disease_subset
parquet files, coverage_report.json, scale_report.json and manifest.json into --out-dir.
This is the only stage that maps gene identifiers; later stages read only these outputs.
"""

import argparse
import hashlib
import json
import logging
import re
import subprocess
import urllib.error
import urllib.request
import zipfile
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pyarrow.compute as pa_compute
import pyarrow.dataset as pa_dataset
import pyarrow.parquet as pa_parquet

from stage0_schemas import (
    ENSEMBL_GENE_ID_PATTERN,
    SUPPORTED_SOURCE_DBS,
    GeneSetSizeBucket,
    GeneSetsSchema,
    GenesSchema,
    InteractionsSchema,
    Manifest,
    OTAssociationsSchema,
    OutputArtifact,
    ResolvedDisease,
    ScaleReport,
    Source,
    SourceFile,
    assert_foreign_key,
    fail,
)

logger = logging.getLogger("ingest")

OPENTARGETS_VERSION = "26.06"
REACTOME_VERSION = "97"
# The functional-interaction (FI) file is versioned by date, separately from Reactome releases.
REACTOME_FI_VERSION = "04142025"

OPENTARGETS_URL = f"https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/{OPENTARGETS_VERSION}/output"
REACTOME_URL = f"https://download.reactome.org/{REACTOME_VERSION}"

# Single files, pinned by sha256 and size.
PINNED_FILES = {
    "reactome_gmt": {
        "url": f"{REACTOME_URL}/ReactomePathways.gmt.zip",
        "sha256": "8c1dbc8578431da5d2d5118262718c60b553a9be3398e93658daa069e4a9afd4",
        "size": 298479,
        "path": "reactome/ReactomePathways.gmt.zip",
    },
    # Gene-level file. The "_PE_All_Levels" variant contains physical-entity labels and must not be used.
    "reactome_ensembl2reactome": {
        "url": f"{REACTOME_URL}/Ensembl2Reactome.txt",
        "sha256": "34852dc2ac258a851bad914708c0045a636161dbbeca604f04343174cca8d7fb",
        "size": 183277348,
        "path": "reactome/Ensembl2Reactome.txt",
    },
    "reactome_fi": {
        "url": f"https://reactome.org/download/tools/ReactomeFIs/FIsInGene_{REACTOME_FI_VERSION}_with_annotations.txt.zip",
        "sha256": "e5c373e3baff9764c07edb36c1adc8b714a3e1b44dfe6f8892bafaeb05fb2335",
        "size": 1182418,
        "path": f"reactome/FIsInGene_{REACTOME_FI_VERSION}_with_annotations.txt.zip",
    },
    # Needed only for the disease name and for EFO -> MONDO cross-references.
    "opentargets_disease": {
        "url": f"{OPENTARGETS_URL}/disease/disease.parquet",
        "sha256": "b328c5e775bc40d954f2e058aa59e8e972f2883a262d096f8dffe365736a1aed",
        "size": 7150787,
        "path": f"opentargets/{OPENTARGETS_VERSION}/disease/disease.parquet",
    },
}

# Spark output directories: part-file names contain random UUIDs, so they cannot be pinned in advance.
# Instead the checksums seen on first download are stored in a local .checksums.json and enforced afterwards.
PINNED_DIRECTORIES = {
    "opentargets_target": {
        "url": f"{OPENTARGETS_URL}/target/",
        "path": f"opentargets/{OPENTARGETS_VERSION}/target",
    },
    "opentargets_association_by_datasource_direct": {
        "url": f"{OPENTARGETS_URL}/association_by_datasource_direct/",
        "path": f"opentargets/{OPENTARGETS_VERSION}/association_by_datasource_direct",
    },
}

# The association table only has a datasource column; this maps it to the Open Targets datatype
# (https://platform-docs.opentargets.org/evidence). An unknown datasource stops the run, because
# Stages 2 and 4 must be able to exclude pathway-derived datatypes reliably.
DATASOURCE_TO_DATATYPE = {
    "gwas_credible_sets": "genetic_association",
    "gene_burden": "genetic_association",
    "eva": "genetic_association",
    "genomics_england": "genetic_association",
    "gene2phenotype": "genetic_association",
    "uniprot_literature": "genetic_association",
    "uniprot_variants": "genetic_association",
    "orphanet": "genetic_association",
    "clingen": "genetic_association",
    "cancer_gene_census": "somatic_mutation",
    "intogen": "somatic_mutation",
    "eva_somatic": "somatic_mutation",
    "cancer_biomarkers": "affected_pathway",
    "crispr_screen": "affected_pathway",
    "crispr": "affected_pathway",
    "reactome": "affected_pathway",
    "europepmc": "literature",
    "expression_atlas": "rna_expression",
    "impc": "animal_model",
    "clinical_precedence": "known_drug",
}

# Placeholder thresholds, not derived from data.
MIN_COVERAGE_PERCENT = 90.0
MIN_GENETIC_ASSOCIATION_GENES = 5

# FI "Direction" codes: an optional left cap (< or |) and right cap (> or |) around a dash.
# ">" / "<" = activation (+1), "|" = inhibition (-1). Each code maps to the edges it implies,
# as (forward, sign) where forward means Gene1 -> Gene2. The bare "-" is an undirected edge.
FI_DIRECTION_CODES = {
    "->": [(True, 1)],
    "-|": [(True, -1)],
    "<-": [(False, 1)],
    "|-": [(False, -1)],
    "<->": [(True, 1), (False, 1)],
    "|-|": [(True, -1), (False, -1)],
    "<-|": [(False, 1), (True, -1)],
    "|->": [(False, -1), (True, 1)],
}

# reactome.org blocks urllib's default User-Agent (HTTP 403).
USER_AGENT = "Mozilla/5.0 (X11; Linux x86_64) targets-from-pathways/stage1_ingest.py"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="stage1_ingest.py",
        description="Stage 1 - ingest & canonicalize. Downloads (or reads from --data-dir) the "
        "version-pinned Open Targets / Reactome sources, maps every gene id to Ensembl, "
        "subsets OT associations to --disease, and writes genes/gene_sets/interactions/"
        "ot_disease_subset parquet + coverage_report.json + scale_report.json + manifest.json "
        "into --out-dir.",
    )
    parser.add_argument("--disease", required=True, help="Disease: EFO / MONDO / other OT-recognized ontology ID.")
    parser.add_argument("--data-dir", required=True, type=Path, help="Source cache directory (downloads land here).")
    parser.add_argument("--out-dir", required=True, type=Path, help="Directory to write this run's Stage 1 files.")
    parser.add_argument(
        "--pathway-db", choices=sorted(SUPPORTED_SOURCE_DBS), default="reactome",
        help="Pathway database (only 'reactome' is supported).",
    )
    parser.add_argument("--no-download", action="store_true", help="Never use the network; --data-dir must be populated.")
    parser.add_argument("--seed", type=int, default=0, help="Random seed, recorded in the manifest (default: 0).")
    parser.add_argument("--min-set-size", type=int, default=10, help="Drop gene sets smaller than this (default: 10).")
    parser.add_argument("--max-set-size", type=int, default=200, help="Drop gene sets larger than this (default: 200).")
    parser.add_argument(
        "--fi-curated-only", action=argparse.BooleanOptionalAction, default=True,
        help="Keep only curated Reactome FI rows (predicted rows carry no direction/sign).",
    )
    parser.add_argument("--fi-min-score", type=float, default=None, help="Drop FI rows with Score below this.")
    parser.add_argument(
        "--allow-low-coverage", action="store_true",
        help="Do not fail on <90%% identifier-mapping coverage or too few genetic_association genes.",
    )
    return parser.parse_args(argv)


def sha256_of_file(path):
    hasher = hashlib.sha256()
    with open(path, "rb") as file:
        while True:
            chunk = file.read(1024 * 1024)
            if not chunk:
                break
            hasher.update(chunk)
    return hasher.hexdigest()


def open_url(url, timeout):
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    return urllib.request.urlopen(request, timeout=timeout)


def download_file(url, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Download to a temporary name so an interrupted download never looks like a complete file.
    partial = destination.with_suffix(destination.suffix + ".part")
    try:
        with open_url(url, timeout=120) as response, open(partial, "wb") as file:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                file.write(chunk)
    except (urllib.error.URLError, TimeoutError) as error:
        partial.unlink(missing_ok=True)
        fail(f"Failed to download {url}: {error}")
    partial.replace(destination)


def list_remote_parquet_files(url):
    """Return the *.snappy.parquet file names linked from an Apache directory-index page."""
    try:
        with open_url(url, timeout=60) as response:
            page = response.read().decode("utf-8", errors="replace")
    except (urllib.error.URLError, TimeoutError) as error:
        fail(f"Failed to list {url}: {error}")
    return re.findall(r'href="([^"?][^"]*\.snappy\.parquet)"', page)


def file_matches(path, expected_sha256, expected_size):
    if not path.exists():
        return False
    if path.stat().st_size != expected_size:
        return False
    return sha256_of_file(path) == expected_sha256


def get_pinned_file(key, data_dir, no_download):
    """Return the local path of a pinned file, downloading it if needed. Never accept a checksum mismatch."""
    pinned = PINNED_FILES[key]
    destination = data_dir / pinned["path"]
    if file_matches(destination, pinned["sha256"], pinned["size"]):
        return destination

    if no_download:
        fail(
            f"[{key}] missing or checksum-invalid at {destination}.\n"
            f"  Obtain it manually:\n"
            f"    curl -L -o {destination} '{pinned['url']}'\n"
            f"  Expected sha256: {pinned['sha256']}\n"
            f"  Expected size:   {pinned['size']} bytes"
        )

    logger.info("[%s] downloading %s -> %s", key, pinned["url"], destination)
    destination.unlink(missing_ok=True)
    download_file(pinned["url"], destination)
    if not file_matches(destination, pinned["sha256"], pinned["size"]):
        actual_sha256 = sha256_of_file(destination)
        destination.unlink(missing_ok=True)
        fail(
            f"[{key}] checksum mismatch after download from {pinned['url']}.\n"
            f"  expected sha256={pinned['sha256']} size={pinned['size']}\n"
            f"  actual   sha256={actual_sha256}"
        )
    return destination


def get_directory_files(key, data_dir, no_download):
    """Return the local parquet part files of a pinned directory, downloading them if needed."""
    remote = PINNED_DIRECTORIES[key]
    local_dir = data_dir / remote["path"]
    checksum_file = local_dir / ".checksums.json"

    if checksum_file.exists():
        checksums = json.loads(checksum_file.read_text())
        paths = []
        for name in sorted(checksums):
            path = local_dir / name
            if not file_matches(path, checksums[name]["sha256"], checksums[name]["size"]):
                fail(
                    f"[{key}] {path} no longer matches {checksum_file}. "
                    f"Delete {local_dir} and re-run to download again, or restore the original files."
                )
            paths.append(path)
        return paths

    if no_download:
        local_files = []
        if local_dir.exists():
            local_files = sorted(local_dir.glob("*.parquet"))
        if not local_files:
            fail(
                f"[{key}] no local files under {local_dir} and --no-download is set.\n"
                f"  Copy every *.snappy.parquet file from {remote['url']} into {local_dir}, then re-run."
            )
        # Files copied in by hand: pin what is there, but completeness cannot be checked offline.
        logger.warning("[%s] --no-download: pinning %d local file(s) in %s without checking the remote listing.",
                       key, len(local_files), local_dir)
        file_names = []
        for path in local_files:
            file_names.append(path.name)
    else:
        logger.info("[%s] listing %s", key, remote["url"])
        file_names = sorted(list_remote_parquet_files(remote["url"]))
        if not file_names:
            fail(f"[{key}] directory listing at {remote['url']} returned no *.snappy.parquet files.")
        local_dir.mkdir(parents=True, exist_ok=True)

    checksums = {}
    paths = []
    for name in file_names:
        path = local_dir / name
        if not path.exists():
            logger.info("[%s] downloading %s", key, name)
            download_file(remote["url"] + name, path)
        checksums[name] = {"sha256": sha256_of_file(path), "size": path.stat().st_size}
        paths.append(path)
    checksum_file.write_text(json.dumps(checksums, indent=2, sort_keys=True))
    return paths


def unzip_single_file(zip_path):
    output_path = zip_path.parent / zip_path.stem
    if output_path.exists():
        return output_path
    with zipfile.ZipFile(zip_path) as archive:
        members = archive.namelist()
        if len(members) != 1:
            fail(f"{zip_path} expected to contain exactly one file, found {members}")
        with archive.open(members[0]) as source, open(output_path, "wb") as destination:
            destination.write(source.read())
    return output_path


@dataclass
class SymbolIndex:
    symbol_to_gene: dict
    synonym_to_genes: dict


def build_symbol_index(genes_df):
    symbol_to_gene = {}
    ambiguous_symbols = set()
    for gene_id, symbol in zip(genes_df["gene_id"], genes_df["symbol"]):
        key = symbol.upper()
        if key in symbol_to_gene and symbol_to_gene[key] != gene_id:
            ambiguous_symbols.add(key)
        else:
            symbol_to_gene[key] = gene_id
    # A symbol shared by two genes cannot be mapped safely.
    for key in ambiguous_symbols:
        del symbol_to_gene[key]

    synonym_to_genes = {}
    for gene_id, synonyms in zip(genes_df["gene_id"], genes_df["synonyms"]):
        for synonym in synonyms:
            key = synonym.upper()
            if key not in synonym_to_genes:
                synonym_to_genes[key] = set()
            synonym_to_genes[key].add(gene_id)
    return SymbolIndex(symbol_to_gene, synonym_to_genes)


def resolve_symbol(raw_symbol, index):
    """Map a symbol to a gene ID: approved symbol first, then a synonym. Returns (gene_id or None, reason)."""
    key = raw_symbol.strip().upper()
    if key in index.symbol_to_gene:
        return index.symbol_to_gene[key], "approved_symbol"
    genes = index.synonym_to_genes.get(key)
    if not genes:
        return None, "unmapped"
    if len(genes) > 1:
        return None, "ambiguous_synonym"
    return list(genes)[0], "synonym"


def read_target_table(target_paths):
    columns = ["id", "approvedSymbol", "biotype", "symbolSynonyms", "obsoleteSymbols", "proteinIds", "transcriptIds"]
    return pa_dataset.dataset(sorted(target_paths), format="parquet").to_table(columns=columns).to_pandas()


def struct_labels(structs):
    # Missing lists arrive as None; a non-empty list arrives as a numpy array, so "if structs" cannot be used.
    labels = []
    if structs is not None:
        for struct in structs:
            labels.append(struct["label"])
    return labels


def build_genes_table(target_df):
    """The Open Targets target table is the authority for gene IDs, symbols and synonyms."""
    all_synonyms = []
    for current_synonyms, obsolete_symbols in zip(target_df["symbolSynonyms"], target_df["obsoleteSymbols"]):
        synonyms = set(struct_labels(current_synonyms)) | set(struct_labels(obsolete_symbols))
        all_synonyms.append(sorted(synonyms))

    genes_df = pd.DataFrame({
        "gene_id": target_df["id"],
        "symbol": target_df["approvedSymbol"],
        "synonyms": all_synonyms,
        "biotype": target_df["biotype"].fillna(""),
    })
    return genes_df.drop_duplicates(subset="gene_id").reset_index(drop=True)


def build_protein_and_transcript_to_gene(target_df):
    """Map Ensembl protein (ENSP) and transcript (ENST) IDs to their gene."""
    id_to_gene = {}
    for gene_id, proteins, transcripts in zip(target_df["id"], target_df["proteinIds"], target_df["transcriptIds"]):
        if proteins is not None:
            for protein in proteins:
                if protein["source"] == "ensembl_PRO":
                    id_to_gene[protein["id"]] = gene_id
        if transcripts is not None:
            for transcript_id in transcripts:
                id_to_gene[transcript_id] = gene_id
    return id_to_gene


@dataclass
class CoverageTracker:
    source: str
    raw: set = field(default_factory=set)
    mapped: set = field(default_factory=set)
    dropped: set = field(default_factory=set)

    def record(self, raw_id, gene_id):
        self.raw.add(raw_id)
        if gene_id is None:
            self.dropped.add(raw_id)
        else:
            self.mapped.add(raw_id)

    def to_entry(self):
        n_raw = len(self.raw)
        n_mapped = len(self.mapped)
        percent_mapped = 0.0
        if n_raw > 0:
            percent_mapped = round(100.0 * n_mapped / n_raw, 2)
        return {
            "source": self.source,
            "raw_gene_count": n_raw,
            "mapped_count": n_mapped,
            "percent_mapped": percent_mapped,
            "dropped_count": n_raw - n_mapped,
            "example_dropped_ids": sorted(self.dropped)[:20],
        }


def parse_gmt(gmt_path):
    """Return [(set_id, set_name, [gene symbols])]. GMT columns: pathway name, Reactome ID, genes..."""
    gene_sets = []
    for line in gmt_path.read_text().splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        symbols = []
        for symbol in fields[2:]:
            if symbol:
                symbols.append(symbol)
        gene_sets.append((fields[1], fields[0], symbols))
    return gene_sets


def build_gene_sets_table(gmt_path, index, source_db, coverage):
    rows = []
    for set_id, set_name, symbols in parse_gmt(gmt_path):
        genes_in_set = set()
        for symbol in symbols:
            gene_id, _ = resolve_symbol(symbol, index)
            coverage.record(symbol, gene_id)
            if gene_id is None or gene_id in genes_in_set:
                continue
            genes_in_set.add(gene_id)
            rows.append([set_id, set_name, source_db, REACTOME_VERSION, gene_id, None, None])
    columns = ["set_id", "set_name", "source_db", "source_version", "gene_id", "hierarchy_level", "parent_id"]
    return pd.DataFrame(rows, columns=columns)


def build_interactions_table(fi_path, index, source_db, coverage, curated_only, min_score):
    """Turn the Reactome FI file into directed, signed gene-gene edges. Swap the interaction source here."""
    fi_df = pd.read_csv(fi_path, sep="\t", dtype=str)

    # pandas pads a short row with NaN instead of raising; a corrupt source row must stop the run.
    required_columns = ["Gene1", "Gene2", "Annotation", "Direction", "Score"]
    malformed = fi_df[fi_df[required_columns].isnull().any(axis=1)]
    if not malformed.empty:
        fail(
            f"{fi_path} has {len(malformed)} malformed row(s) with a missing field "
            f"(e.g. line {malformed.index[0] + 2}: {malformed.iloc[0].to_dict()})."
        )
    fi_df["Score"] = fi_df["Score"].astype(float)

    rows = []
    for gene1, gene2, annotation, direction, score in zip(
        fi_df["Gene1"], fi_df["Gene2"], fi_df["Annotation"], fi_df["Direction"], fi_df["Score"]
    ):
        # Predicted rows have Annotation == "predicted"; everything else is a list of curated relations.
        if annotation == "predicted":
            evidence_type = "predicted"
        else:
            evidence_type = "curated"
        if curated_only and evidence_type != "curated":
            continue
        if min_score is not None and score < min_score:
            continue

        gene_a, _ = resolve_symbol(gene1, index)
        coverage.record(gene1, gene_a)
        gene_b, _ = resolve_symbol(gene2, index)
        coverage.record(gene2, gene_b)
        if gene_a is None or gene_b is None or gene_a == gene_b:
            continue

        if direction == "-":
            rows.append([gene_a, gene_b, False, 0, source_db, REACTOME_FI_VERSION, evidence_type, score])
            continue
        if direction not in FI_DIRECTION_CODES:
            logger.warning("Unrecognized FI Direction code %r for %s/%s - dropping row.", direction, gene1, gene2)
            continue
        for forward, sign in FI_DIRECTION_CODES[direction]:
            if forward:
                rows.append([gene_a, gene_b, True, sign, source_db, REACTOME_FI_VERSION, evidence_type, score])
            else:
                rows.append([gene_b, gene_a, True, sign, source_db, REACTOME_FI_VERSION, evidence_type, score])

    columns = ["gene_a", "gene_b", "directed", "sign", "source_db", "source_version", "evidence_type", "confidence"]
    return pd.DataFrame(rows, columns=columns).drop_duplicates().reset_index(drop=True)


def ensembl2reactome_coverage(e2r_path, genes_df, protein_and_transcript_to_gene):
    """Coverage cross-check only: this file lists pathway members by ENSG/ENSP/ENST ID, not by the symbols
    used in the GMT and FI files, so its rows cannot be merged into gene_sets reliably."""
    known_genes = set(genes_df["gene_id"])
    tracker = CoverageTracker(source="reactome_ensembl2reactome")
    with open(e2r_path) as file:
        for line in file:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 6 or fields[5] != "Homo sapiens":
                continue
            raw_id = fields[0]
            unversioned_id = re.sub(r"\.\d+$", "", raw_id)
            if unversioned_id in known_genes:
                gene_id = unversioned_id
            else:
                gene_id = protein_and_transcript_to_gene.get(unversioned_id)
            tracker.record(raw_id, gene_id)
    return tracker.to_entry()


def load_disease_table(disease_parquet):
    return pa_parquet.read_table(disease_parquet, columns=["id", "name", "dbXRefs"]).to_pandas()


def resolve_disease_id(raw_disease, disease_df):
    """Return (disease_id, name) or (None, None)."""
    disease_id = raw_disease.strip().replace(":", "_")
    match = disease_df.loc[disease_df["id"] == disease_id]
    if not match.empty:
        return match.iloc[0]["id"], match.iloc[0]["name"]

    # Open Targets 26.06 uses MONDO as primary IDs; older EFO IDs survive only in dbXRefs as "EFO:0000685".
    xref = re.sub(r"^([A-Za-z]+)_", r"\1:", disease_id)
    for candidate_id, name, xrefs in zip(disease_df["id"], disease_df["name"], disease_df["dbXRefs"]):
        if xrefs is not None and xref in xrefs:
            return candidate_id, name
    return None, None


def load_disease_association_subset(association_paths, disease_id):
    dataset = pa_dataset.dataset(sorted(association_paths), format="parquet")
    associations = dataset.to_table(
        columns=["targetId", "diseaseId", "aggregationType", "aggregationValue", "associationScore"],
        filter=pa_compute.field("diseaseId") == disease_id,
    ).to_pandas()

    unexpected_types = set(associations["aggregationType"]) - {"datasourceId"}
    if unexpected_types:
        fail(f"aggregationType has unexpected value(s) {sorted(unexpected_types)} (expected only 'datasourceId').")
    unknown_datasources = set(associations["aggregationValue"]) - set(DATASOURCE_TO_DATATYPE)
    if unknown_datasources:
        fail(f"Unrecognized OT datasource id(s) {sorted(unknown_datasources)}: add them to DATASOURCE_TO_DATATYPE.")

    return pd.DataFrame({
        "gene_id": associations["targetId"],
        "disease_id": associations["diseaseId"],
        "datatype_id": associations["aggregationValue"].map(DATASOURCE_TO_DATATYPE),
        "datasource_id": associations["aggregationValue"],
        "score": associations["associationScore"].astype(float),
    })


def resolve_gene(raw_gene, genes_df, index):
    """Return (gene_id, symbol) for a symbol or Ensembl ID, or (None, None)."""
    raw_gene = raw_gene.strip()
    if ENSEMBL_GENE_ID_PATTERN.match(raw_gene):
        gene_id = raw_gene
    else:
        gene_id, _ = resolve_symbol(raw_gene, index)
    match = genes_df.loc[genes_df["gene_id"] == gene_id]
    if gene_id is None or match.empty:
        return None, None
    return gene_id, match.iloc[0]["symbol"]


def size_distribution(sizes):
    counts = {}
    for size in sizes:
        counts[size] = counts.get(size, 0) + 1
    buckets = []
    for size in sorted(counts):
        buckets.append(GeneSetSizeBucket(set_size=size, count=counts[size]))
    return buckets


def count_values(series):
    counts = {}
    for value, count in series.value_counts().items():
        counts[str(value)] = int(count)
    return counts


def apply_size_cap(gene_sets_df, interactions_df, min_set_size, max_set_size):
    """Drop gene sets outside [min, max] size and report sizes and the projected Stage 3 edge count."""
    sizes_before = gene_sets_df.groupby("set_id")["gene_id"].nunique()
    kept_set_ids = sizes_before[(sizes_before >= min_set_size) & (sizes_before <= max_set_size)].index
    kept_gene_sets = gene_sets_df[gene_sets_df["set_id"].isin(kept_set_ids)].reset_index(drop=True)
    sizes_after = sizes_before.loc[kept_set_ids].tolist()

    projected_edges = 0
    for size in sizes_after:
        projected_edges += size * (size - 1) // 2

    report = ScaleReport(
        gene_set_size_distribution_before_cap=size_distribution(sizes_before.tolist()),
        gene_set_size_distribution_after_cap=size_distribution(sizes_after),
        sets_retained=len(kept_set_ids),
        interaction_counts_by_sign=count_values(interactions_df["sign"]),
        interaction_counts_by_evidence_type=count_values(interactions_df["evidence_type"]),
        projected_comembership_edge_count=projected_edges,
    )
    return report, kept_gene_sets


def describe_files(paths):
    files = []
    for path in paths:
        files.append(SourceFile(path=str(path), sha256=sha256_of_file(path), bytes=path.stat().st_size))
    return files


def write_parquet(df, schema, path):
    schema.validate(df).to_parquet(path, index=False)
    return OutputArtifact(path=str(path), sha256=sha256_of_file(path))


def write_json(data, path):
    path.write_text(json.dumps(data, indent=2, sort_keys=True))
    return OutputArtifact(path=str(path), sha256=sha256_of_file(path))


def current_git_commit():
    try:
        result = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True, timeout=10)
        return result.stdout.strip()
    except Exception:
        return "unknown"


def make_run_id(disease, seed, cli_parameters):
    # Deterministic: the same inputs and options always give the same run_id.
    payload = json.dumps({"disease": disease, "seed": seed, "cli": cli_parameters}, sort_keys=True)
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def run_ingest(args):
    data_dir = args.data_dir
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    gmt_zip = get_pinned_file("reactome_gmt", data_dir, args.no_download)
    e2r_path = get_pinned_file("reactome_ensembl2reactome", data_dir, args.no_download)
    fi_zip = get_pinned_file("reactome_fi", data_dir, args.no_download)
    disease_parquet = get_pinned_file("opentargets_disease", data_dir, args.no_download)
    target_paths = sorted(get_directory_files("opentargets_target", data_dir, args.no_download))
    association_paths = sorted(
        get_directory_files("opentargets_association_by_datasource_direct", data_dir, args.no_download)
    )
    gmt_path = unzip_single_file(gmt_zip)
    fi_path = unzip_single_file(fi_zip)

    sources = [
        Source(db="reactome", version=REACTOME_VERSION, files=describe_files([gmt_zip, e2r_path])),
        Source(db="reactome", version=REACTOME_FI_VERSION, files=describe_files([fi_zip])),
        Source(db="opentargets", version=OPENTARGETS_VERSION,
               files=describe_files([disease_parquet] + target_paths + association_paths)),
    ]

    target_df = read_target_table(target_paths)
    genes_df = build_genes_table(target_df)
    index = build_symbol_index(genes_df)

    disease_id, disease_name = resolve_disease_id(args.disease, load_disease_table(disease_parquet))
    if disease_id is None:
        fail(f"Could not resolve --disease {args.disease!r} in Open Targets {OPENTARGETS_VERSION} (checked id and dbXRefs).")

    ot_disease_subset = OTAssociationsSchema.validate(load_disease_association_subset(association_paths, disease_id))
    assert_foreign_key(ot_disease_subset, "gene_id", genes_df)
    genetic_rows = ot_disease_subset[ot_disease_subset["datatype_id"] == "genetic_association"]
    n_genetic_genes = genetic_rows["gene_id"].nunique()
    if n_genetic_genes < MIN_GENETIC_ASSOCIATION_GENES and not args.allow_low_coverage:
        fail(
            f"Disease {disease_id} ({disease_name}) has only {n_genetic_genes} genetic_association gene(s), "
            f"below the {MIN_GENETIC_ASSOCIATION_GENES}-gene floor for GSEA. Pass --allow-low-coverage to override."
        )

    gmt_coverage = CoverageTracker(source="reactome_gmt")
    all_gene_sets = build_gene_sets_table(gmt_path, index, args.pathway_db, gmt_coverage)
    fi_coverage = CoverageTracker(source="reactome_fi")
    interactions_df = build_interactions_table(
        fi_path, index, args.pathway_db, fi_coverage, args.fi_curated_only, args.fi_min_score
    )
    coverage_report = {
        "reactome_gmt": gmt_coverage.to_entry(),
        "reactome_fi": fi_coverage.to_entry(),
        "reactome_ensembl2reactome": ensembl2reactome_coverage(
            e2r_path, genes_df, build_protein_and_transcript_to_gene(target_df)
        ),
    }

    low_coverage = []
    for source_name, entry in coverage_report.items():
        if entry["percent_mapped"] < MIN_COVERAGE_PERCENT:
            low_coverage.append(f"{source_name}={entry['percent_mapped']}%")
    if low_coverage and not args.allow_low_coverage:
        fail(f"Source(s) mapped below the {MIN_COVERAGE_PERCENT}% floor: {'; '.join(low_coverage)}. "
             f"Pass --allow-low-coverage to override.")

    scale_report, gene_sets_df = apply_size_cap(all_gene_sets, interactions_df, args.min_set_size, args.max_set_size)

    assert_foreign_key(gene_sets_df, "gene_id", genes_df)
    assert_foreign_key(interactions_df, "gene_a", genes_df)
    assert_foreign_key(interactions_df, "gene_b", genes_df)
    output_artifacts = [
        write_parquet(genes_df, GenesSchema, out_dir / "genes.parquet"),
        write_parquet(gene_sets_df, GeneSetsSchema, out_dir / "gene_sets.parquet"),
        write_parquet(interactions_df, InteractionsSchema, out_dir / "interactions.parquet"),
        write_parquet(ot_disease_subset, OTAssociationsSchema, out_dir / "ot_disease_subset.parquet"),
        write_json(coverage_report, out_dir / "coverage_report.json"),
        write_json(scale_report.model_dump(mode="json"), out_dir / "scale_report.json"),
    ]

    cli_parameters = {}
    for name, value in vars(args).items():
        if isinstance(value, Path):
            value = str(value)
        cli_parameters[name] = value
    cli_parameters["opentargets_version"] = OPENTARGETS_VERSION
    cli_parameters["reactome_version"] = REACTOME_VERSION
    cli_parameters["reactome_fi_version"] = REACTOME_FI_VERSION

    manifest = Manifest(
        run_id=make_run_id(disease_id, args.seed, cli_parameters),
        git_commit=current_git_commit(),
        created_at=datetime.now(timezone.utc),
        seed=args.seed,
        disease=ResolvedDisease(efo_id=disease_id, name=disease_name, n_associated_genes=ot_disease_subset["gene_id"].nunique()),
        sources=sources,
        cli_parameters=cli_parameters,
        output_artifacts=output_artifacts,
        coverage_report=coverage_report,
        scale_report=scale_report,
    )
    manifest_path = out_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest.model_dump(mode="json"), indent=2, sort_keys=True))
    logger.info("Stage 1 complete: %s", manifest_path)
    return manifest


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    run_ingest(parse_args(argv))


if __name__ == "__main__":
    main()
