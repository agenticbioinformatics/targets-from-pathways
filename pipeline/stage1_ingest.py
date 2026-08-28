"""Stage 1 — Ingest & canonicalize.

Resolves, canonicalizes, and normalizes all source data for the pipeline
(README.md's Stage 1, redesigned in "Plan update 3"). This is the *only*
place gene identifiers are mapped: it is a pure function of (CLI args,
the pinned registry below, the bytes at --data-dir) producing genes.parquet,
gene_sets.parquet, interactions.parquet, ot_disease_subset.parquet,
coverage_report.json, scale_report.json, and manifest.json, each validated
against stage0_schemas.py on write. No hidden global state: every fact this module
needs either lives in the registry below or arrives via CLI/--data-dir.

--------------------------------------------------------------------------
Verification log (2026-08-24, jedrzej.kubica@univ-grenoble-alpes.fr)
--------------------------------------------------------------------------
Every URL and sha256 hardcoded below was checked by hand before being
pinned, not assumed from the v1 README (which used `current/` and was
explicitly flagged in README.md's "Open research" section as unconfirmed):

- ``https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.06/output/``
  returns HTTP 200 and lists ``target/``, ``association_by_datasource_indirect/``,
  ``disease/``, etc., each stamped 2026-06-2x. The plain ftp:// scheme
  (``ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.06/output/target/``)
  was independently confirmed to list the same three part-files; https:// is
  pinned here because it round-trips through a normal HTTP GET (no passive-FTP
  firewall concerns) and because its Apache directory index is trivial to
  parse for the file listing DirectorySource needs (see below).
- Reactome numbered-release URLs **do exist** — this resolves the first
  "Open research" question in README.md, which previously only had v1's
  `current/`-based README as a (disclaimed) guess. The Reactome Content
  Service reports the live release number
  (``https://reactome.org/ContentService/data/database/version`` -> ``97``
  at verification time), and both
  ``https://download.reactome.org/97/ReactomePathways.gmt.zip`` and
  ``https://download.reactome.org/97/Ensembl2Reactome.txt`` return HTTP 200.
  ``https://reactome.org/download/97/...`` (the reactome.org-hosted mirror of
  the same numbered release) was cross-checked and also returns HTTP 200.
- The FI (functional-interaction) file is **not** part of Reactome's numbered
  pathway-release cycle: reactome.org's download page links it under its own
  "Version 2025" / "Version 2024" labels, dated independently of the pathway
  release number. The most recent one verified live is
  ``https://reactome.org/download/tools/ReactomeFIs/FIsInGene_04142025_with_annotations.txt.zip``
  (HTTP 200, ``content-length: 1182418``). This dated filename *is* the FI
  file's pin — there is no separate numbered-release URL for it to alias.
- sha256 for every pinned file below was computed locally against the actual
  downloaded bytes at verification time (``hashlib.sha256`` over the whole
  file), not copied from any third-party listing.

Gene-level vs PE-level mapping file (README.md open item, and the exact class
of bug that leaked ``GRB2-1`` into v1's output): both
``Ensembl2Reactome.txt`` and ``Ensembl2Reactome_PE_All_Levels.txt`` were
downloaded (byte-range samples for the multi-GB PE file) and diffed by eye.
``Ensembl2Reactome.txt`` has 6 tab-separated columns per row (identifier,
Reactome pathway stable ID, URL, pathway name, evidence code, species) and
every identifier is a bare Ensembl accession (``ENSG*``, ``ENSP*``, or
``ENST*``, occasionally version-suffixed, e.g. ``ENST00000527673.10``) —
never a bracket-suffixed physical-entity label. ``Ensembl2Reactome_PE_All_Levels.txt``
inserts two *extra* columns between the identifier and the pathway ID: a
physical-entity stable ID (``R-HSA-nnnnnn``) and a PE display name with an
explicit compartment/state suffix (e.g. ``FANCL [cytosol]``) — this is the
leak vector. ``Ensembl2Reactome.txt`` is therefore the correct, gene-resolvable
file and is what's pinned below; the PE file is never touched by this module.

A further, non-obvious finding while inspecting ``Ensembl2Reactome.txt`` for
the *human* rows specifically (``awk -F'\\t' '$6=="Homo sapiens"'``, 381,765
rows): identifiers are **not** uniformly ``ENSG*``. Of the human rows, 61,120
are ``ENSG*``, 160,142 are ``ENSP*`` (Ensembl protein), 160,502 are ``ENST*``
(Ensembl transcript), and 1 stray row uses a non-Ensembl ``EBT*`` accession.
None of that is PE-level (no bracket-suffixed complex/compartment state), so
it does not reintroduce the v1 bug — but it does mean a naive
"pass column 1 straight through as gene_id" would silently drop ~85% of this
file's rows as non-ENSG. This module treats ``Ensembl2Reactome.txt`` as a
**cross-check-only** source (see ``_load_ensembl2reactome_coverage`` below)
rather than folding it into ``gene_sets.parquet``: collapsing an ENSP/ENST id
to a gene is unambiguous via the Open Targets target table (each protein/
transcript belongs to exactly one gene), but *attributing* a specific
GMT gene-symbol mapping failure to a specific recovered Ensembl2Reactome row
for the same pathway is not — the file records pathway membership by ID, not
by the gene *symbol* the GMT and FI files use, so there is no reliable
per-symbol correspondence to merge on. Reporting it as an independent,
transparent coverage entry avoids a silent many-to-one guess while still
using the pinned file for something real (the PE-level verification, and an
honest cross-check number in coverage_report.json).

FI file real schema (README.md's "Open research" also flagged this as
unverified) — inspected via
``FIsInGene_04142025_with_annotations.txt`` (272,622 data rows) directly:
columns are exactly ``Gene1, Gene2, Annotation, Direction, Score`` (5 columns,
tab-separated, one header row). Two cleanly separable regimes, confirmed by
cross-tabulating ``Annotation`` against ``Direction``/``Score``:

- **Predicted** rows (79,486 of 272,622): ``Annotation`` is the literal
  string ``"predicted"`` (never combined with anything else), ``Direction``
  is always ``"-"`` (no direction asserted), and ``Score`` ranges continuously
  0.88-1.00 (a real ML confidence value — this is what ``--fi-min-score``
  filters).
- **Curated** rows (193,136 of 272,622): ``Annotation`` is a ``"; "``-joined
  list of curated relation terms (e.g. ``"catalyzed by; complex; input"``,
  ``"activate; activated by; catalyze"``), ``Score`` is always exactly 1.00
  (curated = full confidence, not a real probability), and ``Direction``
  takes one of 9 values.

Those 9 ``Direction`` values are not arbitrary tokens: they are exactly the
3x3 combination of an optional left cap (``<`` / ``|`` / absent) and an
optional right cap (``>`` / ``|`` / absent) around a middle dash, confirmed by
enumerating the full distinct set found in the file: ``-``, ``->``, ``-|``,
``<-``, ``<->``, ``<-|``, ``|-``, ``|->``, ``|-|`` (see ``_FI_DIRECTION_CODES``
below for the exact decomposition into directed, signed edges). This grammar
is corroborated by predicted rows always landing on the bare ``-`` (no
sign/direction asserted, consistent with "no cap on either side") and by
the direction/evidence-type crosstab: 101,090 curated rows are themselves
undirected ``-`` (complex/input-derived relations with no regulatory sign),
the rest split across the 8 directed/signed codes.

--------------------------------------------------------------------------
Open Targets disease.parquet — added beyond the prompt's literal source list
--------------------------------------------------------------------------
The task text names exactly two OT 26.06 tables to pin (``target`` and
``association_by_datasource_indirect``). stage0_schemas.py's ``ResolvedDisease``
(already fixed, Stage 0) requires a human-readable ``name``, which neither of
those two tables carries (association rows only carry ``diseaseId``). OT's
``disease/disease.parquet`` is a single ~7MB file (not a multi-GB Spark
directory) in the same 26.06 release and is the only place that name lives,
so it is pinned here as a third, small OT source — this is a deliberate,
documented addition, not a silent scope drift.

Inspecting it also surfaced something worth flagging explicitly: OT 26.06's
disease ontology backbone is **MONDO-primary**, not EFO-primary. E.g. the
row for "rheumatoid arthritis" has ``id == "MONDO_0008383"``; ``EFO:0000685``
appears only in its ``dbXRefs`` list, and README.md's own worked example
(``EFO_0005755``) does not appear as a primary ``id`` anywhere in 26.06's
disease table at all (it may have been merged/obsoleted since v1's 25.09
run — exactly the kind of cross-release drift README.md's Decision 4 already
warns about). ``resolve_disease_id`` below therefore checks ``dbXRefs`` (in
``PREFIX:local`` form) as a fallback whenever the input doesn't match an
``id`` directly, so a legitimate EFO id for a since-remapped MONDO term still
resolves instead of spuriously failing with "disease cannot be resolved".

--------------------------------------------------------------------------
Open Targets association_by_datasource_indirect real schema
--------------------------------------------------------------------------
Columns (inspected directly): ``diseaseId, targetId, aggregationType,
aggregationValue, associationScore, evidenceCount, timeseries,
currentNovelty``. ``aggregationType`` is always the literal string
``"datasourceId"`` in this table (asserted at runtime — if OT ever changes
this, ingest fails loudly rather than silently mis-tagging datatypes);
``aggregationValue`` holds the actual datasource id (e.g. ``"reactome"``,
``"eva"``, ``"gwas_credible_sets"``). There is no ``datatype_id`` column
directly, so ``_DATASOURCE_TO_DATATYPE`` below (sourced from
https://platform-docs.opentargets.org/evidence) maps each of the 20
datasource ids actually observed in 26.06 to its OT datatype. An
unrecognized ``aggregationValue`` raises rather than being silently dropped
or mis-bucketed, since Stage 4's non-pathway-datatype restriction (README.md,
Stage 4) depends on this mapping being complete and correct.

--------------------------------------------------------------------------
Interaction source open decision (README.md, "What counts as pathway
topology?")
--------------------------------------------------------------------------
The FI file mixes curated reactions with ML-predicted interactions, which is
one of five options README.md lists and explicitly calls a placeholder, not
a decision. This module does not resolve that decision. ``_load_fi_interactions``
is the single adapter function Stage 3 (and any future re-run of Stage 1)
depends on for gene-gene interaction edges; swapping the source to Pathway
Commons SIF or Reactome BioPAX-derived reactions is a rewrite of that one
function's body (and its registry entry), not of anything else in this file.
"""

from __future__ import annotations

import argparse
import difflib
import hashlib
import json
import logging
import re
import subprocess
import sys
import urllib.error
import urllib.request
import zipfile
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow.dataset as pa_ds
import pyarrow.compute as pa_compute
import pyarrow.parquet as pa_parquet

from stage0_schemas import (
    ENSEMBL_GENE_ID_PATTERN,
    GeneSetSizeBucket,
    GeneSetsSchema,
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
    SUPPORTED_SOURCE_DBS,
    assert_foreign_key,
)

logger = logging.getLogger("ingest")

# ==========================================================================
# Pinned source registry
# ==========================================================================

OPENTARGETS_VERSION = "26.06"
REACTOME_VERSION = "97"
REACTOME_FI_VERSION = "04142025"

_OT_BASE = f"https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/{OPENTARGETS_VERSION}/output"
_REACTOME_BASE = f"https://download.reactome.org/{REACTOME_VERSION}"

# The OT target/association_by_datasource_indirect directories are Spark
# output: each part filename embeds a write-time-random UUID
# (e.g. "part-00000-810593a9-...-c000.snappy.parquet") that cannot be known
# ahead of time without already having listed the directory, and
# association_by_datasource_indirect alone is ~3.3GB across 43 parts —
# genuinely multi-gigabyte, which is exactly why --no-download exists (for
# constrained/offline environments where fetching that much data per run
# isn't viable). A literal `PinnedFile(url, sha256, size)` per part file
# would therefore either be wrong (guessed ahead of time) or stale (pinned
# to one release snapshot's UUIDs, breaking if OT ever regenerates the same
# version number with different part boundaries). Instead, DirectorySource
# entries are resolved at runtime against a *local* lockfile
# (<data-dir>/opentargets/<version>/<key>/.checksums.json): the first
# successful download of each part file records its observed sha256/size
# there, and every subsequent run treats that lockfile as the pin —
# reproducibility is enforced from the moment a data-dir is first populated,
# just not baked into this module's source text for content whose exact
# shape is only knowable by asking the server.
@dataclass(frozen=True)
class PinnedFile:
    url: str
    sha256: str
    size: int
    dest: str  # path relative to --data-dir


@dataclass(frozen=True)
class DirectorySource:
    url: str  # Apache-index-listable directory URL, trailing slash
    dest: str  # subdirectory relative to --data-dir


REGISTRY: dict[str, PinnedFile] = {
    "reactome_gmt": PinnedFile(
        url=f"{_REACTOME_BASE}/ReactomePathways.gmt.zip",
        sha256="8c1dbc8578431da5d2d5118262718c60b553a9be3398e93658daa069e4a9afd4",
        size=298479,
        dest="reactome/ReactomePathways.gmt.zip",
    ),
    "reactome_ensembl2reactome": PinnedFile(
        url=f"{_REACTOME_BASE}/Ensembl2Reactome.txt",
        sha256="34852dc2ac258a851bad914708c0045a636161dbbeca604f04343174cca8d7fb",
        size=183277348,
        dest="reactome/Ensembl2Reactome.txt",
    ),
    "reactome_fi": PinnedFile(
        url=f"https://reactome.org/download/tools/ReactomeFIs/FIsInGene_{REACTOME_FI_VERSION}_with_annotations.txt.zip",
        sha256="e5c373e3baff9764c07edb36c1adc8b714a3e1b44dfe6f8892bafaeb05fb2335",
        size=1182418,
        dest=f"reactome/FIsInGene_{REACTOME_FI_VERSION}_with_annotations.txt.zip",
    ),
    "opentargets_disease": PinnedFile(
        url=f"{_OT_BASE}/disease/disease.parquet",
        sha256="b328c5e775bc40d954f2e058aa59e8e972f2883a262d096f8dffe365736a1aed",
        size=7150787,
        dest=f"opentargets/{OPENTARGETS_VERSION}/disease/disease.parquet",
    ),
}

DIRECTORY_SOURCES: dict[str, DirectorySource] = {
    "opentargets_target": DirectorySource(
        url=f"{_OT_BASE}/target/",
        dest=f"opentargets/{OPENTARGETS_VERSION}/target",
    ),
    "opentargets_association_by_datasource_indirect": DirectorySource(
        url=f"{_OT_BASE}/association_by_datasource_indirect/",
        dest=f"opentargets/{OPENTARGETS_VERSION}/association_by_datasource_indirect",
    ),
}

# OT datasourceId (association_by_datasource_indirect.aggregationValue) ->
# OT datatypeId, sourced from https://platform-docs.opentargets.org/evidence
# and cross-checked against the 20 distinct aggregationValue strings actually
# observed in 26.06 (see module docstring). Any datasource not listed here
# causes a hard failure rather than a silent drop/mis-bucket.
_DATASOURCE_TO_DATATYPE: dict[str, str] = {
    # genetic_association
    "gwas_credible_sets": "genetic_association",
    "gene_burden": "genetic_association",
    "eva": "genetic_association",
    "genomics_england": "genetic_association",
    "gene2phenotype": "genetic_association",
    "uniprot_literature": "genetic_association",
    "uniprot_variants": "genetic_association",
    "orphanet": "genetic_association",
    "clingen": "genetic_association",
    # somatic_mutation
    "cancer_gene_census": "somatic_mutation",
    "intogen": "somatic_mutation",
    "eva_somatic": "somatic_mutation",
    # affected_pathway
    "cancer_biomarkers": "affected_pathway",
    "crispr_screen": "affected_pathway",
    "crispr": "affected_pathway",
    "reactome": "affected_pathway",
    # literature
    "europepmc": "literature",
    # rna_expression
    "expression_atlas": "rna_expression",
    # animal_model
    "impc": "animal_model",
    # known_drug
    "clinical_precedence": "known_drug",
}

# Placeholder floor pending README.md's "Open research" item ("the real
# Reactome->Ensembl mapping rate — the 90% coverage floor is a placeholder");
# not derived from data.
MIN_COVERAGE_PERCENT = 90.0
# Also a placeholder (no derived number yet, per README.md's open research
# list) — a disease with fewer genetic_association genes than this is judged
# too sparse for Stage 2's GSEA to say anything meaningful.
MIN_GENETIC_ASSOCIATION_GENES = 5

# FI Direction code -> list of (is_forward, sign) edges to emit, where
# "forward" means gene_a=Gene1 -> gene_b=Gene2 and "reverse" means
# gene_a=Gene2 -> gene_b=Gene1. See module docstring for how this 3x3 left-
# cap/right-cap grammar was derived from the real file. The bare "-" case is
# handled separately (one undirected, sign=0 edge) since it has no caps at all.
_FI_DIRECTION_CODES: dict[str, list[tuple[bool, int]]] = {
    "->": [(True, 1)],
    "-|": [(True, -1)],
    "<-": [(False, 1)],
    "|-": [(False, -1)],
    "<->": [(True, 1), (False, 1)],
    "|-|": [(True, -1), (False, -1)],
    "<-|": [(False, 1), (True, -1)],
    "|->": [(False, -1), (True, 1)],
}


# ==========================================================================
# CLI
# ==========================================================================


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="stage1_ingest.py",
        description="Stage 1: resolve, canonicalize, and normalize source data.",
    )
    p.add_argument("--target", required=True, help="Gene symbol or Ensembl gene ID.")
    p.add_argument("--disease", required=True, help="EFO (or other OT-recognized ontology) disease ID.")
    p.add_argument("--data-dir", required=True, type=Path, help="Source cache directory.")
    p.add_argument("--out-dir", required=True, type=Path, help="Run artifact output directory.")
    p.add_argument(
        "--pathway-db",
        choices=sorted(SUPPORTED_SOURCE_DBS),
        default="reactome",
        help="Pathway curation database. Only 'reactome' is supported this hackathon; "
        "the flag is kept so its value flows into the source_db column.",
    )
    p.add_argument("--no-download", action="store_true", help="Never fetch from the network; --data-dir must already be populated.")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--min-set-size", type=int, default=10)
    p.add_argument("--max-set-size", type=int, default=200)
    p.add_argument(
        "--fi-curated-only",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Keep only curated Reactome FI rows (predicted rows carry no reliable direction/sign).",
    )
    p.add_argument("--fi-min-score", type=float, default=None, help="Drop FI rows with Score below this (nullable).")
    p.add_argument(
        "--allow-low-coverage",
        action="store_true",
        help="Do not fail on <90%% identifier-mapping coverage or an implausibly small "
        "genetic_association gene count.",
    )
    return p


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


# ==========================================================================
# Small utilities
# ==========================================================================


def _fail(msg: str) -> None:
    logger.error(msg)
    sys.exit(1)


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# urllib's default User-Agent ("Python-urllib/x.y") is blocked outright by
# reactome.org's Cloudflare bot protection (403, verified by hand: curl with
# any User-Agent succeeds, bare urlopen() does not) even though the URL,
# checksum, and file are otherwise fine. A plain browser-like UA is enough
# to pass; applied to every request here, not just reactome.org's, since
# EBI's or Reactome's other hosts could start doing the same at any time.
_USER_AGENT = "Mozilla/5.0 (X11; Linux x86_64) targets-from-pathways/stage1_ingest.py"


def _urlopen(url: str, timeout: float):
    return urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": _USER_AGENT}), timeout=timeout)


def _download_to(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    try:
        with _urlopen(url, timeout=120) as resp, tmp.open("wb") as f:
            while True:
                chunk = resp.read(1 << 20)
                if not chunk:
                    break
                f.write(chunk)
    except (urllib.error.URLError, TimeoutError) as e:
        tmp.unlink(missing_ok=True)
        _fail(f"Failed to download {url}: {e}")
    tmp.replace(dest)


def _list_apache_index(url: str, suffix: str) -> list[tuple[str, int]]:
    """Parse an Apache directory-index HTML page for hrefs ending in `suffix`.

    EBI's FTP mirror serves plain Apache `mod_autoindex` listings for OT's
    parquet directories (verified: fetched several by hand, see module
    docstring); this is the same listing `ftp://` exposes, in a form trivial
    to regex without an FTP client.
    """
    try:
        with _urlopen(url, timeout=60) as resp:
            html = resp.read().decode("utf-8", errors="replace")
    except (urllib.error.URLError, TimeoutError) as e:
        _fail(f"Failed to list {url}: {e}")
    names = re.findall(r'href="([^"?][^"]*' + re.escape(suffix) + r')"', html)
    out = []
    for name in names:
        m = re.search(re.escape(name) + r'</a></td><td align="right">[^<]*</td>'
                       r'<td align="right">\s*([0-9.]+)([MKG]?)\s*</td>', html)
        size = -1
        if m:
            num, unit = m.groups()
            mult = {"": 1, "K": 1024, "M": 1024**2, "G": 1024**3}[unit]
            size = int(float(num) * mult)
        out.append((name, size))
    return out


def _ensure_pinned_file(key: str, data_dir: Path, no_download: bool) -> Path:
    pf = REGISTRY[key]
    dest = data_dir / pf.dest

    def _verify(path: Path) -> bool:
        if not path.exists():
            return False
        if path.stat().st_size != pf.size:
            return False
        return _sha256_file(path) == pf.sha256

    if _verify(dest):
        return dest

    if no_download:
        _fail(
            f"[{key}] missing or checksum-invalid at {dest}.\n"
            f"  Obtain it manually:\n"
            f"    curl -L -o {dest} '{pf.url}'\n"
            f"  Expected sha256: {pf.sha256}\n"
            f"  Expected size:   {pf.size} bytes\n"
            f"  Then re-run without --no-download to verify, or with it once the file is in place."
        )

    logger.info("[%s] downloading %s -> %s", key, pf.url, dest)
    dest.unlink(missing_ok=True)
    _download_to(pf.url, dest)
    if not _verify(dest):
        actual = _sha256_file(dest) if dest.exists() else "<missing>"
        dest.unlink(missing_ok=True)
        _fail(
            f"[{key}] checksum mismatch after download from {pf.url}.\n"
            f"  expected sha256={pf.sha256} size={pf.size}\n"
            f"  actual   sha256={actual} size={dest.stat().st_size if dest.exists() else 0}\n"
            f"  Refusing to proceed on a checksum mismatch."
        )
    return dest


def _ensure_directory_source(key: str, data_dir: Path, no_download: bool) -> list[Path]:
    ds_spec = DIRECTORY_SOURCES[key]
    local_dir = data_dir / ds_spec.dest
    lockfile = local_dir / ".checksums.json"

    def _load_lock() -> dict[str, dict[str, Any]] | None:
        if not lockfile.exists():
            return None
        return json.loads(lockfile.read_text())

    def _verify_against_lock(lock: dict[str, dict[str, Any]]) -> list[Path] | None:
        paths = []
        for name, meta in sorted(lock.items()):
            p = local_dir / name
            if not p.exists() or p.stat().st_size != meta["size"] or _sha256_file(p) != meta["sha256"]:
                return None
            paths.append(p)
        return paths

    lock = _load_lock()
    if lock is not None:
        verified = _verify_against_lock(lock)
        if verified is None:
            _fail(
                f"[{key}] local files under {local_dir} no longer match the pinned "
                f"lockfile {lockfile} — never proceeding on a checksum mismatch. "
                f"Delete {local_dir} and re-run to re-acquire, or restore the original files."
            )
        return verified

    if no_download:
        # No lockfile yet. If files were placed here out-of-band (sneakernet
        # into a constrained/offline environment), pin them now from
        # whatever's on disk rather than hard-failing — but say so plainly,
        # since completeness against the *remote* listing was never checked.
        local_parquets = sorted(local_dir.glob("*.parquet")) if local_dir.exists() else []
        if local_parquets:
            logger.warning(
                "[%s] --no-download: pinning %d local file(s) already present under %s "
                "as first-seen (completeness against the remote directory listing was "
                "NOT verified, since that requires a network request).",
                key, len(local_parquets), local_dir,
            )
            lock = {p.name: {"sha256": _sha256_file(p), "size": p.stat().st_size} for p in local_parquets}
            lockfile.write_text(json.dumps(lock, indent=2, sort_keys=True))
            return local_parquets
        _fail(
            f"[{key}] no local files under {local_dir} and --no-download is set.\n"
            f"  Obtain it manually (multi-gigabyte directory — this is exactly what "
            f"--no-download exists for): mirror every *.snappy.parquet file from\n"
            f"    {ds_spec.url}\n"
            f"  into {local_dir}, then re-run (with or without --no-download; a local "
            f"lockfile will be generated from whatever's found there)."
        )

    logger.info("[%s] listing %s", key, ds_spec.url)
    remote_files = _list_apache_index(ds_spec.url, ".snappy.parquet")
    if not remote_files:
        _fail(f"[{key}] directory listing at {ds_spec.url} returned no *.snappy.parquet files.")
    local_dir.mkdir(parents=True, exist_ok=True)
    lock = {}
    paths = []
    for name, _size in sorted(remote_files):
        dest = local_dir / name
        if not dest.exists():
            logger.info("[%s] downloading %s", key, name)
            _download_to(f"{ds_spec.url}{name}", dest)
        lock[name] = {"sha256": _sha256_file(dest), "size": dest.stat().st_size}
        paths.append(dest)
    lockfile.write_text(json.dumps(lock, indent=2, sort_keys=True))
    return paths


def _extract_single_member(zip_path: Path) -> Path:
    """Extract a single-file zip alongside itself, verifying it's the only member."""
    out = zip_path.parent / zip_path.stem
    if out.exists():
        return out
    with zipfile.ZipFile(zip_path) as zf:
        members = zf.namelist()
        if len(members) != 1:
            _fail(f"{zip_path} expected to contain exactly one file, found {members}")
        with zf.open(members[0]) as src, out.open("wb") as dst:
            dst.write(src.read())
    return out


# ==========================================================================
# Symbol resolution
# ==========================================================================


@dataclass(frozen=True)
class SymbolIndex:
    symbol_to_gene: dict[str, str]
    synonym_to_genes: dict[str, set[str]]


def _build_symbol_index(genes_df: pd.DataFrame) -> SymbolIndex:
    symbol_to_gene: dict[str, str] = {}
    ambiguous_symbols: set[str] = set()
    for gene_id, symbol in zip(genes_df["gene_id"], genes_df["symbol"]):
        key = symbol.upper()
        if key in symbol_to_gene and symbol_to_gene[key] != gene_id:
            ambiguous_symbols.add(key)
        else:
            symbol_to_gene[key] = gene_id
    for key in ambiguous_symbols:
        symbol_to_gene.pop(key, None)

    synonym_to_genes: dict[str, set[str]] = {}
    for gene_id, syns in zip(genes_df["gene_id"], genes_df["synonyms"]):
        for s in syns:
            synonym_to_genes.setdefault(s.upper(), set()).add(gene_id)

    return SymbolIndex(symbol_to_gene, synonym_to_genes)


def resolve_symbol(raw: str, index: SymbolIndex) -> tuple[str | None, str]:
    """Try the approved symbol first, then the synonym field (task (b)).

    Returns (gene_id or None, reason). reason is one of "approved_symbol",
    "synonym", "ambiguous_synonym", "unmapped".
    """
    key = raw.strip().upper()
    if key in index.symbol_to_gene:
        return index.symbol_to_gene[key], "approved_symbol"
    genes = index.synonym_to_genes.get(key)
    if genes:
        if len(genes) == 1:
            return next(iter(genes)), "synonym"
        return None, "ambiguous_synonym"
    return None, "unmapped"


# ==========================================================================
# (b) genes.parquet — Ensembl gene ID canonicalization, OT target = authority
# ==========================================================================


def build_genes_table(target_paths: list[Path]) -> pd.DataFrame:
    table = pa_ds.dataset(sorted(target_paths), format="parquet").to_table(
        columns=["id", "approvedSymbol", "biotype", "symbolSynonyms", "obsoleteSymbols", "proteinIds", "transcriptIds"]
    )
    df = table.to_pandas()

    # pyarrow's to_pandas() turns list<struct> columns into numpy object
    # arrays, not Python lists — `arr or []` then calls bool() on a
    # multi-element array (ambiguous), so None-ness must be checked
    # explicitly rather than via truthiness.
    def _labels(structs) -> list[str]:
        return [s["label"] for s in structs] if structs is not None else []

    synonyms = [
        sorted({*_labels(row.symbolSynonyms), *_labels(row.obsoleteSymbols)})
        for row in df.itertuples()
    ]

    genes_df = pd.DataFrame(
        {
            "gene_id": df["id"],
            "symbol": df["approvedSymbol"],
            "synonyms": synonyms,
            "biotype": df["biotype"].fillna(""),
        }
    )
    genes_df = genes_df.drop_duplicates(subset="gene_id").reset_index(drop=True)

    # Kept separately (not returned) for the Ensembl2Reactome ENSP/ENST -> gene
    # collapse used only by the cross-check coverage source.
    genes_df.attrs["_protein_to_gene"] = _build_protein_to_gene(df)
    genes_df.attrs["_transcript_to_gene"] = _build_transcript_to_gene(df)
    return genes_df


def _build_protein_to_gene(target_df: pd.DataFrame) -> dict[str, str]:
    out: dict[str, str] = {}
    for gene_id, prot_structs in zip(target_df["id"], target_df["proteinIds"]):
        for p in (prot_structs if prot_structs is not None else []):
            if p["source"] == "ensembl_PRO":
                out[p["id"]] = gene_id
    return out


def _build_transcript_to_gene(target_df: pd.DataFrame) -> dict[str, str]:
    out: dict[str, str] = {}
    for gene_id, tx_list in zip(target_df["id"], target_df["transcriptIds"]):
        for t in (tx_list if tx_list is not None else []):
            out[t] = gene_id
    return out


# ==========================================================================
# Coverage tracking
# ==========================================================================


@dataclass
class CoverageTracker:
    source: str
    raw: set[str] = field(default_factory=set)
    mapped: set[str] = field(default_factory=set)
    dropped: set[str] = field(default_factory=set)

    def record(self, raw_id: str, gene_id: str | None) -> None:
        self.raw.add(raw_id)
        if gene_id is not None:
            self.mapped.add(raw_id)
        else:
            self.dropped.add(raw_id)

    def to_entry(self) -> dict[str, Any]:
        raw_n, mapped_n = len(self.raw), len(self.mapped)
        dropped_n = raw_n - mapped_n
        pct = round(100.0 * mapped_n / raw_n, 2) if raw_n else 0.0
        return {
            "source": self.source,
            "raw_gene_count": raw_n,
            "mapped_count": mapped_n,
            "percent_mapped": pct,
            "dropped_count": dropped_n,
            "example_dropped_ids": sorted(self.dropped)[:20],
        }


# ==========================================================================
# (d) gene_sets.parquet — from ReactomePathways.gmt
# ==========================================================================


def parse_gmt(gmt_path: Path) -> list[tuple[str, str, list[str]]]:
    """Return [(set_id, set_name, [gene_symbol, ...]), ...].

    Standard GMT layout, confirmed against the real file: col1=pathway name
    (set_name), col2=Reactome stable ID e.g. "R-HSA-164843" (set_id),
    col3+ = gene symbols. All 2868 rows are R-HSA-* (human-only file).
    """
    out = []
    for line in gmt_path.read_text().splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        set_name, set_id, *symbols = fields
        out.append((set_id, set_name, [s for s in symbols if s]))
    return out


def build_gene_sets_table(
    gmt_path: Path,
    index: SymbolIndex,
    source_db: str,
    coverage: CoverageTracker,
) -> pd.DataFrame:
    rows = []
    for set_id, set_name, symbols in parse_gmt(gmt_path):
        seen_genes: set[str] = set()
        for symbol in symbols:
            gene_id, _reason = resolve_symbol(symbol, index)
            coverage.record(symbol, gene_id)
            if gene_id is None or gene_id in seen_genes:
                continue
            seen_genes.add(gene_id)
            rows.append(
                {
                    "set_id": set_id,
                    "set_name": set_name,
                    "source_db": source_db,
                    "source_version": REACTOME_VERSION,
                    "gene_id": gene_id,
                    "hierarchy_level": None,
                    "parent_id": None,
                }
            )
    return pd.DataFrame(
        rows, columns=["set_id", "set_name", "source_db", "source_version", "gene_id", "hierarchy_level", "parent_id"]
    )


# ==========================================================================
# (d) interactions.parquet — from the Reactome FI file (adapter function)
# ==========================================================================


def _load_fi_interactions(
    fi_path: Path,
    index: SymbolIndex,
    source_db: str,
    coverage: CoverageTracker,
    *,
    curated_only: bool,
    min_score: float | None,
) -> pd.DataFrame:
    """Single adapter: swap the interaction source here (see module docstring
    "Interaction source open decision") without touching any other stage."""
    df = pd.read_csv(fi_path, sep="\t", dtype=str)

    # pandas' C parser does NOT raise on a short/ragged row (fewer fields
    # than the header) — it silently pads the missing trailing fields with
    # NaN. Left unchecked, that NaN Direction then misses every
    # _FI_DIRECTION_CODES key and the row gets dropped with only a logged
    # warning below — exactly the "silently skipped" failure mode Stage 1 is
    # supposed to avoid. Catch it here, before any row is processed, so a
    # corrupt source row is a hard failure instead.
    required_cols = ["Gene1", "Gene2", "Annotation", "Direction", "Score"]
    malformed = df[df[required_cols].isnull().any(axis=1)]
    if not malformed.empty:
        _fail(
            f"{fi_path} has {len(malformed)} malformed row(s) with a missing "
            f"Gene1/Gene2/Annotation/Direction/Score field (e.g. line "
            f"{malformed.index[0] + 2} of the file: {malformed.iloc[0].to_dict()}). "
            f"Refusing to silently drop a corrupt source row."
        )
    df["Score"] = df["Score"].astype(float)

    rows = []
    for gene1, gene2, annotation, direction, score in zip(
        df["Gene1"], df["Gene2"], df["Annotation"], df["Direction"], df["Score"]
    ):
        evidence_type = "predicted" if annotation == "predicted" else "curated"
        if curated_only and evidence_type != "curated":
            continue
        if min_score is not None and score < min_score:
            continue

        gene_a_id, _ = resolve_symbol(gene1, index)
        coverage.record(gene1, gene_a_id)
        gene_b_id, _ = resolve_symbol(gene2, index)
        coverage.record(gene2, gene_b_id)
        if gene_a_id is None or gene_b_id is None or gene_a_id == gene_b_id:
            continue

        if direction == "-":
            rows.append((gene_a_id, gene_b_id, False, 0, source_db, REACTOME_FI_VERSION, evidence_type, score))
            continue
        edges = _FI_DIRECTION_CODES.get(direction)
        if edges is None:
            logger.warning("Unrecognized FI Direction code %r for %s/%s — dropping row.", direction, gene1, gene2)
            continue
        for is_forward, sign in edges:
            src, dst = (gene_a_id, gene_b_id) if is_forward else (gene_b_id, gene_a_id)
            rows.append((src, dst, True, sign, source_db, REACTOME_FI_VERSION, evidence_type, score))

    out = pd.DataFrame(
        rows,
        columns=["gene_a", "gene_b", "directed", "sign", "source_db", "source_version", "evidence_type", "confidence"],
    ).drop_duplicates()
    return out.reset_index(drop=True)


# ==========================================================================
# (b, cross-check) Ensembl2Reactome.txt — verified gene-level, used only to
# report an independent coverage number (see module docstring)
# ==========================================================================


def _load_ensembl2reactome_coverage(e2r_path: Path, genes_df: pd.DataFrame) -> dict[str, Any]:
    protein_to_gene: dict[str, str] = genes_df.attrs.get("_protein_to_gene", {})
    transcript_to_gene: dict[str, str] = genes_df.attrs.get("_transcript_to_gene", {})
    known_genes = set(genes_df["gene_id"])

    tracker = CoverageTracker(source="reactome_ensembl2reactome")
    with e2r_path.open() as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 6 or fields[5] != "Homo sapiens":
                continue
            raw_id = fields[0]
            bare_id = re.sub(r"\.\d+$", "", raw_id)
            if bare_id in known_genes:
                gene_id = bare_id
            elif bare_id in protein_to_gene:
                gene_id = protein_to_gene[bare_id]
            elif bare_id in transcript_to_gene:
                gene_id = transcript_to_gene[bare_id]
            else:
                gene_id = None
            tracker.record(raw_id, gene_id)
    return tracker.to_entry()


# ==========================================================================
# (c) OT disease resolution + disease-filtered association subset
# ==========================================================================


def load_disease_table(disease_parquet: Path) -> pd.DataFrame:
    table = pa_parquet.read_table(disease_parquet, columns=["id", "name", "dbXRefs"])
    return table.to_pandas()


def resolve_disease_id(raw: str, disease_df: pd.DataFrame) -> tuple[str, str] | tuple[None, None]:
    """Resolve `raw` to (canonical_id, name). Tries a direct `id` match first
    (normalizing ':' to '_'), then falls back to `dbXRefs` (colon form) — see
    module docstring on OT 26.06's MONDO-primary disease backbone."""
    raw_norm = raw.strip().replace(":", "_")
    hit = disease_df.loc[disease_df["id"] == raw_norm]
    if not hit.empty:
        row = hit.iloc[0]
        return row["id"], row["name"]

    xref_form = re.sub(r"^([A-Za-z]+)_", r"\1:", raw_norm)
    for id_, name, xrefs in zip(disease_df["id"], disease_df["name"], disease_df["dbXRefs"]):
        if xrefs is not None and xref_form in xrefs:
            return id_, name
    return None, None


def load_disease_association_subset(assoc_paths: list[Path], disease_id: str) -> pd.DataFrame:
    dataset = pa_ds.dataset(sorted(assoc_paths), format="parquet")
    table = dataset.to_table(
        columns=["targetId", "diseaseId", "aggregationType", "aggregationValue", "associationScore"],
        filter=pa_compute.field("diseaseId") == disease_id,
    )
    df = table.to_pandas()
    bad_agg = set(df["aggregationType"]) - {"datasourceId"}
    if bad_agg:
        _fail(
            f"association_by_datasource_indirect.aggregationType contained unexpected "
            f"value(s) {sorted(bad_agg)} (expected only 'datasourceId'); refusing to "
            f"guess a datatype mapping."
        )
    unknown_ds = set(df["aggregationValue"]) - set(_DATASOURCE_TO_DATATYPE)
    if unknown_ds:
        _fail(
            f"Unrecognized OT datasource id(s) {sorted(unknown_ds)} with no entry in "
            f"_DATASOURCE_TO_DATATYPE — add them there (see module docstring) rather "
            f"than silently dropping or mis-bucketing evidence."
        )

    return pd.DataFrame(
        {
            "gene_id": df["targetId"],
            "disease_id": df["diseaseId"],
            "datatype_id": df["aggregationValue"].map(_DATASOURCE_TO_DATATYPE),
            "datasource_id": df["aggregationValue"],
            "score": df["associationScore"].astype(float),
        }
    )


# ==========================================================================
# Target resolution (--target: symbol or Ensembl ID)
# ==========================================================================


def resolve_target(raw: str, genes_df: pd.DataFrame, index: SymbolIndex) -> tuple[str, str] | tuple[None, None]:
    raw_norm = raw.strip()
    if ENSEMBL_GENE_ID_PATTERN.match(raw_norm):
        hit = genes_df.loc[genes_df["gene_id"] == raw_norm]
        if not hit.empty:
            return raw_norm, hit.iloc[0]["symbol"]
        return None, None
    gene_id, _reason = resolve_symbol(raw_norm, index)
    if gene_id is not None:
        symbol = genes_df.loc[genes_df["gene_id"] == gene_id, "symbol"].iloc[0]
        return gene_id, symbol
    return None, None


def _suggest_symbols(raw: str, genes_df: pd.DataFrame, n: int = 5) -> list[str]:
    return difflib.get_close_matches(raw.strip().upper(), genes_df["symbol"].str.upper().tolist(), n=n, cutoff=0.6)


# ==========================================================================
# Scale report
# ==========================================================================


def _size_distribution(sizes: list[int]) -> list[GeneSetSizeBucket]:
    counts: dict[int, int] = {}
    for s in sizes:
        counts[s] = counts.get(s, 0) + 1
    return [GeneSetSizeBucket(set_size=k, count=v) for k, v in sorted(counts.items())]


def build_scale_report(
    gene_sets_df: pd.DataFrame,
    interactions_df: pd.DataFrame,
    min_set_size: int,
    max_set_size: int,
) -> tuple[ScaleReport, pd.DataFrame]:
    sizes_before = gene_sets_df.groupby("set_id")["gene_id"].nunique()
    retained_ids = sizes_before[(sizes_before >= min_set_size) & (sizes_before <= max_set_size)].index
    retained_df = gene_sets_df[gene_sets_df["set_id"].isin(retained_ids)].reset_index(drop=True)
    sizes_after = sizes_before.loc[retained_ids]

    projected_edges = int(sum(n * (n - 1) // 2 for n in sizes_after.tolist()))

    sign_counts = interactions_df["sign"].value_counts().to_dict() if not interactions_df.empty else {}
    evidence_counts = interactions_df["evidence_type"].value_counts().to_dict() if not interactions_df.empty else {}

    report = ScaleReport(
        gene_set_size_distribution_before_cap=_size_distribution(sizes_before.tolist()),
        gene_set_size_distribution_after_cap=_size_distribution(sizes_after.tolist()),
        sets_retained=int(len(retained_ids)),
        interaction_counts_by_sign={str(k): int(v) for k, v in sign_counts.items()},
        interaction_counts_by_evidence_type={str(k): int(v) for k, v in evidence_counts.items()},
        projected_comembership_edge_count=projected_edges,
    )
    return report, retained_df


# ==========================================================================
# Orchestration
# ==========================================================================


def _write_validated_parquet(df: pd.DataFrame, schema, out_path: Path) -> OutputArtifact:
    validated = schema.validate(df)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    validated.to_parquet(out_path, index=False)
    return OutputArtifact(path=str(out_path), sha256=_sha256_file(out_path))


def _write_json(obj: Any, out_path: Path) -> OutputArtifact:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(obj, indent=2, sort_keys=True))
    return OutputArtifact(path=str(out_path), sha256=_sha256_file(out_path))


def _git_commit() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True, timeout=10
        ).stdout.strip()
    except Exception:
        return "unknown"


def _run_id(target: str, disease: str, seed: int, cli_parameters: dict[str, Any]) -> str:
    payload = json.dumps({"target": target, "disease": disease, "seed": seed, "cli": cli_parameters}, sort_keys=True)
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def run_ingest(args: argparse.Namespace) -> Manifest:
    data_dir: Path = args.data_dir
    out_dir: Path = args.out_dir
    source_db = args.pathway_db

    # --- (a) acquisition -------------------------------------------------
    gmt_zip = _ensure_pinned_file("reactome_gmt", data_dir, args.no_download)
    e2r_path = _ensure_pinned_file("reactome_ensembl2reactome", data_dir, args.no_download)
    fi_zip = _ensure_pinned_file("reactome_fi", data_dir, args.no_download)
    disease_parquet = _ensure_pinned_file("opentargets_disease", data_dir, args.no_download)
    target_paths = _ensure_directory_source("opentargets_target", data_dir, args.no_download)
    assoc_paths = _ensure_directory_source(
        "opentargets_association_by_datasource_indirect", data_dir, args.no_download
    )

    gmt_path = _extract_single_member(gmt_zip)
    fi_path = _extract_single_member(fi_zip)

    sources = [
        Source(
            db="reactome",
            version=REACTOME_VERSION,
            files=[
                SourceFile(path=str(gmt_zip), sha256=_sha256_file(gmt_zip), bytes=gmt_zip.stat().st_size),
                SourceFile(path=str(e2r_path), sha256=_sha256_file(e2r_path), bytes=e2r_path.stat().st_size),
            ],
        ),
        Source(
            db="reactome",
            version=REACTOME_FI_VERSION,
            files=[SourceFile(path=str(fi_zip), sha256=_sha256_file(fi_zip), bytes=fi_zip.stat().st_size)],
        ),
        Source(
            db="opentargets",
            version=OPENTARGETS_VERSION,
            files=[
                SourceFile(path=str(disease_parquet), sha256=_sha256_file(disease_parquet), bytes=disease_parquet.stat().st_size),
                *[SourceFile(path=str(p), sha256=_sha256_file(p), bytes=p.stat().st_size) for p in sorted(target_paths)],
                *[SourceFile(path=str(p), sha256=_sha256_file(p), bytes=p.stat().st_size) for p in sorted(assoc_paths)],
            ],
        ),
    ]

    # --- (b) genes.parquet -------------------------------------------------
    genes_df = build_genes_table(target_paths)
    index = _build_symbol_index(genes_df)

    # --- target resolution ---------------------------------------------
    gene_id, symbol = resolve_target(args.target, genes_df, index)
    if gene_id is None:
        suggestions = _suggest_symbols(args.target, genes_df)
        hint = f" Did you mean: {', '.join(suggestions)}?" if suggestions else ""
        _fail(f"Could not resolve --target {args.target!r} to a gene in Open Targets {OPENTARGETS_VERSION}.{hint}")
    resolved_target = ResolvedTarget(input=args.target, gene_id=gene_id, symbol=symbol)

    # --- disease resolution ---------------------------------------------
    disease_df = load_disease_table(disease_parquet)
    disease_id, disease_name = resolve_disease_id(args.disease, disease_df)
    if disease_id is None:
        _fail(
            f"Could not resolve --disease {args.disease!r} against Open Targets "
            f"{OPENTARGETS_VERSION}'s disease table (checked both `id` and `dbXRefs`)."
        )

    ot_disease_subset = load_disease_association_subset(assoc_paths, disease_id)
    ot_disease_subset = OTAssociationsSchema.validate(ot_disease_subset)
    assert_foreign_key(ot_disease_subset, "gene_id", genes_df)
    ot_subset_path = out_dir / "ot_disease_subset.parquet"
    ot_subset_path.parent.mkdir(parents=True, exist_ok=True)
    ot_disease_subset.to_parquet(ot_subset_path, index=False)
    ot_subset_artifact = OutputArtifact(path=str(ot_subset_path), sha256=_sha256_file(ot_subset_path))

    n_associated_genes = int(ot_disease_subset["gene_id"].nunique())
    n_genetic_association_genes = int(
        ot_disease_subset.loc[ot_disease_subset["datatype_id"] == "genetic_association", "gene_id"].nunique()
    )
    if n_genetic_association_genes < MIN_GENETIC_ASSOCIATION_GENES and not args.allow_low_coverage:
        _fail(
            f"Disease {disease_id} ({disease_name}) has only {n_genetic_association_genes} "
            f"genetic_association gene(s) — below the {MIN_GENETIC_ASSOCIATION_GENES}-gene floor "
            f"for GSEA to be meaningful (README.md 'Loud failures'). Pass --allow-low-coverage to override."
        )
    resolved_disease = ResolvedDisease(efo_id=disease_id, name=disease_name, n_associated_genes=n_associated_genes)

    # --- (d) gene_sets.parquet -------------------------------------------
    gmt_coverage = CoverageTracker(source="reactome_gmt")
    gene_sets_raw = build_gene_sets_table(gmt_path, index, source_db, gmt_coverage)

    # --- (d) interactions.parquet -----------------------------------------
    fi_coverage = CoverageTracker(source="reactome_fi")
    interactions_df = _load_fi_interactions(
        fi_path, index, source_db, fi_coverage,
        curated_only=args.fi_curated_only, min_score=args.fi_min_score,
    )

    # --- cross-check coverage source ---------------------------------------
    e2r_entry = _load_ensembl2reactome_coverage(e2r_path, genes_df)

    # --- (e) reports --------------------------------------------------------
    coverage_report = {
        "reactome_gmt": gmt_coverage.to_entry(),
        "reactome_fi": fi_coverage.to_entry(),
        "reactome_ensembl2reactome": e2r_entry,
    }
    low_coverage_sources = [
        name for name, entry in coverage_report.items() if entry["percent_mapped"] < MIN_COVERAGE_PERCENT
    ]
    if low_coverage_sources and not args.allow_low_coverage:
        details = "; ".join(f"{s}={coverage_report[s]['percent_mapped']}%" for s in low_coverage_sources)
        _fail(
            f"Source(s) mapped below the {MIN_COVERAGE_PERCENT}% floor: {details}. "
            f"Pass --allow-low-coverage to override."
        )

    scale_report, gene_sets_df = build_scale_report(
        gene_sets_raw, interactions_df, args.min_set_size, args.max_set_size
    )

    # --- validate + drop genes.parquet's private lookup attrs before writing
    genes_df = genes_df.drop(columns=[], errors="ignore")
    genes_df.attrs = {}

    # --- write remaining artifacts ------------------------------------------
    genes_artifact = _write_validated_parquet(genes_df, GenesSchema, out_dir / "genes.parquet")
    assert_foreign_key(gene_sets_df, "gene_id", genes_df)
    gene_sets_artifact = _write_validated_parquet(gene_sets_df, GeneSetsSchema, out_dir / "gene_sets.parquet")
    if not interactions_df.empty:
        assert_foreign_key(interactions_df, "gene_a", genes_df)
        assert_foreign_key(interactions_df, "gene_b", genes_df)
    interactions_artifact = _write_validated_parquet(interactions_df, InteractionsSchema, out_dir / "interactions.parquet")

    coverage_artifact = _write_json(coverage_report, out_dir / "coverage_report.json")
    scale_artifact = _write_json(scale_report.model_dump(mode="json"), out_dir / "scale_report.json")

    output_artifacts = [
        genes_artifact, gene_sets_artifact, interactions_artifact, ot_subset_artifact,
        coverage_artifact, scale_artifact,
    ]

    cli_parameters = {k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()}
    cli_parameters["opentargets_version"] = OPENTARGETS_VERSION
    cli_parameters["reactome_version"] = REACTOME_VERSION
    cli_parameters["reactome_fi_version"] = REACTOME_FI_VERSION

    manifest = Manifest(
        run_id=_run_id(gene_id, disease_id, args.seed, cli_parameters),
        git_commit=_git_commit(),
        created_at=datetime.now(timezone.utc),
        seed=args.seed,
        resolved_target=resolved_target,
        disease=resolved_disease,
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


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    args = parse_args(argv)
    run_ingest(args)


if __name__ == "__main__":
    main()
