"""Stage 6 — Contextual annotation & filtering.

Takes Stage 5's ranked ``candidate_scores.tsv`` and attaches, per candidate
gene, an Open Targets **tractability bucket** and a **safety flag**, then
computes a single ``composite_score`` as a configurable weighted average.

Inputs (all resolved from ``--manifest``'s directory unless overridden):

- ``--candidates``  — Stage 5's ``candidate_scores.tsv`` (``gene``,
  ``topology_score``, optionally ``rwr_score``).
- ``--graph-dir``   — holds Stage 4's ``graph_metadata.json``, whose
  ``genetic_evidence_score`` ``{gene: float}`` map supplies the
  ``genetic_evidence`` component. Optional: a Stage-3-only pipeline never
  ran Stage 4, so this key may be missing — the component is then dropped
  and its weight redistributed (logged).
- ``--target-parquet`` — the Open Targets ``target`` parquet (a directory
  of part files on a real run). Resolved from ``manifest.sources`` by
  default: the ``opentargets`` source files whose path matches
  ``*target*.parquet``. Only three columns are ever read — ``id``,
  ``tractability``, ``safetyLiabilities``.

**No tissue-expression data of any kind is used.** This pipeline is
deliberately pathway- and association-evidence-based only (README.md plan
update 1, hypothesis H10 dropped). Stage 6 does not read
``expression``/``baselineExpression``/GTEx or any tissue column, has no
``--expression`` flag, and emits no expression-derived column.

Tractability bucket (from ``tractability`` = ``list<struct<modality, id,
value>>``; any modality counts):

- ``clinical``  — a ``value == true`` flag whose ``id`` names a clinical
  stage ("Approved Drug", "Phase 1 Clinical", "Advanced Clinical", ...).
- ``discovery`` — some other ``value == true`` flag (structure, ligand,
  pocket, localisation, druggable family).
- ``unknown``  — no ``value == true`` flag, or the gene is absent from the
  target parquet.

Safety flag (from ``safetyLiabilities``, a list):

- ``has_liabilities`` — the list is non-empty (one or more curated
  liabilities).
- ``unknown``         — the list is empty, or the gene is absent.

**Open Targets safety coverage is sparse, and this stage treats it that
way.** There is deliberately no ``safe`` value: OT records known
liabilities or says nothing, never "this gene is safe". A gene with no
liabilities annotation is ``unknown`` (a neutral 0.5 in the composite),
*never* rewarded as if confirmed safe.

Composite (``--weights`` = ``k:v,...``; default
``topology:0.3,rwr:0.3,genetic_evidence:0.2,tractability:0.1,safety:0.1``):
``composite_score = sum(w_k * component_k) / sum(w_k)`` over the components
that have data, where

- ``topology``         = min-max of ``topology_score`` across the candidate
  set (an unbounded proximity number, only meaningful ranked within the
  set);
- ``rwr``              = min-max of ``rwr_score`` across the candidate set;
- ``genetic_evidence`` = ``genetic_evidence_score`` as-is (already an
  absolute [0, 1] Open Targets evidence score);
- ``tractability``     = ``{clinical: 1.0, discovery: 0.5, unknown: 0.0}``;
- ``safety``           = ``{has_liabilities: 0.0, unknown: 0.5}``.

Every component is in [0, 1], so ``composite_score`` is too. Weights need
not sum to 1 (the divisor is the sum of the weights actually used). A
weighted component with no data (``rwr`` when Stage 5 ran topology-only;
``genetic_evidence`` with no Stage 4 metadata) is dropped and the rest
renormalised, with a log line.

**Output** — ``candidates_annotated.tsv`` (or ``--out``), validated against
``schemas.AnnotatedCandidatesSchema``, ranked by ``composite_score``
descending.

Per-candidate **evidence trace** (so Stage 8's report is interpretable, not
just a final number):

- the raw component scores it stands on — ``topology_score``,
  ``rwr_score`` (when Stage 5 produced it), ``genetic_evidence_score``
  (when Stage 4 ran);
- the annotation buckets and their evidence — ``tractability``,
  ``safety``, ``n_safety_liabilities``;
- ``composite_breakdown`` — ``k=contribution|...``, each component's
  *weighted contribution* to ``composite_score`` (the terms sum to it), so
  a reader sees exactly what drove the rank;
- ``composite_weights`` — the renormalised weights actually used
  (``k:fraction,...``, sum 1), repeated on every row so one row is a
  self-contained explanation.

What the trace does **not** carry, by design:

- *which shared pathways / interactions* put a candidate near the target —
  that provenance is not in this stage's inputs (Stage 5 collapses the
  graph to scores), and re-deriving Stage 3's pathway union here would
  duplicate it. It is Stage 8's drill-down, assembled from
  ``gene_sets.parquet`` / ``interactions.parquet`` at render time (see
  PROMPTS.md Stage 8).
- *which Open Targets datatype* produced ``genetic_evidence_score`` — Stage
  4 records only the final score; surfacing the winning datatype is a
  Stage 4 follow-up.

No rows are filtered out; "& filtering" in the stage name is left to Stage
8's UI, which can threshold on any of these columns.

Production note: a cleaner design would have Stage 1 emit a
``target_annotations.parquet`` artifact so this stage reads only Stage 1
outputs (the pipeline's usual discipline). Reading the OT ``target`` parquet
directly here is a deliberate prototype shortcut — the extraction is a
dozen lines and can move upstream later without changing this stage's
output.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import pandas as pd
import pyarrow.compute as pc
import pyarrow.dataset as pa_ds

from stage0_schemas import (
    AnnotatedCandidatesSchema,
    CandidateScoresSchema,
    ENSEMBL_GENE_ID_PATTERN,
    Manifest,
)

logger = logging.getLogger("annotate_context")

DEFAULT_WEIGHTS = "topology:0.3,rwr:0.3,genetic_evidence:0.2,tractability:0.1,safety:0.1"
WEIGHTABLE_COMPONENTS = ("topology", "rwr", "genetic_evidence", "tractability", "safety")

# Ordinal proxies for the two categorical annotations (see module docstring).
TRACTABILITY_SCORE = {"clinical": 1.0, "discovery": 0.5, "unknown": 0.0}
SAFETY_SCORE = {"has_liabilities": 0.0, "unknown": 0.5}

# tractability `id` values that mean "in the clinic" (matched case-insensitively).
_CLINICAL_TRACTABILITY_IDS = {
    "approved drug",
    "phase 1 clinical",
    "phase 2 clinical",
    "phase 3 clinical",
    "advanced clinical",
    "advanced clinical trials",
}


# ==========================================================================
# CLI
# ==========================================================================


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="stage6_annotate_context.py",
        description="Stage 6: attach Open Targets tractability/safety to Stage 5's candidates "
        "and compute a configurable composite score.",
    )
    p.add_argument("--manifest", required=True, type=Path, help="Path to Stage 1's manifest.json.")
    p.add_argument(
        "--candidates",
        type=Path,
        default=None,
        help="Stage 5's candidate_scores.tsv. Default: alongside --manifest.",
    )
    p.add_argument(
        "--graph-dir",
        type=Path,
        default=None,
        help="Directory with Stage 4's graph_metadata.json (for genetic_evidence_score). "
        "Default: --manifest's directory.",
    )
    p.add_argument(
        "--target-parquet",
        type=Path,
        default=None,
        help="Open Targets `target` parquet file or directory. Default: resolved from "
        "manifest.sources (the opentargets *target*.parquet files).",
    )
    p.add_argument(
        "--weights",
        default=DEFAULT_WEIGHTS,
        help=f"Comma-separated component:weight pairs. Default: {DEFAULT_WEIGHTS}. "
        f"Components: {', '.join(WEIGHTABLE_COMPONENTS)}. Weights need not sum to 1.",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output TSV path. Default: <candidates dir>/candidates_annotated.tsv.",
    )
    return p


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    return build_arg_parser().parse_args(argv)


def _fail(msg: str) -> None:
    logger.error(msg)
    sys.exit(1)


def parse_weights(raw: str) -> dict[str, float]:
    weights: dict[str, float] = {}
    for token in (t.strip() for t in raw.split(",") if t.strip()):
        if ":" not in token:
            _fail(f"--weights entry {token!r} is not `component:weight`.")
        key, _, value = token.partition(":")
        key = key.strip().lower()
        if key not in WEIGHTABLE_COMPONENTS:
            _fail(f"--weights has unknown component {key!r}; valid: {list(WEIGHTABLE_COMPONENTS)}.")
        try:
            w = float(value)
        except ValueError:
            _fail(f"--weights entry {token!r} has a non-numeric weight.")
        if w < 0:
            _fail(f"--weights {key!r} is negative ({w}).")
        weights[key] = w
    if not weights or sum(weights.values()) <= 0:
        _fail("--weights resolved to nothing positive.")
    return weights


# ==========================================================================
# Input resolution
# ==========================================================================


def load_manifest(manifest_path: Path) -> Manifest:
    return Manifest.model_validate(json.loads(manifest_path.read_text()))


def _resolve_beside_manifest(path: Path, manifest_path: Path, what: str) -> Path:
    """Return ``path`` if it exists, else the same basename next to the
    manifest, else fail — Stage 1's out-dir paths are relative to whatever
    cwd it ran in, which may not be this process's."""
    if path.exists():
        return path
    beside = manifest_path.parent / path.name
    if beside.exists():
        return beside
    _fail(f"{what} not found at {path} (also checked {beside}).")


def resolve_target_parquet(
    manifest: Manifest, manifest_path: Path, override: Path | None
) -> list[Path]:
    if override is not None:
        if override.is_dir():
            parts = sorted(override.glob("*.parquet"))
            if not parts:
                _fail(f"--target-parquet {override} is a directory with no .parquet files.")
            return parts
        if not override.exists():
            _fail(f"--target-parquet {override} does not exist.")
        return [override]

    # Real OT layout is a `target/` *directory* of `part-*.parquet` files, so
    # "target" lives in the path, not the basename — match on the full path,
    # and exclude the sibling disease/association parquet trees.
    candidates: list[Path] = []
    for source in manifest.sources:
        if source.db != "opentargets":
            continue
        for f in source.files:
            path_l = str(f.path).lower()
            if not path_l.endswith(".parquet") or "target" not in path_l:
                continue
            if "disease" in path_l or "association" in path_l:
                continue
            candidates.append(Path(f.path))
    if not candidates:
        _fail(
            "manifest.sources lists no opentargets *target*.parquet file — pass --target-parquet "
            "explicitly (the Open Targets `target` parquet directory)."
        )
    return [_resolve_beside_manifest(p, manifest_path, "target parquet part") for p in candidates]


def load_candidates(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    return CandidateScoresSchema.validate(df)


def load_genetic_evidence(graph_dir: Path) -> dict[str, float] | None:
    meta_path = graph_dir / "graph_metadata.json"
    if not meta_path.exists():
        logger.warning("No graph_metadata.json in %s — 'genetic_evidence' component unavailable.", graph_dir)
        return None
    scores = json.loads(meta_path.read_text()).get("genetic_evidence_score")
    if not scores:
        logger.warning(
            "graph_metadata.json has no genetic_evidence_score (Stage 4 not run?) — "
            "'genetic_evidence' component unavailable."
        )
        return None
    return {g: float(s) for g, s in scores.items()}


# ==========================================================================
# Open Targets tractability / safety
# ==========================================================================


def _bucket_tractability(entries) -> str:
    try:
        empty = entries is None or len(entries) == 0
    except TypeError:  # a scalar null rather than a list
        empty = True
    if empty:
        return "unknown"
    true_ids = [str(e["id"]).lower() for e in entries if bool(e["value"])]
    if not true_ids:
        return "unknown"
    if any(tid in _CLINICAL_TRACTABILITY_IDS for tid in true_ids):
        return "clinical"
    return "discovery"


def _safety(entries) -> tuple[str, int]:
    try:
        n = len(entries) if entries is not None else 0
    except TypeError:  # a scalar null (NaN) rather than a list
        n = 0
    return ("has_liabilities" if n > 0 else "unknown"), int(n)


def build_annotations(target_parquet_paths: list[Path], genes: list[str]) -> pd.DataFrame:
    """One row per gene in ``genes``: tractability bucket, safety flag,
    liability count. Genes absent from the target parquet come back as
    unknown/unknown/0 — never dropped, never penalised."""
    wanted = set(genes)
    table = pa_ds.dataset([str(p) for p in target_parquet_paths], format="parquet").to_table(
        columns=["id", "tractability", "safetyLiabilities"],
        filter=pc.field("id").isin(list(wanted)),  # pushed down: only the candidate rows are read
    )
    raw = table.to_pandas()

    rows = []
    for r in raw.itertuples(index=False):
        flag, n = _safety(r.safetyLiabilities)
        rows.append(
            {
                "gene": r.id,
                "tractability": _bucket_tractability(r.tractability),
                "safety": flag,
                "n_safety_liabilities": n,
            }
        )
    annotated = pd.DataFrame(rows, columns=["gene", "tractability", "safety", "n_safety_liabilities"])
    annotated["n_safety_liabilities"] = annotated["n_safety_liabilities"].astype(int)

    missing = sorted(wanted - set(annotated["gene"]))
    if missing:
        filler = pd.DataFrame(
            {
                "gene": missing,
                "tractability": "unknown",
                "safety": "unknown",
                "n_safety_liabilities": 0,
            }
        )
        annotated = pd.concat([annotated, filler], ignore_index=True)
    return annotated


# ==========================================================================
# Composite score
# ==========================================================================


def _minmax(series: pd.Series) -> pd.Series:
    lo, hi = series.min(), series.max()
    if hi == lo:
        return pd.Series(0.0, index=series.index)
    return (series - lo) / (hi - lo)


def compute_composite(
    df: pd.DataFrame, weights: dict[str, float]
) -> tuple[pd.Series, pd.Series, dict[str, float]]:
    """Returns (composite_score, composite_breakdown, used_fractions).

    ``used_fractions`` is the weight of each *available* component
    renormalised to sum to 1 (a component with no data column is dropped and
    its weight redistributed). ``composite_score`` is the sum of the
    per-component **weighted contributions** ``fraction_k * normalised_k``;
    ``composite_breakdown`` is the ``k=contribution`` string of exactly
    those terms, so it adds up to ``composite_score`` and shows what built
    the number, not just the number."""
    components: dict[str, pd.Series] = {}
    if "topology" in weights:
        components["topology"] = _minmax(df["topology_score"])
    if "rwr" in weights and "rwr_score" in df.columns:
        components["rwr"] = _minmax(df["rwr_score"])
    if "genetic_evidence" in weights and "genetic_evidence_score" in df.columns:
        components["genetic_evidence"] = df["genetic_evidence_score"].astype(float).clip(0.0, 1.0)
    if "tractability" in weights:
        components["tractability"] = df["tractability"].map(TRACTABILITY_SCORE).astype(float)
    if "safety" in weights:
        components["safety"] = df["safety"].map(SAFETY_SCORE).astype(float)

    dropped = [k for k in weights if k not in components]
    if dropped:
        logger.warning("Dropping weighted component(s) with no data: %s — weights renormalised.", dropped)

    total = sum(weights[k] for k in components)
    if total <= 0:
        _fail("Every weighted component was unavailable — nothing to score.")
    fractions = {k: weights[k] / total for k in components}

    contributions = {k: fractions[k] * components[k] for k in components}
    composite = sum(contributions.values())
    breakdown = pd.Series(
        [
            "|".join(f"{k}={contributions[k].iloc[i]:.3f}" for k in fractions)
            for i in range(len(df))
        ],
        index=df.index,
    )
    return composite.clip(0.0, 1.0), breakdown, fractions


# ==========================================================================
# Orchestration
# ==========================================================================


def run(args: argparse.Namespace) -> Path:
    weights = parse_weights(args.weights)
    manifest = load_manifest(args.manifest)

    candidates_path = args.candidates or (args.manifest.parent / "candidate_scores.tsv")
    candidates_path = _resolve_beside_manifest(candidates_path, args.manifest, "candidate_scores.tsv")
    graph_dir = args.graph_dir or args.manifest.parent
    target_parquets = resolve_target_parquet(manifest, args.manifest, args.target_parquet)

    df = load_candidates(candidates_path)
    bad = [g for g in df["gene"] if not ENSEMBL_GENE_ID_PATTERN.match(g)]
    if bad:
        _fail(f"candidate_scores.tsv has non-Ensembl gene id(s): {bad[:10]}")
    logger.info("Annotating %d candidates from %s.", len(df), candidates_path)

    evidence = load_genetic_evidence(graph_dir)
    if evidence is not None:
        df["genetic_evidence_score"] = df["gene"].map(evidence).astype(float)
        n_missing = int(df["genetic_evidence_score"].isna().sum())
        if n_missing:
            logger.warning("%d candidate(s) missing from genetic_evidence_score — filled with 0.0.", n_missing)
            df["genetic_evidence_score"] = df["genetic_evidence_score"].fillna(0.0)

    annotations = build_annotations(target_parquets, df["gene"].tolist())
    df = df.merge(annotations, on="gene", how="left")

    n_tract = int((df["tractability"] != "unknown").sum())
    n_safety = int((df["safety"] != "unknown").sum())
    logger.info(
        "OT annotation: tractability %d/%d, safety liabilities %d/%d.",
        n_tract, len(df), n_safety, len(df),
    )

    composite, breakdown, used_fractions = compute_composite(df, weights)
    weights_str = ",".join(f"{k}:{v:.3f}" for k, v in used_fractions.items())
    df["composite_score"] = composite
    df["composite_breakdown"] = breakdown
    df["composite_weights"] = weights_str  # constant per run — makes each row's trace self-contained

    df = df.sort_values("composite_score", ascending=False, kind="stable").reset_index(drop=True)

    ordered = [
        "gene", "topology_score", "rwr_score", "genetic_evidence_score",
        "tractability", "safety", "n_safety_liabilities",
        "composite_score", "composite_breakdown", "composite_weights",
    ]
    df = df[[c for c in ordered if c in df.columns]]
    df = AnnotatedCandidatesSchema.validate(df)

    out_path = args.out or (candidates_path.parent / "candidates_annotated.tsv")
    df.to_csv(out_path, sep="\t", index=False)
    logger.info("Stage 6 complete (weights %s): %s", weights_str, out_path)
    return out_path


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    run(parse_args(argv))


if __name__ == "__main__":
    main()
