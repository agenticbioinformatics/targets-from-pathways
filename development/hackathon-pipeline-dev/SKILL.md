---
name: hackathon-pipeline-dev
description: >-
  Turn a problem statement plus a stage-by-stage pipeline description into a
  modular, reproducible, tested bioinformatics pipeline under time pressure
  — schema-first contracts, one independently-runnable script per stage
  chained by output path, synthetic example data, per-stage plus
  end-to-end tests, and a single static-HTML deliverable. Works for any
  input data and any databases (sequence, variant, expression, network,
  annotation, imaging). Use when you have a clear pipeline spec and need to
  build it. Triggers on: "build this pipeline", "chain these analysis
  steps", "turn this notebook/spec into a real pipeline", "ingest <public
  DB> and produce <ranked output / model / report>".
---

# Building a bioinformatics pipeline from a spec

You have two things: a **problem statement** and a **pipeline description**
(an ordered list of stages, each with an input and an output). Turn them
into working, reproducible, tested code. Apply these patterns in order —
each is something that costs rework when skipped.

## 1. Lock the plan in a README before writing stage code

Write a short design doc first: the question being answered, every input
data source **with its exact version**, the stage list (one line each: what
it consumes, what it emits), the hypotheses or claims being tested, and the
**open decisions you are deliberately not resolving yet** (method choices,
thresholds, cutoffs). If several approaches are plausible, score them on
scientific rigour × data availability × time budget and pick one. When the
plan changes later — it will — append a dated/numbered "Plan update" rather
than editing history; the reversal reasoning is the most reusable artifact.

## 2. Stage 0 is a schema module, validated on read *and* write

Before stage 1, define the inter-stage contracts in one module:

- one table schema per artifact that crosses a stage boundary
  (`pandera`/`pandera`-like), `strict` (no unexpected columns);
- one schema for the run manifest / config (`pydantic`-like);
- a foreign-key / referential-integrity helper.

Every stage validates each input right after reading it and each output
right before writing it. Prose comments describing columns drift and fail
silently. Carry a single **canonical identifier** end to end (gene/variant/
sample/feature ID), pattern-checked, and a **provenance column** naming the
source, even if only one source is used now.

## 3. One script per stage, chained by output path

`stageN_<name>` scripts, each:

- takes the previous stage's **output directory or manifest path** as a CLI
  argument — no shared in-memory state, each independently runnable;
- writes into that same run directory so the next stage finds it;
- has an `argparse` (or equivalent) parser with a `description` and
  per-argument `help`, so `--help` *is* the documentation;
- validates inputs (file exists, expected columns, non-empty) and fails
  with a specific one-line error + non-zero exit, never a bare traceback;
- guards the obvious edge cases you meet while reading the data
  (empty input, divide-by-zero, all-NaN column, single-row group).

Add a thin `run_pipeline` orchestrator that runs the stages as
**subprocesses** — a crash in one prints which stage failed and its exact
command, and stops the chain. Give it a flag to start from a later stage's
input so the expensive first stage can be skipped on re-runs.

## 4. Synthetic example data that bypasses the un-runnable stage

Some stage — usually ingestion — needs large downloads, credentials, or a
cluster, and can't run in a quick test. So:

- write a generator that hand-builds a **tiny, schema-valid** copy of that
  stage's *output* (a handful of records, sized to exercise every
  downstream branch, including at least one deliberate negative control and
  one edge case). Commit it. Downstream stages then run with no network.
- Give the un-runnable stage its own tests that monkeypatch its
  acquisition/IO layer to serve tiny committed fixtures.
- Match the real formats exactly (same columns, dtypes, delimiters, file
  layout) so swapping in real data is a no-op.
- Keep synthetic values **biologically plausible** for what the stage
  measures — a method that "works" on data that couldn't occur doesn't
  validate anything.

## 5. Classes of gotcha to check for every stage

- **Statistical methods have minimum input sizes.** Enrichment tests,
  permutation nulls, distribution fits, clustering, dimensionality
  reduction — each silently returns garbage or crashes below some N of
  samples/features. Size synthetic fixtures above that floor, or the test
  for that stage can't run.
- **Circularity.** If you weight/score using evidence derived from source
  X and the thing being scored *is* source X, you are validating X against
  itself. Restrict to an explicit non-circular subset and comment why.
- **Version drift.** Pin exact release numbers of every database and tool;
  never a `latest`/`current` pointer. Record a checksum. Verify each URL
  resolves before hardcoding it.
- **Identifier hygiene.** Map to one canonical ID space at a single point,
  with a format check that rejects near-misses (wrong ID type, versioned
  accession, physical-entity vs gene). Prefer the mapping file at the
  granularity you actually need.
- **Coordinate / orientation correctness.** 0- vs 1-based, half-open vs
  closed, strand, reference build — never shortcut these; an off-by-one
  here is invisible until the biology is wrong.
- **Determinism.** Fix and record seeds; sort before writing; make repeat
  runs byte-identical modulo timestamps so a diff means a real change.

## 6. Test each stage as it lands, then end-to-end

Per stage: a test with **hand-computed expected values** on a tiny fixture
(the exact edge list, the exact score, the exact rank) — not just "runs
without error". Then one end-to-end test that executes every real stage
script in order and re-validates every hand-off against its schema (no
missing/renamed columns, canonical IDs and provenance preserved, the final
output populated for every input). That end-to-end test is what catches a
schema break on demo morning. Keep the whole suite network-free and fast
(seconds); give every script a `--test`/example mode that exercises it on a
minimal slice.

## 7. Deliver a single self-contained HTML report

For the final artifact, prefer a script that reads the finished run
directory and writes ONE self-contained `.html` (inline CSS/JS, no server,
no external assets) over a web app — zero install, shareable as a file,
opens from a `file://` URL. Make it degrade gracefully: only the core
result table plus the manifest are required; a missing optional input
narrows the report with a visible note rather than crashing. Show entities
by human-readable label, keep the stable ID as the key and the anchor.

## 8. Reorganize once, in bulk, at the end

Do renames and directory moves after the structure has settled, not
piecemeal during development. `git mv` (never delete-and-recreate), sweep
every reference in the same change, and keep the docs split: a root README
with a **Quick start** only, a `docs/` with the detailed per-stage
reference, and a `development/` with the plan, the prompt/decision log, and
this skill.
