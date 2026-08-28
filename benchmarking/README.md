# Benchmark validation

> Ground-truth check for the [Targets from Pathways](../README.md)
> pipeline. Pulled out of the main pipeline into its own module: it *drives*
> the pipeline (Stages 1–6) rather than being a stage of it, it has its own
> human-curated input file, and its results must be read differently from
> the pipeline's own output.

## Purpose

Score a small, literature-curated set of known resistance / compensation
gene pairs (e.g. EGFR→MET) and report where each known `alternative_target`
lands in the pipeline's ranked candidate list — descriptive rank / percentile
recovery against a random (ideally degree-matched) background, **not** a
significance test. Addresses **H8** (the method recovers clinically
established alternative targets better than chance).

This is the check that catches pipeline bugs before any biological
conclusion is trusted. It runs last, on Day 3, and it gates the demo.

## Inputs

- **`benchmarking/resistance_pairs.tsv`** — the curated pair list. Columns:
  `original_target`, `alternative_target`, `disease_or_context`,
  `source_citation` (a real, checked PMID per row). **This file must be
  authored by a human, not generated** — it is the single file the entire
  validation rests on. Curating it is Day 1 literature work done in
  parallel with pipeline coding; if it is missing, the benchmark stops and
  says so rather than inventing cases or citations.
- The full pipeline (Stages 1–6) and its inputs — the benchmark runs the
  pipeline once per known pair.

## Holdout — avoiding self-validation

For each pair, Stage 2 (`pipeline/stage2_gsea_discovery.py`) is run with
`--benchmark-holdout-file` set to exclude the pair's `alternative_target`
from the disease-associated seed genes, so a resistance-pair gene that
happens to overlap Stage 2's own evidence does not validate itself. The
holdout file lists genes one per line (`#` comments allowed) and resolves
them to Ensembl IDs via Stage 1's `genes.parquet` — symbol, synonym, or
bare Ensembl ID. See `example_data/benchmark_holdout_example.txt` for the
format, and the "Benchmark holdout" note in the root README's *Running
Stage 2* section for a worked example.

## Reporting discipline

With n < 20 curated pairs the result is descriptive only — e.g.
`"X of N known pairs ranked in the top 5%"`. **Do not compute or report a
p-value against this benchmark.** The small-n limitation is real and
permanent (the curated set is deliberately 5–15 cases); it is a sanity
check, not a powered evaluation.

## Runtime note

The benchmark re-runs the pipeline many times. GSEA is memoized across
pairs that share a disease so enrichment is not recomputed unnecessarily.
Any change that makes Stage 2 heavier (e.g. a permutation null for
annotation-bias correction — see *Open decisions* in the root README)
multiplies through the benchmark's repeated runs and should be weighed
against that.

## See also

- [`PROMPTS.md`](PROMPTS.md) — the implementation prompt for
  `benchmark_validate.py`.
- [`../README.md`](../README.md) — the pipeline this validates.
