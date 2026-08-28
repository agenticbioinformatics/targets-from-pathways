# Benchmark validation — implementation prompt

Ready-to-use prompt for an agentic coding assistant, matched to the Day 3
plan in [`../README.md`](../README.md). Self-contained — hand it to a fresh
agent session with no other context beyond this repo. See
[`README.md`](README.md) for what this module is and why it lives outside
the main pipeline.

---

## Benchmark validation (`benchmarking/benchmark_validate.py`)

**Implement:**
> Write `benchmarking/benchmark_validate.py` that takes a small hand-curated
> TSV of known resistance/compensation gene pairs (columns:
> `original_target`, `alternative_target`, `disease_or_context`,
> `source_citation`) — this file **must be curated by a human, not generated
> here** — it is the single file the entire validation rests on, and every
> `source_citation` needs a real, checked PMID. Assume it already exists at
> `benchmarking/resistance_pairs.tsv`; if it does not, stop and say so
> rather than inventing cases or citations. Run the full pipeline (Stages
> 1–6) for each known pair with `--benchmark-holdout-file` set to exclude
> the `alternative_target` gene from Stage 2's seed genes, memoizing the
> disease GSEA across pairs that share a disease so the benchmark does not
> re-run enrichment unnecessarily. Record the rank and percentile of
> `alternative_target` in the resulting candidate list. Print a clear
> disclaimer that with n < 20 cases results must be reported descriptively
> (e.g. "X of N known pairs ranked in top 5%"), not as a significance test —
> do not compute or report a p-value against this benchmark.

**Test:**
> Add a pytest test using a tiny synthetic pipeline setup (small graph,
> small association table) with one deliberately "easy" known pair (short
> graph distance) and confirm `benchmark_validate.py` reports it near the
> top of the ranking with the correct percentile calculation.

**Wire to the report:**
> Ensure `benchmark_validate.py` writes its summary as both a TSV and a
> short human-readable text block, so the Stage 7 report (output/report +
> web UI) can embed the validation summary directly without re-parsing raw
> output.
