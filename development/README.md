# How this pipeline was built

This directory is the paper trail for **Targets from Pathways**.

| file | what it is |
|---|---|
| [`PROMPTS.md`](PROMPTS.md) | The implementation prompts, one per pipeline stage, matched to the Day 1/2/3 hackathon plan. Each is self-contained: hand it to a fresh agent session and it builds that stage. |
| [`prompt_history`](prompt_history) | Running plain-text log of every request made during development and a summary of what changed in response. |
| [`hackathon-pipeline-dev/SKILL.md`](hackathon-pipeline-dev/SKILL.md) | A transferable agentic skill — generic to any bioinformatics pipeline built from a problem statement plus a stage-by-stage description. |

## The process, in short

1. **Plan first.** A senior-bioinformatician brainstorm scored ~10 research
   hypotheses on rigour / data availability / 3-day feasibility, merged the
   survivors into one modular pipeline, and put it through several rounds of
   professorial review. That produced a "Pipeline overview" + "Plan
   updates" + "Open decisions" document, written before any stage code —
   since trimmed out of the root README into the project's git history, and
   carried forward as the per-stage prompts in [`PROMPTS.md`](PROMPTS.md).

2. **Stage 0 is the contract, and it's enforced.** `stage0_schemas.py`
   defines every inter-stage table (`pandera`) and the manifest
   (`pydantic`). Every stage calls `SCHEMA.validate(df)` on **read and
   write** — prose comments describing columns were explicitly rejected as
   insufficient. `tests/test_end_to_end.py` re-validates every hand-off, so
   a renamed column anywhere fails one test.

3. **One script per stage, chained by output path.** `stageN_*.py` each
   take the previous stage's output directory (or its `manifest.json`) as a
   CLI argument and are independently runnable. `run_pipeline.py`
   orchestrates them as subprocesses; `stage7_report.py` renders the run.

4. **A synthetic example that bypasses Stage 1.** Stage 1 needs multi-GB
   pinned Open Targets / Reactome downloads, so `example_data/build_gsea_example.py`
   hand-writes a tiny, schema-valid Stage 1 output (`example_data/example_run/`)
   that Stages 2–7 consume offline. Stage 1's own code path is covered by
   `tests/test_ingest_*.py` against `tests/fixtures/`.

5. **Test as you go.** Each stage landed with `tests/test_<stage>.py`
   (hand-computed expected values on tiny fixtures), then the end-to-end
   schema-conformance test. ~100 tests, no network.

6. **The plan changed six times.** Multi-database support was added then
   withdrawn (Reactome only); the backtrack-free-walk scoring method and
   its C-extension dependency were added then removed; benchmark validation
   was promoted into the pipeline, then split back out into
   [`../benchmarking/README.md`](../benchmarking/README.md) as its own
   driver; the web app became a static HTML report. Each reversal was
   recorded as a numbered "Plan update" (in the README's history) rather
   than edited away, and is summarised in [`prompt_history`](prompt_history).

7. **Renames happened late and in bulk.** `schemas.py` → `stage0_schemas.py`
   etc. once the stage numbering settled; `example_data/stage1_run/` →
   `example_data/example_run/`; Stage 8 → Stage 7 after the benchmark split.
   `git mv` throughout, references swept with the change.

See `prompt_history` for the blow-by-blow.
