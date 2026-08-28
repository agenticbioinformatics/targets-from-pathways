# Targets from Pathways - User guide

Reference for each stage of the pipeline: setup, run command, expected
output, flags, and the unit tests.

Throughout, `data/` is the source cache (`--data-dir`) and `out/` is the
run directory the stages write into — substitute your own if necessary.

## Full run

One (target, disease) all the way to an HTML report. `run_pipeline.py`
chains Stages 1–6 as subprocesses; `stage7_report.py` renders the HTML report. Every stage script has `--help`.

```
python3 pipeline/run_pipeline.py \
    --target NOD2 --disease MONDO_0005265 \
    --data-dir data/ --out-dir out/ --report
```

Offline (no real Open Targets /
Reactome data), start from the checked-in synthetic Stage 1 output:

```
python3 example_data/build_gsea_example.py

python3 pipeline/run_pipeline.py --from-manifest example_data/example_run/manifest.json --report
```

## Running Stage 0 (`pipeline/stage0_schemas.py`)

Stage 0 is a library, not a script — there is nothing to invoke. It defines
the inter-stage data contracts (the `pandera` `DataFrameSchema` objects and
the `Manifest` pydantic model) that every later stage imports and calls as
`SCHEMA.validate(df)` / `Manifest.model_validate(...)` on both read and
write. All other `pipeline/stageN_*.py` scripts `from stage0_schemas import
...`; keep `pipeline/` importable (run them as `python3 pipeline/stageN_*.py`
from the repo root, which puts `pipeline/` on `sys.path` automatically).

**Setup:** `pip install pandas pandera pydantic numpy` (or `pip install -r requirements.txt`).

**Check it:**

```
python3 -m pytest tests/test_schemas.py
```

`tests/test_schemas.py` asserts each schema *rejects* a deliberately
malformed frame (wrong dtype, missing column, out-of-vocabulary value, and
an HGNC symbol where an Ensembl gene ID is required — the `GRB2-1` bug
class) with an error naming the offending column, as loudly as it accepts a
valid one.

## Running Stage 1 (`pipeline/stage1_ingest.py`)

**Setup:** `pip install pandas pyarrow pandera pydantic` (or `pip install -r requirements.txt`).

**Run command** (the entry point of the whole pipeline — takes a target
gene and a disease ID, nothing upstream):

```
python3 pipeline/stage1_ingest.py --target NOD2 --disease MONDO_0005265 --data-dir data/ --out-dir out/
```

On the first run this downloads the version-pinned source files (Open
Targets Platform 26.06; Reactome release 97 + the `04142025` FI file) into
`--data-dir` and verifies each against a hardcoded sha256, then resolves
every gene identifier to an Ensembl gene ID using the Open Targets `target`
parquet as the authority. Later runs reuse `--data-dir` as a cache. Key
flags:

- `--no-download` — never touch the network; `--data-dir` must already be
  populated (see the error message it prints for the exact mirror layout).
- `--allow-low-coverage` — do not fail when identifier-mapping coverage is
  below 90% or the disease's `genetic_association` gene count is
  implausibly small. Off by default, so an unusable run aborts early.
- `--fi-curated-only` (default on) — keep only curated Reactome functional-
  interaction rows; predicted rows carry no reliable direction/sign.
- `--min-set-size` / `--max-set-size` — gene-set size cap applied before
  Stage 2.

**Expected output** — seven artifacts in `--out-dir`, each schema-validated
on write: `genes.parquet` (the canonical Ensembl-keyed gene table),
`gene_sets.parquet` and `interactions.parquet` (the two database-agnostic,
`source_db`-tagged Reactome tables), `ot_disease_subset.parquet` (OT
associations subset to the disease), `coverage_report.json` (per-source
identifier-mapping rate and example dropped IDs), `scale_report.json`
(gene-set size distribution and projected Stage 3 edge count), and
`manifest.json`. Every later stage reads only these, resolved through
`manifest.json` — never a raw source file.

There is no tiny checked-in example run for this stage: `stage1_ingest.py`
needs the real pinned sources (or the test fixtures). The
`example_data/example_run/` directory is synthesised directly by
`python3 example_data/build_gsea_example.py`, which writes a hand-built,
schema-valid Stage 1 output *without* running ingest, so Stages 2–4 have
something to consume offline. The ingest logic itself — pinned-download
verification, PE→gene mapping, coverage reporting, manifest determinism —
is covered by:

```
python3 -m pytest tests/test_ingest_identifier_mapping.py tests/test_ingest_manifest_determinism.py tests/test_ingest_schema_conformance.py
```

which patch the acquisition layer to serve the tiny fixtures under
`tests/fixtures/` (see `tests/conftest.py`) so nothing touches the network.

## Running Stage 2 (`pipeline/stage2_gsea_discovery.py`)

**Setup:** `pip install blitzgsea statsmodels pandas pandera pydantic pyarrow` (or `pip install -r requirements.txt`).

**Run command**, against Stage 1's output for a real (target, disease) pair:

```
python3 pipeline/stage2_gsea_discovery.py --manifest out/manifest.json
```

`pipeline/stage2_gsea_discovery.py` takes no other input — it resolves `gene_sets.parquet`,
`ot_disease_subset.parquet`, and the target gene entirely from the manifest,
and writes `disease_pathways.tsv` next to it.

A tiny hand-built example (18 genes, 5 gene sets including one deliberate
ancestor/child near-duplicate pair and one scattered negative control) is
checked in under `example_data/` — regenerate it with
`python3 example_data/build_gsea_example.py`, then:

```
python3 pipeline/stage2_gsea_discovery.py --manifest example_data/example_run/manifest.json
```

**Expected output** (`example_data/example_run/disease_pathways.tsv`,
columns per `schemas.DiseasePathwaysSchema`; exact `pval`/`fdr` floats jitter
slightly run-to-run — blitzgsea's internal calibration isn't bit-for-bit
deterministic even with a fixed seed — but the three pathways returned, their
order, and `contains_target` are stable):

```
set_id     gene_set                  source_db  pval       fdr        contains_target
R-HSA-101  Child Broad Signaling     reactome   ~7e-05     ~3e-04     False
R-HSA-300  Unrelated Weak Pathway    reactome   ~4e-04     ~9e-04     False
R-HSA-200  Target Signaling Pathway  reactome   ~8e-04     ~1e-03     True
```

The log line `GSEA: 18 seed genes, 4/5 gene sets tested, 3 significant`
confirms the mechanics that matter: `R-HSA-100` ("Ancestor Broad
Signaling", Jaccard 0.83 against its child `R-HSA-101`) is the one set
dropped by collapsing before testing (`Collapsed 1 near-duplicate gene
set(s)`), so 4 of 5 are tested; `R-HSA-400` ("Scattered Noise Pathway", no
coherent enrichment signal) is tested but excluded by the FDR filter, so 3
of 4 are significant.

**Benchmark holdout** (`--benchmark-holdout-file`, for benchmark validation
— see [`benchmarking/`](../benchmarking/README.md)): excludes genes from the
disease-associated seed set before GSEA runs, so a
literature-curated resistance-pair gene that happens to overlap Stage 2's
own evidence doesn't validate itself. Genes are listed one per line
(`#` comments allowed) and resolved to Ensembl IDs via Stage 1's
`genes.parquet` — symbol, synonym, or bare Ensembl ID, the same resolution
`pipeline/stage1_ingest.py` uses for `--target` — never matched directly against
`ot_disease_subset` (which carries no symbol column at all). Example:

```
python3 pipeline/stage2_gsea_discovery.py --manifest example_data/example_run/manifest.json \
    --benchmark-holdout-file example_data/benchmark_holdout_example.txt
```

`example_data/benchmark_holdout_example.txt` holds `OLDG17`, a *synonym*
(not the approved symbol) for the gene at the core of `R-HSA-300`. Excluding
it drops `R-HSA-300` to 4 genes — one below blitzgsea's minimum tested-set
size — so it disappears from `disease_pathways.tsv` entirely (`GSEA tested 3
pathways`, down from 4) rather than merely losing significance.

`python3 -m pytest tests/test_gsea_discovery.py` runs the unit-level GSEA
tests (`tests/test_gsea_discovery.py`, 8 tests): a deliberately-enriched
pathway is returned significant, an individually-significant pathway is
correctly excluded once BH-FDR is applied across the tested library, non-
enriched pathways are tested but excluded from output, and a near-duplicate
ancestor/child pair collapses to the smaller child before GSEA ever sees the
ancestor.

## Running Stage 3 (`pipeline/stage3_build_graph.py`)

**Setup:** `pip install networkx scipy pandas pandera pydantic pyarrow` (or `pip install -r requirements.txt`).

**Run command**, after Stage 1 and Stage 2 have produced `manifest.json`/`disease_pathways.tsv`:

```
python3 pipeline/stage3_build_graph.py --manifest out/manifest.json
```

Reads `gene_sets.parquet` and `interactions.parquet` via the manifest, and
`disease_pathways.tsv` alongside it by default (override with
`--disease-pathways`). Writes `graph_weight.npz`, `graph_sign.npz`,
`graph_gene_index.json`, and `graph_metadata.json` next to the manifest —
not GraphML; see the module docstring for why.

Against the checked-in example (regenerate with
`python3 example_data/build_gsea_example.py` if needed):

```
python3 pipeline/stage3_build_graph.py --manifest example_data/example_run/manifest.json
```

**Expected output** (log line):

```
INFO build_graph: Graph: 17 nodes, 136 edges (4 pathways in union, 2 interaction edges in scope, 1 dropped out of scope).
```

The one dropped interaction row (`G12 -> TARGET`) is deliberate: `G12`
belongs only to `R-HSA-100`, the redundant ancestor pathway Stage 2's
collapsing excluded from the union — proof that interaction edges really
are scoped to the pathway-gene pool, not added genome-wide.

`python3 -m pytest tests/test_build_graph.py` runs the unit-level graph
tests (4 tests): the exact directed/signed edge list (with `source_db`
tags) for a hand-built pathway + interaction fixture, the co-membership
weight formula for both a small and a large (overlapping) pathway, and that
a gene alone in its own pathway survives as an isolated node.

## Running Stage 4 (`pipeline/stage4_genetic_evidence_weights.py`)

**Setup:** same as Stage 3, plus `numpy` (already required transitively).

**Run command**, after Stage 3 has produced the graph artifacts:

```
python3 pipeline/stage4_genetic_evidence_weights.py --manifest out/manifest.json
```

Reads `ot_disease_subset.parquet` via the manifest and Stage 3's graph
artifacts alongside it by default (override with `--graph-dir`/`--out-dir`).
`--datatypes` (default `genetic_association,known_drug`) is an *inclusion*
list — `affected_pathway` and any literature/pathway-derived datatype are
excluded simply by not being named, which matters concretely: OT 26.06 maps
its own `reactome` evidence datasource to `affected_pathway` (see
`pipeline/stage1_ingest.py`'s `_DATASOURCE_TO_DATATYPE`), so including it would let
Reactome-derived evidence weight a graph whose topology already *is*
Reactome pathway structure.

Writes **exactly Stage 3's four filenames**, not a parallel set, so
`stage5_score_candidates.py` (Stage 5) can point `--graph` at either a Stage 3 or
a Stage 4 directory and load it through the same code path:
`graph_weight.npz` is *replaced* with the genetic-evidence-derived edge
weight (same sparsity pattern as Stage 3's structural weight — no edge
added or removed, only re-weighted); `graph_sign.npz`/`graph_gene_index.json`
are copied through unchanged; `graph_metadata.json` is extended (not
replaced) with a `genetic_evidence_score` gene->score mapping. Plus a
separate, human-facing `gene_weights.tsv` that Stage 5 never reads. See the
attribute-name comment at the top of `pipeline/stage4_genetic_evidence_weights.py` for the
exact contract.

Against the checked-in example:

```
python3 pipeline/stage4_genetic_evidence_weights.py --manifest example_data/example_run/manifest.json
```

**Expected output** (`example_data/example_run/gene_weights.tsv`, top rows):

```
gene_id           genetic_evidence_score
ENSG00000000002   0.98
ENSG00000000003   0.95
ENSG00000000004   0.92
```

`--edge-weight-mode product` (default `avg`) changes only the values written
into `graph_weight.npz`; e.g. the target-`ENSG...002` edge is
`avg(0.60, 0.98) = 0.79` by default or `0.60 * 0.98 = 0.588` with `product`.

## Running Stage 5 (`pipeline/stage5_score_candidates.py`)

**Setup:** `pip install networkx scipy numpy pandas pandera pydantic` (or `pip install -r requirements.txt`).

**Run command**, against a Stage 3 *or* Stage 4 graph directory (the two
share one format — see Stage 4 above):

```
python3 pipeline/stage5_score_candidates.py --graph-dir out/
```

`--graph-dir` is the only required argument: it reads `graph_weight.npz`,
`graph_sign.npz`, `graph_gene_index.json` and `graph_metadata.json` from
there, takes the target node from the metadata's `resolved_target.gene_id`
(override with `--target`), and writes `candidate_scores.tsv` alongside
(override with `--out`).

`--method` is a comma-separated subset of `{topology, rwr}`, **default
`topology`**. `topology_score` is always emitted; `rwr_score` appears only
when `rwr` is requested (`--method topology,rwr` or `--method rwr`). RWR
knobs: `--restart-prob` (default 0.5; PageRank `alpha = 1 - restart_prob`)
and `--node-weight-mix` (default 0.5 — on a Stage 4 graph, extra restart
mass spread over genes by `genetic_evidence_score`; the target always keeps
the largest share, so 0.5 leaves it ~2/3 of the restart mass; ignored on a
Stage 3 graph).

Against the checked-in example (after Stages 1–4 have run, so
`example_data/example_run/` holds a **Stage 4** graph):

```
python3 pipeline/stage5_score_candidates.py --graph-dir example_data/example_run --method topology,rwr
```

**Expected output** (`example_data/example_run/candidate_scores.tsv`, top
rows; columns per `schemas.CandidateScoresSchema`, ranked by `rwr_score`
since `rwr` was requested — else by `topology_score`):

```
gene             topology_score  rwr_score
ENSG00000000004  3.4545          0.0769
ENSG00000000002  2.6480          0.0717
ENSG00000000003  2.6197          0.0703
ENSG00000000009  3.2522          0.0688
ENSG00000000005  2.5616          0.0675
```

The five weak "unrelated pathway" genes (`ENSG...017`–`021`) land at the
bottom with scores near zero — they reach the target only by long paths and
share none of its pathway neighbourhood, which is the negative-control
behaviour the example is built to show. On the **Stage 3** graph (run
`stage3_build_graph.py` but not `stage4_...` into the directory) the default
`topology`-only ranking instead leads with `ENSG...002`, the target's
direct interaction partner (`TARGET -> G2`, confidence 0.9).

The `gene` column is an **Ensembl gene ID** — the project's one canonical
identifier (Decision 2), copied verbatim from `graph_gene_index.json`, not
mapped. It is the same namespace as `genes.gene_id`,
`ot_disease_subset.gene_id` and Stage 4's `genetic_evidence_score` keys, so
Stage 6 joins `candidate_scores.gene` directly onto Open Targets
tractability/safety with no ID-mapping step; `CandidateScoresSchema`'s
`^ENSG\d{11}$` check (plus an explicit guard in `run()`) makes a namespace
slip fail loudly rather than silently drop rows on the join. HGNC symbols
are a display-time lookup (`genes.parquet`) in Stages 6/7, never written
here.

`python3 -m pytest tests/test_score_candidates.py` runs the unit-level
tests (10 tests): hand-computed `topology_score` values on a small
cycle+chain graph (both at the function level and end-to-end through a
default `run()` with no `--method`), an isolated node scoring 0, RWR
reconstructed as a full stationary distribution (non-negative, sums to 1)
that favours the target's neighbours over far-away chain nodes, the
node-weight blend pulling mass toward a high-evidence gene, and the
output-column contract (`topology_score` always, `rwr_score` only when
requested, target row excluded, unknown `--method` exits non-zero).

## Running Stage 6 (`pipeline/stage6_annotate_context.py`)

**Setup:** `pip install pandas pyarrow pandera pydantic` (or `pip install -r requirements.txt`).

**Run command**, after Stage 5 has produced `candidate_scores.tsv`:

```
python3 pipeline/stage6_annotate_context.py --manifest out/manifest.json
```

Resolves, all next to `--manifest` by default: `candidate_scores.tsv`
(Stage 5), `graph_metadata.json` (Stage 4's `genetic_evidence_score` map —
optional), and the Open Targets `target` parquet (from
`manifest.sources` — the `opentargets` `*target*.parquet` files; override
with `--target-parquet`). Writes `candidates_annotated.tsv` alongside the
candidates file.

Per candidate it adds:

- **`tractability`** — `clinical` (an Open Targets `value==true` flag for a
  clinical stage), `discovery` (some other true flag), or `unknown` (no
  true flag, or the gene is absent from the `target` parquet).
- **`safety`** — `has_liabilities` (one or more `safetyLiabilities`
  entries) or `unknown`. **There is no `safe` value**: Open Targets records
  known liabilities or says nothing, and its safety coverage is sparse, so
  a gene with no annotation is `unknown` (a neutral 0.5 in the composite),
  never rewarded as if confirmed safe. `n_safety_liabilities` carries the
  count.
- **`composite_score`** ∈ [0, 1] — `--weights` (default
  `topology:0.3,rwr:0.3,genetic_evidence:0.2,tractability:0.1,safety:0.1`)
  weighted-averaged over: min-max-scaled `topology_score` and `rwr_score`
  (proximity numbers only meaningful ranked within the candidate set),
  as-is `genetic_evidence_score`, and the tractability
  (`clinical`=1/`discovery`=0.5/`unknown`=0) and safety
  (`has_liabilities`=0/`unknown`=0.5) ordinals. Weights need not sum to 1;
  a component with no data (`rwr` on a topology-only Stage 5 run,
  `genetic_evidence` with no Stage 4) is dropped and the rest
  renormalised.
- **`composite_breakdown`** — `k=contribution|...`, each component's
  *weighted contribution* to `composite_score` (the terms sum to it), so
  the trace shows what built the rank, not just the number.
- **`composite_weights`** — the renormalised weights actually used
  (`k:fraction,...`, sum 1), repeated on every row so one row is a
  self-contained explanation for the Stage 7 report.

Per-candidate **evidence trace**: the passthrough score columns
(`topology_score` / `rwr_score` / `genetic_evidence_score`), the annotation
buckets (`tractability` / `safety` / `n_safety_liabilities`), and
`composite_breakdown` + `composite_weights` together explain every rank.
Two provenance items are deliberately **not** here — *which shared
pathways/interactions* (not in this stage's inputs; Stage 7 assembles it
from `gene_sets.parquet` at drill-down time) and *which Open Targets
datatype* produced the genetic-evidence score (a Stage 4 follow-up — Stage
4 records only the final score).

**No tissue-expression data of any kind is read or used** — only `id`,
`tractability` and `safetyLiabilities` are pulled from the `target`
parquet, by design (plan update 1, H10 dropped). No rows are filtered out.

Against the checked-in example (after Stages 1–5):

```
python3 pipeline/stage6_annotate_context.py --manifest example_data/example_run/manifest.json
```

**Expected output** (`example_data/example_run/candidates_annotated.tsv`,
top rows; columns per `schemas.AnnotatedCandidatesSchema`, ranked by
`composite_score`):

```
gene             topology_score  genetic_evidence_score  tractability  safety           n_safety_liabilities  composite_score
ENSG00000000004  3.4545          0.92                    discovery     has_liabilities  1                     0.834
ENSG00000000002  2.6480          0.98                    clinical      has_liabilities  2                     0.805
ENSG00000000003  2.6197          0.95                    discovery     unknown          0                     0.790
ENSG00000000009  3.2522          0.77                    unknown       unknown          0                     0.754
ENSG00000000005  2.5616          0.89                    unknown       unknown          0                     0.712
```

`ENSG...004` stays on top despite a safety liability because its proximity
and genetic-evidence terms dominate; `ENSG...017`, clinically tractable but
far from the target and weakly associated, lands near the bottom
(`composite_score` ~0.16) — tractability alone can't lift a low-proximity
candidate. The example's synthetic `ot_target_subset.parquet` covers 6 of
the 17 genes; the other 11 are `unknown`/`unknown` and ranked purely on the
score terms, never penalised for missing annotation.

`python3 -m pytest tests/test_annotate_context.py` runs the unit-level
tests (9 tests): tractability bucketing and the safety flag/count, weight
parsing and rejection of junk, a hand-computed composite for a known
two-row input (with `composite_breakdown` contributions that sum back to
the score), a second hand-computed composite driven end-to-end from a
small fixture tractability/safety table, component-drop-and-renormalise
when `rwr`/`genetic_evidence` data is absent, and — the case the spec
calls out — a gene with an empty `safetyLiabilities` list and a gene absent
from the `target` parquet both coming out `unknown` rather than defaulting
to a numeric or `safe` value, in `build_annotations` and end-to-end.

## Running Stage 7 (`pipeline/run_pipeline.py` + `pipeline/stage7_report.py`)

**Setup:** `pip install pandas pyarrow pandera pydantic` (or `pip install -r requirements.txt`).
No web framework — `stage7_report.py` writes one self-contained HTML file
using only the standard library for templating.

**Orchestrator** — `run_pipeline.py` chains Stages 1–6 for one
(target, disease):

```
python3 pipeline/run_pipeline.py --target NOD2 --disease MONDO_0005265 \
    --data-dir data/ --out-dir out/ --report
```

Each stage runs as its own subprocess (the same CLI as its "Running Stage
N" section above). If one exits non-zero, `run_pipeline.py` prints which
stage failed, its exit code and the exact command, and stops the run
non-zero — the stage's own message (e.g. Stage 1's "Could not resolve
`--target 'FOO'`. Did you mean: …?") is right above it, not buried under a
traceback. `--report` additionally writes `out/report.html` at the end.

Offline, against the checked-in synthetic Stage 1 output (no real Open
Targets / Reactome data needed) — this regenerates Stages 2–6 in place:

```
python3 pipeline/run_pipeline.py --from-manifest example_data/example_run/manifest.json --report
```

**Report only** — `stage7_report.py` turns a finished run directory into
the HTML report without re-running anything:

```
python3 pipeline/stage7_report.py --run-dir out/
```

**What the report contains** (single file, opens in any browser):

- header — target (HGNC symbol + Ensembl), disease, seed, Reactome / Open
  Targets versions, whether the graph is Stage 3 (structural) or Stage 4
  (genetic-evidence-weighted);
- a **sortable** ranked table — `symbol`, `gene` (Ensembl), `composite`,
  `topology`, `rwr`, `genetic ev.`, `tractability`, `safety` (click a
  header to sort; `rwr`/`genetic ev.` columns appear only if those stages
  ran);
- a **click-to-expand evidence trace** per candidate — its
  `composite_breakdown` + weights, the pathways it shares with the target
  (and the disease-relevant pathways it is merely a member of), its Open
  Targets datatypes + scores (with the ones feeding `genetic_evidence_score`
  flagged), and the signed interaction edges touching it (interactions to
  non-graph genes are omitted, matching Stage 3's scoping);
- a **static** "Benchmark validation" panel — the text of
  `benchmarking/benchmark_summary.txt`, or, until that exists,
  `benchmarking/benchmark_summary.example.txt` shown with a clear
  "EXAMPLE — not a real validation run" label. Never recomputed here; see
  [`benchmarking/`](../benchmarking/README.md).

Genes are displayed by HGNC symbol; the Ensembl gene ID stays the internal
key (its own column, and the anchor id for each detail block). Symbols are
resolved from `genes.parquet` at render time only.

Only `manifest.json` and `candidates_annotated.tsv` are required; every
other artifact is optional. A partial run directory still produces a
report — missing trace inputs (`gene_sets.parquet`, `interactions.parquet`,
`ot_disease_subset.parquet`, …) show that section as unavailable, a
Stage-3-only run (no Stage 4) just drops the `rwr`/`genetic ev.` columns,
no `genes.parquet` falls the symbol column back to Ensembl ids, and no
benchmark summary shows a "not run yet" note — a header line lists whatever
was missing.

**Expected output** (`example_data/example_run/report.html`, top of the
ranked table — same ordering as Stage 6's `candidates_annotated.tsv`):

```
#  symbol  gene             composite  topology  rwr     genetic ev.  tractability  safety
1  GENE04  ENSG00000000004  0.834      3.455     0.0769  0.920        discovery     has_liabilities
2  GENE02  ENSG00000000002  0.805      2.648     0.0717  0.980        clinical      has_liabilities
3  GENE03  ENSG00000000003  0.790      2.620     0.0703  0.950        discovery     unknown
```

Expanding `GENE04` shows: shared pathways *Target Signaling Pathway* and
*Scattered Noise Pathway*; also a member of *Child Broad Signaling*;
`genetic_association` = 0.92 (feeds the genetic-evidence score); interaction
`GENE04 → GENE09` (inhibiting, curated, confidence 0.60).

`python3 -m pytest tests/test_stage7_report.py tests/test_run_pipeline.py`
runs the tests (11): `evidence_trace`'s shared-pathway / graph-scoping /
datatype logic on a hand-built run, the `feeds` flag going false when Stage
4 did not run, the `benchmark_summary` real→example→none fallback,
`build_report` producing a self-contained symbol-labelled HTML with one
detail block per candidate and the benchmark panel, a partial run
directory (only `manifest.json` + `candidates_annotated.tsv`) still
rendering with the optional columns/trace sections marked unavailable, a
`run_pipeline --from-manifest --report` smoke test that Stages 2–6 chain
and produce every expected artifact, and a stage-failure producing a
helpful "which stage / which command / exit code" message rather than a
traceback.
