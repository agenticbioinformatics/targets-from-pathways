# Targets from Pathways - User guide

Reference for each stage of the pipeline: setup, run command, expected
output, flags, and the unit tests.

Throughout, `data/` is the source cache (`--data-dir`) and `out/` is the
run directory the stages write into — substitute your own if necessary.

## Full run

One disease all the way to an HTML report. `run_pipeline.py`
chains Stages 1–6 as subprocesses; `stage7_report.py` renders the HTML report. Every stage script has `--help`.

```
python3 pipeline/run_pipeline.py \
    --disease MONDO_0005265 \
    --data-dir data/ --out-dir out/ --report
```

There is no target gene: the disease's genes with Open Targets
`genetic_association` evidence are the **seed genes** (Stage 2 writes them to
`disease_seed_genes.tsv`), the graph is built around them (Stage 3), and
Stage 5 propagates from them to rank every non-seed graph gene.

Any option `run_pipeline.py` does not recognise is forwarded to
`stage1_ingest.py`, so a run that needs one of Stage 1's flags is just
`… --report --allow-low-coverage` (see `stage1_ingest.py --help`).

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

`python3 -m pytest tests/test_run_directory.py` (6 tests) checks run-directory
integrity on copies of the example run: a copied run reads its own Stage 1
files even when the manifest lists the original's, an artifact changed after
Stage 1 is refused, a benchmark holdout gene's own evidence never feeds Stage
4 (and a re-run without a holdout clears it), and re-running Stage 4 gives
the same graph.

## Running Stage 1 (`pipeline/stage1_ingest.py`)

**Setup:** `pip install pandas pyarrow pandera pydantic` (or `pip install -r requirements.txt`).

**Run command** (the entry point of the whole pipeline — takes a disease
ID, nothing upstream):

```
python3 pipeline/stage1_ingest.py --disease MONDO_0005265 --data-dir data/ --out-dir out/
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
`manifest.json` — never a raw source file. A stage looks for each file next
to the manifest first (the listed path is relative to wherever Stage 1 ran,
so in a copied or moved run directory it could point at another run), and
refuses a file whose sha256 no longer matches the manifest.

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

**Run command**, against Stage 1's output for a real disease:

```
python3 pipeline/stage2_gsea_discovery.py --manifest out/manifest.json
```

`pipeline/stage2_gsea_discovery.py` takes no other input — it resolves `gene_sets.parquet`
and `ot_disease_subset.parquet` entirely from the manifest, and writes three
files next to it: `disease_pathways.tsv` (the significant pathways),
`disease_seed_genes.tsv` (`gene_id`, `score`: the genes with
`genetic_association` evidence — the GSEA ranking signature — minus any
benchmark holdout) and `benchmark_holdout_genes.tsv` (`gene_id`: the
resolved holdout, header only without `--benchmark-holdout-file`). Stages 3
and 5 use `disease_seed_genes.tsv` as the seed genes; Stage 4 reads the
holdout.

A tiny hand-built example (30 genes of which 18 are seeds, 7 gene sets
including one deliberate ancestor/child near-duplicate pair, one scattered
negative control and one pathway with no seed gene) is
checked in under `example_data/` — regenerate it with
`python3 example_data/build_gsea_example.py`, then:

```
python3 pipeline/stage2_gsea_discovery.py --manifest example_data/example_run/manifest.json
```

**Expected output** (`example_data/example_run/disease_pathways.tsv`,
columns per `schemas.DiseasePathwaysSchema`; exact `pval`/`fdr` floats jitter
slightly run-to-run — blitzgsea's internal calibration isn't bit-for-bit
deterministic even with a fixed seed — but the three pathways returned and
their order are stable):

```
set_id     gene_set                  source_db  pval       fdr
R-HSA-101  Child Broad Signaling     reactome   ~8e-05     ~3e-04
R-HSA-300  Unrelated Weak Pathway    reactome   ~4e-04     ~9e-04
R-HSA-200  Signaling Pathway A       reactome   ~8e-04     ~1e-03
```

The log line `GSEA: 18 seed genes, 4/7 gene sets tested, 3 significant`
confirms the mechanics that matter: `R-HSA-100` ("Ancestor Broad
Signaling", Jaccard 0.86 against its child `R-HSA-101`) is dropped by
collapsing before testing (`Collapsed 1 near-duplicate gene set(s)`);
`R-HSA-500` (one seed) and `R-HSA-600` (no seed) overlap the signature too
little to be tested; `R-HSA-400` ("Scattered Noise Pathway", no coherent
enrichment signal) is tested but excluded by the FDR filter, so 3 of 4 are
significant. `disease_seed_genes.tsv` lists the 18 seed genes, strongest
first.

**Benchmark holdout** (`--benchmark-holdout-file`, for benchmark
validation): excludes genes from the seed set — both the GSEA signature and
`disease_seed_genes.tsv` — so a held-out gene is never a seed and can only
come back as a *candidate* in Stage 5, rather than validating itself. Stage 4
also ignores the held-out genes' own Open Targets evidence (read from
`benchmark_holdout_genes.tsv`), so it weights neither their edges nor Stage
6's `genetic_evidence` term. Genes are listed one per line
(`#` comments allowed) and resolved to Ensembl IDs via Stage 1's
`genes.parquet` — symbol, synonym, or bare Ensembl ID (`stage1_ingest.resolve_gene`) — never matched directly against
`ot_disease_subset` (which carries no symbol column at all). Example:

```
python3 pipeline/stage2_gsea_discovery.py --manifest example_data/example_run/manifest.json \
    --benchmark-holdout-file example_data/benchmark_holdout_example.txt
```

`example_data/benchmark_holdout_example.txt` holds `OLDG17`, a *synonym*
(not the approved symbol) for the gene at the core of `R-HSA-300`. Excluding
it drops `R-HSA-300` to 4 genes — one below blitzgsea's minimum tested-set
size — so it disappears from `disease_pathways.tsv` entirely (`3/7 gene sets
tested`, down from 4) rather than merely losing significance. G17 is no
longer a seed (17 seeds), but `R-HSA-300` still contains seeds (G18–G21), so
it stays in the Stage 3 graph and G17 is ranked as a candidate in Stage 5
(4th of 8 with `--method topology` on the Stage 3 graph).

`python3 -m pytest tests/test_gsea_discovery.py` runs the unit-level GSEA
tests (`tests/test_gsea_discovery.py`, 8 tests): a deliberately-enriched
pathway is returned significant, an individually-significant pathway is
correctly excluded once BH-FDR is applied across the tested library, non-
enriched pathways are tested but excluded from output, and a near-duplicate
ancestor/child pair collapses to the smaller child before GSEA ever sees the
ancestor.

## Running Stage 3 (`pipeline/stage3_build_graph.py`)

**Setup:** `pip install networkx scipy pandas pandera pydantic pyarrow` (or `pip install -r requirements.txt`).

**Run command**, after Stage 1 and Stage 2 have produced `manifest.json`, `disease_pathways.tsv` and `disease_seed_genes.tsv`:

```
python3 pipeline/stage3_build_graph.py --manifest out/manifest.json
```

Reads `gene_sets.parquet` and `interactions.parquet` via the manifest, and
`disease_pathways.tsv` and `disease_seed_genes.tsv` alongside it by default
(override with `--disease-pathways` / `--seed-genes`). The graph's nodes are
all genes of the union of (a) the Stage 2 disease pathways and (b) every
size-capped pathway that contains at least one seed gene — so each seed
gene, its pathway-mates, and every gene on the disease pathways are in the
graph. A Reactome interaction with no direction (`-` in the FI file) is
added in both directions; a directed one only in its own direction. Writes `graph_weight.npz`, `graph_sign.npz`,
`graph_gene_index.json`, and `graph_metadata.json` next to the manifest —
not GraphML; see the module docstring for why.

Against the checked-in example (regenerate with
`python3 example_data/build_gsea_example.py` if needed):

```
python3 pipeline/stage3_build_graph.py --manifest example_data/example_run/manifest.json
```

**Expected output** (log line):

```
INFO build_graph: Genes on the 3 disease-relevant pathways: 19.
INFO build_graph: Genes on the 6 pathways containing a seed gene: 25.
INFO build_graph: Graph: 25 nodes, 258 edges (6 pathways in the union; 4 interaction edges added, 1 out of scope).
INFO build_graph: Seed genes in the graph: 18/18.
```

The two gene counts overlap (here every disease pathway also contains
seeds), so the graph's node count is the size of their union, not their sum.

Six of the seven pathways contain a seed gene, so they are all in the graph
— including `R-HSA-100` (dropped by Stage 2's collapsing, but it contains
seeds) and `R-HSA-500` (too few seeds for GSEA). Only `R-HSA-600`, the
seedless pathway, is left out: its genes `N29`–`N33` are not nodes, and the
one dropped interaction row (`N29 -> N22`) is dropped for that reason —
proof that interaction edges really are scoped to the pathway-gene pool,
not added genome-wide. `graph_metadata.json` records `n_seed_pathways`,
`n_disease_pathway_genes`, `n_seed_pathway_genes`, `n_seed_genes_in_graph`
and `n_seed_genes_total`.

On a real disease with many genetic-association genes, rule (b) can pull in
a large share of Reactome; check `scale_report.json`'s projected edge count
and tighten `--max-set-size` in Stage 1 if the graph gets too big.

`python3 -m pytest tests/test_build_graph.py` runs the unit-level graph
tests (6 tests): the exact directed/signed edge list (with `source_db`
tags) for a hand-built pathway + interaction fixture, an undirected
interaction coming out symmetric, the co-membership
weight formula for both a small and a large (overlapping) pathway, that
a seed gene alone in its own pathway survives as an isolated node, and that
the pathway union is exactly disease pathways + seed-containing pathways (a
pathway with neither stays out).

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
`pipeline/stage1_ingest.py`'s `DATASOURCE_TO_DATATYPE`), so including it would let
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
separate, human-facing `gene_weights.tsv` that Stage 5 never reads.
Genes listed in Stage 2's `benchmark_holdout_genes.tsv` (next to the
manifest; override with `--holdout-genes`) score 0.0 whatever their Open
Targets rows say, and `graph_metadata.json` records them under
`stage4_genetic_evidence.holdout_genes`. See the
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

`ENSG...022` is not a seed but scores 0.5 from its `known_drug` row; the
other non-seed genes score 0.

`--edge-weight-mode product` (default `avg`) changes only the values written
into `graph_weight.npz`; e.g. the `ENSG...001`→`ENSG...002` edge is
`avg(0.60, 0.98) = 0.79` by default or `0.60 * 0.98 = 0.588` with `product`.
Note that with `product`, every edge touching a gene without evidence gets
weight 0, which cuts most non-seed genes off from RWR/BFWalk propagation.

## Running Stage 5 (`pipeline/stage5_score_candidates.py`)

**Setup:** `pip install networkx scipy numpy pandas pandera pydantic` (or `pip install -r requirements.txt`).

**Run command**, against a Stage 3 *or* Stage 4 graph directory (the two
share one format — see Stage 4 above):

```
python3 pipeline/stage5_score_candidates.py --graph-dir out/
```

`--graph-dir` is the only required argument: it reads `graph_weight.npz`,
`graph_sign.npz`, `graph_gene_index.json` and `graph_metadata.json` from
there, reads the seed genes from `disease_seed_genes.tsv` in the same
directory (override with `--seed-genes`, any TSV with a `gene_id` column),
and writes `candidate_scores.tsv` alongside (override with `--out`). Seed
genes that are not graph nodes are dropped (logged); every graph gene that
is not a seed is a candidate, so seeds never appear in the ranking.

`--method` is a comma-separated subset of `{topology, rwr, bfwalk}`,
**default `topology`**. `topology_score` is always emitted; `rwr_score` /
`bfwalk_score` appear only when `rwr` / `bfwalk` is requested (e.g.
`--method topology,rwr` or `--method rwr,bfwalk`). The output is ranked by
the last requested method in the order topology, rwr, bfwalk. All three
start from the seed set:

- `topology` = (1 + number of seed neighbours) / (1 + weighted
  shortest-path distance from the nearest seed, each edge costing
  1 / weight); unreachable genes score 0. A zero-weight edge (Stage 4,
  both genes without evidence) costs as much as the weakest positive-weight
  edge, so it is never cheaper than an edge with evidence.
- `rwr` = personalised PageRank restarting uniformly on the seeds.
  `--restart-prob` (default 0.5; PageRank `alpha = 1 - restart_prob`).

`bfwalk` runs [BFWalk](https://github.com/jedrzejkubica/BFWalk)
(backtrack-free network propagation with in-degree normalisation) with the
seed genes as its seeds, on the same directed weighted graph. BFWalk needs edge weights in
`]0, 1]`, so zero-weight edges (Stage 4 genes without evidence) are dropped —
as in RWR, they carry no flow — and the remaining weights are divided by the
largest one. Knobs: `--bfwalk-alpha` (attenuation coefficient, default 0.5)
and `--bfwalk-threads` (default 0 = all cores). BFWalk is not on PyPI yet:
`requirements.txt` installs it from TestPyPI (see the top-level README); it
will be available via `pip install bfwalk` shortly. Stage 6 carries
`bfwalk_score` through and can use it as a composite component
(`--weights ...,bfwalk:0.3`; it is not in the default weights), and Stage 7
shows it as a `bfwalk` column.

Against the checked-in example (after Stages 1–4 have run, so
`example_data/example_run/` holds a **Stage 4** graph):

```
python3 pipeline/stage5_score_candidates.py --graph-dir example_data/example_run --method topology,rwr
```

**Expected output** (`example_data/example_run/candidate_scores.tsv`, all
7 candidates; columns per `schemas.CandidateScoresSchema`, ranked by
`rwr_score` since `rwr` was requested — else by `topology_score`):

```
gene             topology_score  rwr_score
ENSG00000000022  5.5287          0.0207
ENSG00000000023  4.2752          0.0127
ENSG00000000024  2.3020          0.0068
ENSG00000000025  0.9866          0.0023
ENSG00000000026  0.4615          0.0011
ENSG00000000027  0.4615          0.0011
ENSG00000000028  0.4615          0.0011
```

`N22`/`N23` share two pathways with 10–12 strong seeds and lead;
`N24` shares one pathway with six seeds; `N25`–`N28` share a pathway with
a single mid-strength seed, and `N25` is lifted above its pathway-mates by
its interaction with seed `G2`. The 18 seeds are not in the table.

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
tests (21 tests): hand-computed `topology_score` values on a small
two-seed cycle+chain graph (both at the function level and end-to-end
through a default `run()` with no `--method`), an isolated node scoring 0,
a zero-weight edge never being cheaper than an edge with evidence, an
all-seed graph writing an empty ranking, `--restart-prob` outside ]0, 1[
exiting non-zero,
RWR restarting on the seeds and reconstructed as a full stationary
distribution (sums to 1), BFWalk following edge direction and ignoring
zero-weight edges, `--seed-genes` overriding the default seed file (and a
seed list without a `score` column), seeds outside the graph dropped, no
seed in the graph / missing seed file exiting non-zero, and the
output-column contract (`topology_score` always, `rwr_score`/`bfwalk_score`
only when requested, seed rows excluded, unknown `--method` exits non-zero).

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
ENSG00000000022  5.5287          0.5                     clinical      has_liabilities  1                     0.800
ENSG00000000023  4.2752          0.0                     unknown       unknown          0                     0.453
ENSG00000000024  2.3020          0.0                     discovery     unknown          0                     0.296
ENSG00000000025  0.9866          0.0                     clinical      unknown          0                     0.199
ENSG00000000026  0.4615          0.0                     unknown       unknown          0                     0.050
```

`ENSG...022` stays on top despite a safety liability because its proximity
terms and its `known_drug` evidence dominate; `ENSG...025`, clinically
tractable but close to only one seed, stays below `ENSG...023`/`024` —
tractability alone can't lift a low-proximity candidate. Candidates are
non-seed genes, so their `genetic_evidence_score` can only come from
datatypes other than `genetic_association` (here `known_drug`). Genes absent
from the synthetic `ot_target_subset.parquet` are `unknown`/`unknown` and
ranked purely on the score terms, never penalised for missing annotation.

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

**Orchestrator** — `run_pipeline.py` chains Stages 1–6 for one disease:

```
python3 pipeline/run_pipeline.py --disease MONDO_0005265 \
    --data-dir data/ --out-dir out/ --report
```

Each stage runs as its own subprocess (the same CLI as its "Running Stage
N" section above). If one exits non-zero, `run_pipeline.py` prints which
stage failed, its exit code and the exact command, and stops the run
non-zero — the stage's own message (e.g. Stage 1's "Could not resolve
`--disease 'FOO'` …") is right above it, not buried under a traceback. `--report` additionally writes `out/report.html` at the end.

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

- header — disease, number of seed genes and candidates, random seed, Reactome / Open
  Targets versions, whether the graph is Stage 3 (structural) or Stage 4
  (genetic-evidence-weighted);
- a **sortable** ranked table — `symbol`, `gene` (Ensembl), `composite`,
  `topology`, `rwr`, `genetic ev.`, `tractability`, `safety` (click a
  header to sort; `rwr`/`genetic ev.` columns appear only if those stages
  ran);
- a **click-to-expand evidence trace** for each of the top 100 candidates
  (`--n-traces`; the table still lists every candidate) — its
  `composite_breakdown` + weights, the graph pathways it shares with seed
  genes (with the number of seeds in each; and the disease-relevant pathways
  without seeds it is a member of), its Open
  Targets datatypes + scores (with the ones feeding `genetic_evidence_score`
  flagged; none for a benchmark holdout gene), and the signed interaction edges touching it (edges with a seed
  are badged "with seed"; interactions to non-graph genes are omitted,
  matching Stage 3's scoping).

Genes are displayed by HGNC symbol; the Ensembl gene ID stays the internal
key (its own column, and the anchor id for each detail block). Symbols are
resolved from `genes.parquet` at render time only.

Only `manifest.json` and `candidates_annotated.tsv` are required; every
other artifact is optional. A partial run directory still produces a
report — missing trace inputs (`gene_sets.parquet`, `interactions.parquet`,
`ot_disease_subset.parquet`, …) show that section as unavailable, a
Stage-3-only run (no Stage 4) just drops the `rwr`/`genetic ev.` columns,
and no `genes.parquet` falls the symbol column back to Ensembl ids — a
header line lists whatever was missing.

**Expected output** (`example_data/example_run/report.html`, top of the
ranked table — same ordering as Stage 6's `candidates_annotated.tsv`):

```
#  symbol  gene             composite  topology  rwr     genetic ev.  tractability  safety
1  GENE22  ENSG00000000022  0.800      5.529     0.0207  0.500        clinical      has_liabilities
2  GENE23  ENSG00000000023  0.453      4.275     0.0127  0.000        unknown       unknown
3  GENE24  ENSG00000000024  0.296      2.302     0.0068  0.000        discovery     unknown
```

Expanding `GENE25` shows: pathways shared with seed genes *Seed-Adjacent
Pathway* (1 seed gene); no Open Targets rows for this disease; interaction
`GENE02 → GENE25` (activating, curated, confidence 0.80, badged "with seed").

`python3 -m pytest tests/test_stage7_report.py tests/test_run_pipeline.py`
runs the tests (12): `evidence_trace`'s seed-pathway / graph-scoping /
datatype logic on a hand-built run, the `feeds` flag going false when Stage
4 did not run or the gene is a benchmark holdout, `build_report` producing a
self-contained symbol-labelled HTML with one detail block per candidate,
detail blocks (and symbol links) only for the top `--n-traces`, a partial run
directory (only `manifest.json` + `candidates_annotated.tsv`) still
rendering with the optional columns/trace sections marked unavailable, a
`run_pipeline --from-manifest --report` smoke test that Stages 2–6 chain
and produce every expected artifact, and a stage-failure producing a
helpful "which stage / which command / exit code" message rather than a
traceback.
