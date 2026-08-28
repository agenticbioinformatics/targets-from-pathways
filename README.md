# Target Rescue via Pathways (v2) — Hackathon Plan

> Fresh redesign of the pathway-based alternative-target tool, building on the
> lessons and open questions from the original [targets-from-pathways](../README.md)
> project (Open Targets Hackathon, Oct 2025), but not bound to its specific
> method choices.

## Introduction

When a therapy acting on a specific molecular target stops working for a
patient, clinicians and researchers need alternative or compensatory targets
that could restore treatment efficacy, but systematically finding them is
hard without patient-specific multi-omics data. Biological pathways encode a
public, disease-independent map of how genes functionally relate to a target,
so our core research question is whether pathway topology alone — without
patient-specific omics — can provide a credible, orthogonal line of evidence
for alternative-target discovery. We sharpen this beyond simple pathway
overlap by adding a network-diffusion (physics-inspired) propagation model —
standard random-walk-with-restart — alongside topological distance, and by
explicitly cross-validating pathway evidence against independent,
non-pathway-derived
genetic evidence to avoid circular reasoning — a rigor gap the original
project's own scoring formula did not address. The tool uses only public
data: Open Targets Platform for disease-target genetic associations,
tractability, and safety, and Reactome for pathway gene sets and
directional, signed functional interactions; no patient-specific data is
required. Given a target of interest and a disease,
the pipeline finds disease-relevant pathways, builds a **pathway-based
gene-gene graph** spanning the target's and the disease's pathway
neighborhoods — a directed graph with positive/negative (activating/
inhibiting) edges sourced from Reactome's functional-interaction
annotations — incorporates independent genetic evidence directly as node and
edge weights on that graph, scores candidate genes by weighted topological
and diffusion proximity to the target, and annotates candidates with
tractability and safety information so results are actionable, not just
statistically interesting. The pipeline runs on Reactome alone for the
hackathon, but every stage consumes database-agnostic normalized tables
carrying a `source_db` column, so a second curation (WikiPathways, OmniPath,
KEGG) can later be added as a single adapter without touching downstream
stages — cross-database comparison is a documented stretch goal rather than
core scope. What's novel here relative to the original project is the
explicit orthogonality discipline (excluding pathway-derived Open Targets
evidence from the weighting step), the directed/signed graph representation,
end-to-end identifier
canonicalisation with mapping-coverage reporting, and the addition of a
literature-curated benchmark of known clinical resistance mechanisms (e.g.
EGFR→MET) to sanity-check the method against ground truth before trusting
its output. The expected deliverable is
a small web app: a researcher enters a target and a disease and gets back a
ranked, evidence-annotated list of candidate alternative targets. Known
limitations going in: single-database (Reactome) coverage and curation bias,
uncorrected in this version because no second curation is run; sparse
Open Targets safety annotation coverage for many genes, a benchmark set small
enough (5–15 cases) that validation results must be read as descriptive, not
statistically powered, and the deliberate exclusion of patient-specific
expression data entirely — no tissue-expression filter is included, keeping
the method strictly pathway- and association-evidence-based. This tool is
meant to complement, not replace, genetic and clinical evidence in target
prioritization.

## Pipeline overview

| Stage | Purpose | Addresses |
|---|---|---|
| 1. Ingest & canonicalize | Take target gene + disease (EFO ID); resolve source data by version-pinned, checksummed download (or validate an existing `--data-dir`); **canonicalize every gene identifier to Ensembl gene ID** using the Open Targets `target` parquet as the ID authority; subset OT associations to the disease; normalize Reactome gene sets and functional interactions into two database-agnostic tables; emit a manifest plus mapping-coverage and graph-scale reports. | Setup, H1 |
| 2. Disease-pathway discovery | GSEA (FDR-corrected) of disease-associated genes against Reactome gene sets, over size-capped and redundancy-collapsed sets → disease-relevant pathway list. | H1 |
| 3. Pathway-based gene-gene graph construction | Build a **directed, signed** pathway-based gene-gene graph (`networkx.DiGraph`) from the union of disease-relevant and target-relevant pathway co-membership — with co-membership edges down-weighted by pathway size — plus directional, positive/negative (activating/inhibiting) relation edges from Stage 1's normalized Reactome functional-interaction table; DB versions and seed carried from the manifest. | H1, H2 |
| 4. Genetic-evidence weighting | Compute Open Targets evidence restricted to non-pathway-derived datatypes (e.g. `genetic_association`, `known_drug`, explicitly excluding `affected_pathway`/literature-pathway sources) and map the resulting scores onto the Stage 3 graph as **node weights and edge weights** — so genetic evidence directly informs Stage 5's scoring rather than only re-ranking its output after the fact, while preserving orthogonality with Stages 2–3. | H4 |
| 5. Candidate scoring | `--method` defaults to **topology score only** (co-membership / shortest path / branch-convergence, direction- and weight-aware) — the cheapest, deterministic method, run for every candidate by default. Random-walk-with-restart (RWR) is opt-in via `--method` (e.g. `--method topology,rwr`), not run unless requested. Both score from the target node and are able to consume Stage 4's node/edge weights and Stage 3's edge directions/signs. | H1, H2 |
| 6. Contextual annotation & filtering | Attach Open Targets tractability bucket and safety flags; compute a configurable composite score. No tissue-expression filtering — the method stays strictly pathway- and association-evidence-based. | H5, H6 |
| 7. Benchmark validation | Score a small literature-curated set of known resistance/compensation gene pairs (held out of Stage 2's seed genes); report descriptive rank/percentile recovery against a random background, not a significance test. | H8 |
| 8. Output/report | Ranked table (per-candidate scores, evidence trace, tractability/safety) plus a web UI: target + disease in, ranked list out. | Deliverable |

Hypotheses H9 (target-modality generalization) and H10 (disease-context
expression weighting) were evaluated and explicitly dropped from hackathon
scope — H10 was removed entirely (not kept as an optional filter) per the
first plan update below. H3 (backtrack-free / non-backtracking walk beats
standard RWR) was in scope through plan update 4, then dropped in plan
update 5: Stage 5 now ships topology + standard RWR only, with no
non-backtracking-walk variant and no external walk-simulation dependency.
H7 (cross-pathway-DB consensus) was briefly
promoted into core scope, then returned to **stretch goal** status in plan
update 3: the hackathon build runs Reactome only, but Stage 1's normalized
`source_db`-tagged tables keep the seam open so a second curation can be
added later as one adapter.

### Review changelog (Phase 6, converged after 3 of 5 iterations)

- **Iter 1:** merged target-pathway context into graph construction; merged
  tractability/safety/expression into one annotation stage; added FDR
  correction requirement to GSEA; required benchmark genes be held out of
  the seed-gene list; required using existing libraries for GSEA and
  diffusion rather than reimplementing; required pinned DB versions and
  fixed seeds.
- **Iter 2:** caught that Open Targets' aggregated association score is
  partly pathway-derived — restricted Stage 5 to non-pathway datatypes only,
  to avoid validating pathway evidence against itself.
- **Iter 3:** no substantive issues found — converged, stopped early.

### Plan update (user-directed, post-review)

- **Stage 5 (candidate scoring):** added backtrack-free walk (BFW /
  non-backtracking random-walk-with-restart) as a third scoring method
  alongside topology and standard RWR, using the existing
  [GBA-centrality](https://github.com/jedrzejkubica/GBA-centrality) tool
  rather than reimplementing non-backtracking walk logic — consistent with
  the earlier review requirement to reuse existing libraries for graph
  algorithms. *(Reversed by plan update 5 — BFW and the GBA-centrality
  dependency were removed from scope.)*
- **Stage 4 (was "genetic-evidence integration"), reordered ahead of
  candidate scoring:** changed from post-hoc re-ranking of Stage 5's output
  to computing node weights (and optional edge weights) from Open Targets
  association scores, still restricted to non-pathway datatypes for
  orthogonality, which now feed directly into Stage 5's topology and
  diffusion computations.
- **Stage 6 (contextual annotation & filtering):** removed the optional
  GTEx tissue-expression filter entirely; H10 is dropped from scope rather
  than kept as an off-by-default option.

### Plan update 2 (user-directed)

- **Stage 3, renamed "pathway-based gene-gene graph construction":** the
  graph is now explicitly directed and signed (positive/negative,
  activating/inhibiting edges), sourced from Reactome's functional-
  interaction annotations and, where selected, KEGG's relation types.
  "Gene-gene graph" is renamed "pathway-based gene-gene graph" throughout to
  make the pathway provenance explicit.
- **Stages 1–3 (and, by extension, 2 and 7):** added multi-database support
  — `--pathway-db` now accepts `reactome`, `kegg`, or `reactome,kegg`, so
  GSEA and graph construction can run on either database alone or their
  union. H7 (cross-pathway-DB consensus) is promoted from Backup to core
  scope: Stage 7's benchmark now reports validation results per database
  configuration side by side, which is the actual comparison H7 called for.
- **Stage 4 (genetic-evidence weighting):** edge weighting is no longer
  optional — both node and edge weights are computed from Open Targets data
  by default (previously edge weights defaulted to off).

### Plan update 3 (user-directed, Stage 1 redesign)

Three decisions, and a rewrite of Stage 1 to match.

**Decision 1 — Reactome only for the hackathon.** Multi-database support is
withdrawn from core scope (see H7 note above). This is a deliberate scope
trade, not a retreat: it removes one GSEA run and the combined run from
Stage 2, removes the need for a `MultiDiGraph` with parallel conflicting
edges in Stage 3 (cross-DB sign conflicts cannot arise with one source, so a
plain `DiGraph` suffices), cuts Stage 7 from ~45 full pipeline runs to
~10–15, and cuts Stage 8's side-by-side database toggle. The recovered time
is explicitly reallocated to the degree-matched null model and the
non-pathway baseline comparison, without which the ranked output cannot be
distinguished from a list of network hubs. The `source_db` column and the
adapter seam are retained throughout even though only one value is ever
emitted, so adding a source later is a new file rather than a refactor.

**Decision 2 — Ensembl gene ID is the canonical internal identifier.** Open
Targets is the primary evidence source and is natively Ensembl-keyed, and
Ensembl IDs are stable across symbol renames. Reactome symbols are mapped in
at Stage 1 and mapped back out to symbols only for display in Stages 6/8.

**Decision 3 — Stage 1 may download its own inputs**, from a version-pinned
registry with sha256 verification and a `--no-download` opt-out, replacing
the previous print-instructions-only rule. Reproducibility is then enforced
rather than documented, and a new team member gets a working data directory
in one command.

**Decision 4 — pin Open Targets Platform release 26.06** (released
2026-06-24), rather than the 25.09 release used by v1. Pin the release
*number*: `.../platform/latest/` is a moving pointer and would silently
change results between runs, which is the opposite of what the registry is
for. Note the consequence — v1's published PTGS2 / rheumatic-disease numbers
were computed on 25.09 and are therefore a loose sanity reference, not a
reproducible comparison; do not treat a difference against them as a bug
without first checking whether the underlying associations changed.

**Stage 1, renamed "Ingest & canonicalize":** promoted from a file-existence
checker to the resolve → normalize → canonicalize stage. It is now the only
place gene identifiers are mapped, which is where that work belongs — v1's
published top-10 contains `GRB2-1`, a Reactome *physical entity* identifier
leaking into a gene list, which is exactly the class of bug a single
canonicalisation point prevents. Stage 1 now also emits a mapping-coverage
report (fail the run below a coverage floor) and a graph-scale report
projecting Stage 3's edge count from the gene-set size distribution, so an
infeasible graph is caught on Day 1 morning rather than on Day 2 evening.
Its outputs are two database-agnostic tables (`gene_sets`, `interactions`)
plus a canonical `genes` table, a disease-subset association table, and a
manifest; downstream stages read only these, never raw source files.

### Plan update 4 (user-directed)

- **Stage 5 (candidate scoring):** `--method` now defaults to topology
  score only; RWR and BFW are opt-in, not run by default. Topology is
  deterministic and cheap (no permutations, no walk simulation), so it's
  the right always-on baseline; RWR and BFW are heavier (RWR needs a
  converged stationary distribution, BFW needs GBA-centrality's walk
  simulation) and are requested explicitly via `--method` when a
  diffusion-based comparison is actually wanted, rather than computed on
  every run whether or not anything downstream uses them. *(Superseded in
  part by plan update 5 — BFW was subsequently removed; the
  topology-default / opt-in split still applies to RWR.)*

### Plan update 5 (user-directed)

- **Stage 5 (candidate scoring):** backtrack-free walk (BFW /
  non-backtracking RWR) and its external
  [GBA-centrality](https://github.com/jedrzejkubica/GBA-centrality)
  dependency are **removed from scope entirely**, reversing the Stage 5
  bullet in "Plan update (user-directed, post-review)". GBA-centrality
  needs a compiled C extension and its own git submodule, a build step out
  of proportion to a hackathon prototype and a third scoring method the
  core hypotheses (H1/H2) don't require. Stage 5 ships **topology (default)
  + standard RWR (opt-in)** only; `--method` accepts `topology`, `rwr`, or
  `topology,rwr`. H3 is dropped (see the hypotheses note above).

## Open decisions

Nothing in this section is agreed. Each item names the placeholder the plan
currently assumes, and each **requires research and a decision owner** before
Day 1 code is final.

### 1. What counts as "pathway topology"?

The project claims pathways only, which is what rules out Open Targets'
`interaction` table — that is protein-protein interaction data and carries
neither sign nor direction. But the same objection applies to the source the
plan currently uses: Reactome's Functional Interaction network (Wu et al.
2010) merges curated reactions with **machine-learning-predicted interactions
trained on PPI, co-expression and domain-domain features**. Using
`FIsInGene_*` unfiltered breaks the pathways-only claim exactly as a PPI
network would.

*Suggested rule, if adopted:* an edge is admissible iff it is derivable from
a curated Reactome reaction or complex — no PPI assay, co-expression, text
mining, or ML prediction.

Reactome's reaction graph, projected to genes, would supply direction and
sign from curation alone:

| Reactome relation | Edge | Sign |
|---|---|---|
| A input to reaction R, B output of R | A → B | + |
| A catalyses R, B output of R | A → B | + |
| A positive / negative regulator of R | A → outputs(R) | + / − |
| A, B in the same complex | A — B | 0 |

*Options:* curated-only FI (**current placeholder**) · Pathway Commons SIF
filtered to `datasource = reactome` · Reactome BioPAX Level 3, sign from
`Control.controlType` · co-membership only, no interaction edges.

*Trade-off:* reaction-derived graphs are far sparser than FI's ~230k edges —
less hub domination, but lower candidate recall and some benchmark pairs may
become unreachable.

### 2. Correcting Stage 2 for annotation bias

GSEA assumes exchangeable gene labels. Well-studied genes carry both more
Open Targets evidence *and* more Reactome annotation, so enrichment inflates
systematically rather than randomly — hardest for the large, well-curated
pathways least informative about any specific disease. v1 returned 250
pathways and 3,259 candidates, which is barely discriminative. The Jaccard
collapsing already in Stage 2 fixes redundancy, not bias.

*Options:* BH-FDR alone (**current placeholder**) · disease-label permutation
null — run the identical enrichment for ~200 other EFO diseases with
comparable gene counts and report an empirical specificity p per pathway ·
a cheaper cached "promiscuity score", the fraction of a reference disease
panel for which a pathway is significant, excluding the top decile.

The permutation null materially increases Stage 2 runtime, which matters for
Stage 7's repeated runs.

### 3. Capping the target-pathway union in Stage 3

Stage 3 builds from (disease-relevant pathways) ∪ (target-containing
pathways). The size cap constrains the first term only. If the target sits in
one very large pathway — PTGS2 is in several — that single membership pulls
much of Reactome into the graph.

*Options:* no cap (**current placeholder**) · apply the same
`--min/--max-set-size` to the target side · keep only the *k* smallest
pathways containing the target · exclude Reactome's top-level hierarchy
roots. Not mutually exclusive.

### Open research

Needed before Stage 1 can be finished. All are measurable or checkable on
Day 1 — Stage 1's `coverage_report.json` and `scale_report.json` are designed
to answer the last three directly.

- Do Reactome **numbered-release URLs** exist for `ReactomePathways.gmt` and
  the FI file? If not, pin by recorded sha256 rather than by URL.
- The FI file's **real schema** — and is curated cleanly separable from
  predicted? Assumed to be `Gene1, Gene2, Annotation, Direction, Score` with
  curated at `Score == 1.0`; unverified.
- Does **Pathway Commons SIF retain BioPAX `controlType`** (ACTIVATION /
  INHIBITION)? This alone decides whether that option is viable at all.
- Is the Reactome mapping file **PE-level or gene-level**, and what is a
  correct PE→gene collapse? (v1 leaked `GRB2-1` by getting this wrong.)
- The real **Reactome→Ensembl mapping rate** — the 90% coverage floor is a
  placeholder.
- The real **set-size distribution and projected edge count** — the
  `--max-set-size` of 200 is GSEA convention, not a derived number. The
  Jaccard threshold of 0.7 is likewise a placeholder.

## Development plan

Sized for a small (2–4 person) mixed bio+CS team on laptop-scale compute.

**Day 1 — Core graph pipeline**
- Stage 0 (`pipeline/stage0_schemas.py`: the inter-stage table contracts, written once and
  validated on read *and* write by every stage — ~45 min, and the single
  highest-leverage hour of the three days)
- Stage 1 (ingest & canonicalize: pinned download, Ensembl canonicalisation,
  normalized Reactome tables, manifest, coverage + scale reports)
- Stage 2 (GSEA disease-pathway discovery, FDR-corrected)
- Stage 3 (directed, signed pathway-based gene-gene graph construction,
  versioned/seeded)
- In parallel, no code required: curate the Stage 7 benchmark TSV of known
  resistance pairs with verified PMIDs. This is literature work, it gates
  Day 3, and it is the one file that must not be agent-generated.
- Milestone: for a toy (target, disease) pair, produce a disease-relevant
  pathway list and a constructed pathway-based gene-gene graph, with Stage
  1's coverage report showing >90% identifier mapping for every source.

**Day 2 — Weighting and scoring**
- Stage 4 (genetic-evidence node/edge weighting, non-pathway datatypes only)
- Stage 5 (topology + RWR candidate scoring, both weight-aware)
- Stage 6 (tractability + safety annotation, composite score)
- Milestone: full ranked, annotated candidate list on a real example (e.g.
  rheumatic disease / PTGS2, reusing the original project's example as a
  known-reasonable sanity check).

**Day 3 — Validation, integration, demo**
- Stage 7 (benchmark validation against curated resistance cases)
- Stage 8 (output/report + minimal web UI)
- End-to-end integration test, documentation pass, demo rehearsal buffer.

## Testing steps

- **Per-stage sanity checks**: a tiny synthetic fixture (~15–20 genes, 3–5
  pathways with a hand-constructed, known-correct answer, including at least
  one deliberate positive and one negative edge) run through each stage
  independently, checking exact expected output (e.g. Stage 3's graph has
  the expected directed, signed edges, Stage 5's topology score matches a
  hand-computed value).
- **Identifier-mapping check**: assert Stage 1 maps a clean HGNC symbol, a
  deprecated symbol via the Open Targets synonym field, and records an
  unmappable symbol in the dropped list with the correct count — the check
  that would have caught v1's `GRB2-1`.
- **Manifest determinism check**: the same inputs produce a byte-identical
  manifest modulo timestamp, and every recorded artifact sha256 matches the
  file on disk.
- **End-to-end integration run**: full pipeline on a real (target, disease)
  pair (e.g. PTGS2 / rheumatic disease, EFO_0005755) checked for a
  non-degenerate ranked output and sane runtime.
- **Ground-truth validation**: Stage 7's benchmark run against the
  literature-curated resistance/compensation pairs — this is the check that
  catches pipeline bugs before any
  biological conclusion is trusted, and the one place small-n results must
  be reported descriptively (e.g. "N of 10 known pairs ranked in top 5%"),
  not as a statistical test.

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
python3 pipeline/stage1_ingest.py --target NOD2 --disease MONDO_0005265 --data-dir ./test_data --out-dir ./test_out
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
python3 pipeline/stage2_gsea_discovery.py --manifest test_out/manifest.json
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

Log lines confirm the mechanics that matter: `R-HSA-100` ("Ancestor Broad
Signaling", Jaccard 0.83 against its child `R-HSA-101`) is dropped by the
collapsing step before testing (`5 sets before, 4 sets after`), and
`R-HSA-400` ("Scattered Noise Pathway", no coherent enrichment signal) is
tested but excluded by the FDR filter (`3/4 tested pathways pass`).

**Benchmark holdout** (`--benchmark-holdout-file`, for Stage 7): excludes
genes from the disease-associated seed set before GSEA runs, so a
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
python3 pipeline/stage3_build_graph.py --manifest test_out/manifest.json
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

**Expected output** (log lines):

```
INFO build_graph: Pathway union: 3 disease-relevant + 2 target-containing = 4 unique pathways (1 overlap both).
INFO build_graph: Co-membership: 4 pathways contributed edges (of 4 in the union; ...), inducing 68 unordered gene pairs (136 directed edges).
INFO build_graph: Interactions: 2 rows in scope (both genes in the pathway union) added/updated as edges; 1 row(s) out of scope dropped (...); 0 duplicate (...) pair(s) collapsed ...
INFO build_graph: Pathway-based gene-gene graph: 17 nodes, 136 directed edges.
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
python3 pipeline/stage4_genetic_evidence_weights.py --manifest test_out/manifest.json
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
python3 pipeline/stage5_score_candidates.py --graph-dir test_out/
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
are a display-time lookup (`genes.parquet`) in Stages 6/8, never written
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
python3 pipeline/stage6_annotate_context.py --manifest ./run1/manifest.json
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
  self-contained explanation for Stage 8's report.

Per-candidate **evidence trace**: the passthrough score columns
(`topology_score` / `rwr_score` / `genetic_evidence_score`), the annotation
buckets (`tractability` / `safety` / `n_safety_liabilities`), and
`composite_breakdown` + `composite_weights` together explain every rank.
Two provenance items are deliberately **not** here — *which shared
pathways/interactions* (not in this stage's inputs; Stage 8 assembles it
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
