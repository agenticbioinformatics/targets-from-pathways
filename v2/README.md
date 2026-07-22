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
overlap by adding network-diffusion (physics-inspired) propagation models —
standard random-walk-with-restart and a backtrack-free (non-backtracking)
walk variant — alongside topological distance, and by explicitly
cross-validating pathway evidence against independent, non-pathway-derived
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
the backtrack-free walk scoring method, end-to-end identifier
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
| 5. Candidate scoring | Swappable `--method`: topology score (co-membership / shortest path / branch-convergence, direction- and weight-aware), random-walk-with-restart (RWR), or backtrack-free walk (BFW, a non-backtracking RWR variant, via [GBA-centrality](https://github.com/jedrzejkubica/GBA-centrality)) — all from the target node, all able to consume Stage 4's node/edge weights and Stage 3's edge directions/signs. | H1, H2, H3 |
| 6. Contextual annotation & filtering | Attach Open Targets tractability bucket and safety flags; compute a configurable composite score. No tissue-expression filtering — the method stays strictly pathway- and association-evidence-based. | H5, H6 |
| 7. Benchmark validation | Score a small literature-curated set of known resistance/compensation gene pairs (held out of Stage 2's seed genes); report descriptive rank/percentile recovery against a random background, not a significance test. | H8 |
| 8. Output/report | Ranked table (per-candidate scores, evidence trace, tractability/safety) plus a web UI: target + disease in, ranked list out. | Deliverable |

Hypotheses H9 (target-modality generalization) and H10 (disease-context
expression weighting) were evaluated and explicitly dropped from hackathon
scope — H10 was removed entirely (not kept as an optional filter) per the
first plan update below. H7 (cross-pathway-DB consensus) was briefly
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
  algorithms.
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
- Stage 0 (`schemas.py`: the inter-stage table contracts, written once and
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
- Stage 5 (topology + RWR + BFW candidate scoring, all weight-aware)
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
