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
tractability, and safety, and Reactome (optionally combined with KEGG) for
pathway gene sets and directional, signed functional interactions; no
patient-specific data is required. Given a target of interest and a disease,
the pipeline finds disease-relevant pathways, builds a **pathway-based
gene-gene graph** spanning the target's and the disease's pathway
neighborhoods — a directed graph with positive/negative (activating/
inhibiting) edges sourced from Reactome's (and optionally KEGG's) relation
annotations — incorporates independent genetic evidence directly as node and
edge weights on that graph, scores candidate genes by weighted topological
and diffusion proximity to the target, and annotates candidates with
tractability and safety information so results are actionable, not just
statistically interesting. Because the graph-construction and pathway-
discovery stages are parameterized by pathway database, the same pipeline
can be run on Reactome alone, KEGG alone, or their union, letting us compare
results across databases rather than committing to one curation's biases.
What's novel here relative to the original project is the explicit
orthogonality discipline (excluding pathway-derived Open Targets evidence
from the weighting step), the directed/signed graph representation, the
backtrack-free walk scoring method, the cross-database comparison, and the
addition of a literature-curated benchmark of known clinical resistance
mechanisms (e.g. EGFR→MET) to sanity-check the method against ground truth
before trusting its output. The expected deliverable is
a small web app: a researcher enters a target and a disease and gets back a
ranked, evidence-annotated list of candidate alternative targets. Known
limitations going in: pathway database coverage and curation bias, sparse
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
| 1. Ingest | Take target gene + disease (EFO ID) + `--pathway-db` config (`reactome`, `kegg`, or `reactome,kegg`); fetch/cache Open Targets and the selected pathway database(s)' data. | Setup |
| 2. Disease-pathway discovery | GSEA (FDR-corrected) of disease-associated genes against pathway gene sets, run per selected database and, when both are selected, also on their combined gene sets → disease-relevant pathway list(s), tagged by source DB. | H1, H7 |
| 3. Pathway-based gene-gene graph construction | Build a **directed, signed** pathway-based gene-gene graph from the union of disease-relevant and target-relevant pathway co-membership plus directional, positive/negative (activating/inhibiting) relation edges parsed from Reactome's functional-interaction annotations and/or KEGG's relation types (activation/inhibition/expression/repression), depending on `--pathway-db`; pin DB version(s), fix random seeds. | H1, H2, H7 |
| 4. Genetic-evidence weighting | Compute Open Targets evidence restricted to non-pathway-derived datatypes (e.g. `genetic_association`, `known_drug`, explicitly excluding `affected_pathway`/literature-pathway sources) and map the resulting scores onto the Stage 3 graph as **node weights and edge weights** — so genetic evidence directly informs Stage 5's scoring rather than only re-ranking its output after the fact, while preserving orthogonality with Stages 2–3. | H4 |
| 5. Candidate scoring | Swappable `--method`: topology score (co-membership / shortest path / branch-convergence, direction- and weight-aware), random-walk-with-restart (RWR), or backtrack-free walk (BFW, a non-backtracking RWR variant, via [GBA-centrality](https://github.com/jedrzejkubica/GBA-centrality)) — all from the target node, all able to consume Stage 4's node/edge weights and Stage 3's edge directions/signs. | H1, H2, H3 |
| 6. Contextual annotation & filtering | Attach Open Targets tractability bucket and safety flags; compute a configurable composite score. No tissue-expression filtering — the method stays strictly pathway- and association-evidence-based. | H5, H6 |
| 7. Benchmark validation | Score a small literature-curated set of known resistance/compensation gene pairs (held out of Stage 2's seed genes) once per `--pathway-db` configuration (Reactome-only, KEGG-only, combined); report descriptive rank/percentile recovery vs. a random background per configuration, side by side, not a significance test. | H8, H7 |
| 8. Output/report | Ranked table (per-candidate scores, evidence trace, tractability/safety, source pathway DB) plus a web UI: target + disease + pathway-DB selection in, ranked list out. | Deliverable |

Hypotheses H9 (target-modality generalization) and H10 (disease-context
expression weighting) were evaluated and explicitly dropped from hackathon
scope — H10 was removed entirely (not kept as an optional filter) per the
first plan update below. H7 (cross-pathway-DB consensus), originally scored
**Backup**, has since been promoted into core scope per a later plan update:
the pipeline now runs on Reactome, KEGG, or their combination so results can
be compared across databases directly.

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

## Development plan

Sized for a small (2–4 person) mixed bio+CS team on laptop-scale compute.

**Day 1 — Core graph pipeline**
- Stage 1 (ingest Open Targets + Reactome, and KEGG if time allows)
- Stage 2 (GSEA disease-pathway discovery, FDR-corrected, per database)
- Stage 3 (directed, signed pathway-based gene-gene graph construction,
  versioned/seeded; start with Reactome only, add KEGG once the Reactome
  path works end to end)
- Milestone: for a toy (target, disease) pair, produce a disease-relevant
  pathway list and a constructed pathway-based gene-gene graph, at least for
  Reactome.

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
- **Multi-database check**: run Stages 1–3 with `--pathway-db reactome`,
  `--pathway-db kegg`, and `--pathway-db reactome,kegg` on the same toy
  fixture, and confirm the combined-mode graph is the expected union (no
  duplicated or dropped edges) rather than testing only one database and
  assuming the others work.
- **End-to-end integration run**: full pipeline on a real (target, disease)
  pair (e.g. PTGS2 / rheumatic disease, EFO_0005755) checked for a
  non-degenerate ranked output and sane runtime, for at least the
  Reactome-only configuration.
- **Ground-truth validation**: Stage 7's benchmark run against the
  literature-curated resistance/compensation pairs, once per pathway-DB
  configuration — this is the check that catches pipeline bugs before any
  biological conclusion is trusted, and the one place small-n results must
  be reported descriptively (e.g. "N of 10 known pairs ranked in top 5%"),
  not as a statistical test.
