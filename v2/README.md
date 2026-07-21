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
overlap by adding a network-diffusion (physics-inspired) propagation model
alongside topological distance, and by explicitly cross-validating pathway
evidence against independent, non-pathway-derived genetic evidence to avoid
circular reasoning — a rigor gap the original project's own scoring formula
did not address. The tool uses only public data: Open Targets Platform for
disease-target genetic associations, tractability, and safety, and Reactome
for pathway gene sets and functional interactions; no patient-specific data
is required. Given a target of interest and a disease, the pipeline finds
disease-relevant pathways, builds a gene-gene graph spanning the target's and
the disease's pathway neighborhoods, scores candidate genes by topological
and diffusion proximity to the target, re-weights by independent genetic
evidence, and annotates candidates with tractability and safety information
so results are actionable, not just statistically interesting. What's novel
here relative to the original project is the explicit orthogonality
discipline (excluding pathway-derived Open Targets evidence from the
validation step) and the addition of a literature-curated benchmark of known
clinical resistance mechanisms (e.g. EGFR→MET) to sanity-check the method
against ground truth before trusting its output. The expected deliverable is
a small web app: a researcher enters a target and a disease and gets back a
ranked, evidence-annotated list of candidate alternative targets. Known
limitations going in: pathway database coverage and curation bias, sparse
Open Targets safety annotation coverage for many genes, a benchmark set small
enough (5–15 cases) that validation results must be read as descriptive, not
statistically powered, and the deliberate exclusion of patient-specific
expression data from the default pipeline (available as an optional,
off-by-default filter). This tool is meant to complement, not replace,
genetic and clinical evidence in target prioritization.

## Pipeline overview

| Stage | Purpose | Addresses |
|---|---|---|
| 1. Ingest | Take target gene + disease (EFO ID) + pathway DB config; fetch/cache Open Targets and Reactome data. | Setup |
| 2. Disease-pathway discovery | GSEA (FDR-corrected) of disease-associated genes against pathway gene sets → disease-relevant pathway list. | H1 |
| 3. Pathway graph construction | Build a gene-gene graph from the union of disease-relevant and target-relevant pathway co-membership (+ optional Reactome functional-interaction edges); pin DB version, fix random seeds. | H1, H2 |
| 4. Candidate scoring | Swappable `--method`: topology score (co-membership / shortest path / branch-convergence) or diffusion score (personalized PageRank / random-walk-with-restart from the target node, via an existing library). | H1, H2, H3 |
| 5. Genetic-evidence integration | Re-weight/flag candidates using Open Targets evidence restricted to non-pathway-derived datatypes (e.g. `genetic_association`, `known_drug`), explicitly excluding `affected_pathway`/literature-pathway sources, to preserve orthogonality with Stages 2–4. | H4 |
| 6. Contextual annotation & filtering | Attach Open Targets tractability bucket and safety flags; compute a configurable composite score. Optional, off-by-default tissue-expression filter (GTEx) down-weights candidates not expressed in the disease-relevant tissue. | H5, H6, H10 |
| 7. Benchmark validation | Score a small literature-curated set of known resistance/compensation gene pairs (held out of Stage 2's seed genes); report descriptive rank/percentile recovery vs. a random background, not a significance test. | H8 |
| 8. Output/report | Ranked table (per-candidate scores, evidence trace, tractability/safety) plus a web UI: target + disease in, ranked list out. | Deliverable |

Hypotheses H7 (cross-pathway-DB consensus) and H9 (target-modality
generalization) were evaluated and explicitly dropped/deferred from hackathon
scope — see the brainstorm changelog below.

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

## Development plan

Sized for a small (2–4 person) mixed bio+CS team on laptop-scale compute.

**Day 1 — Core graph pipeline**
- Stage 1 (ingest) + Stage 2 (GSEA disease-pathway discovery, FDR-corrected)
- Stage 3 (pathway graph construction, versioned/seeded)
- Milestone: for a toy (target, disease) pair, produce a disease-relevant
  pathway list and a constructed gene-gene graph.

**Day 2 — Scoring and evidence integration**
- Stage 4 (topology + diffusion candidate scoring, both methods runnable)
- Stage 5 (genetic-evidence integration, non-pathway datatypes only)
- Stage 6 (tractability + safety annotation, composite score; expression
  filter only if time allows)
- Milestone: full ranked, annotated candidate list on a real example (e.g.
  rheumatic disease / PTGS2, reusing the original project's example as a
  known-reasonable sanity check).

**Day 3 — Validation, integration, demo**
- Stage 7 (benchmark validation against curated resistance cases)
- Stage 8 (output/report + minimal web UI)
- End-to-end integration test, documentation pass, demo rehearsal buffer.

## Testing steps

- **Per-stage sanity checks**: a tiny synthetic fixture (~15–20 genes, 3–5
  pathways with a hand-constructed, known-correct answer) run through each
  stage independently, checking exact expected output (e.g. Stage 3's graph
  has the expected edges, Stage 4's topology score matches a hand-computed
  value).
- **End-to-end integration run**: full pipeline on a real (target, disease)
  pair (e.g. PTGS2 / rheumatic disease, EFO_0005755) checked for a
  non-degenerate ranked output and sane runtime.
- **Ground-truth validation**: Stage 7's benchmark run against the
  literature-curated resistance/compensation pairs — this is the check that
  catches pipeline bugs before any biological conclusion is trusted, and the
  one place small-n results must be reported descriptively (e.g. "N of 10
  known pairs ranked in top 5%"), not as a statistical test.
