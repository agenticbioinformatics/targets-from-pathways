# Implementation Prompts

Ready-to-use prompts for an agentic coding assistant, organized by pipeline
stage and matched to the Day 1/2/3 plan in README.md. Each prompt is
self-contained — hand it to a fresh agent session with no other context
beyond this repo.

---

## Day 1

### Stage 0 — Shared schemas (prerequisite, ~45 min)

**Implement:**
> Write `schemas.py` defining the inter-stage data contracts for this
> pipeline as `pydantic` models (for the JSON manifest) and `pandera`
> DataFrameSchemas (for the tabular artifacts). Every later stage must
> validate its inputs on read and its outputs on write against these
> schemas — prose comments describing column names are not sufficient and
> have repeatedly failed on projects of this shape. Define at minimum:
>
> - `GenesSchema`: `gene_id` (str, Ensembl gene ID, unique, primary key),
>   `symbol` (str, HGNC approved), `synonyms` (list[str]), `biotype` (str).
> - `GeneSetsSchema`: `set_id` (str), `set_name` (str), `source_db` (str,
>   currently always `reactome`), `source_version` (str), `gene_id` (str,
>   FK to genes), `hierarchy_level` (int, nullable), `parent_id` (str,
>   nullable). One row per (set, gene) pair.
> - `InteractionsSchema`: `gene_a`, `gene_b` (str, FK to genes), `directed`
>   (bool), `sign` (int in {-1, 0, +1}), `source_db` (str),
>   `source_version` (str), `evidence_type` (str in {`curated`,
>   `predicted`}), `confidence` (float, nullable).
> - `OTAssociationsSchema`: `gene_id`, `disease_id`, `datatype_id`,
>   `datasource_id`, `score` (float in [0, 1]).
> - `Manifest`: run_id, git commit, created_at, seed, resolved target
>   (`{input, gene_id, symbol}`), disease (`{efo_id, name,
>   n_associated_genes}`), a list of sources (each with db name, version,
>   and per-file `{path, sha256, bytes}`), all CLI parameter values, all
>   output artifact paths with sha256, and the coverage/scale report
>   summaries.
>
> Note in the module docstring that `source_db` is retained as a column even
> though only `reactome` is emitted in this version, so that a second
> curation can be added later as one adapter without a downstream refactor.
> Do not add KEGG, WikiPathways or OmniPath support now.

**Test:**
> Add a pytest test asserting each schema rejects a deliberately malformed
> frame (wrong dtype, missing column, `sign` value of 2, a `gene_id` that
> looks like an HGNC symbol rather than an Ensembl ID) with a clear error
> naming the offending column.

### Stage 1 — Ingest & canonicalize

**Implement:**
> Write a standalone Python module `ingest.py` that resolves, canonicalizes
> and normalizes all source data for the pipeline, validating every output
> against `schemas.py` (Stage 0). CLI args: `--target` (gene symbol or
> Ensembl ID), `--disease` (EFO ID), `--data-dir` (source cache),
> `--out-dir` (run artifacts), `--pathway-db` (accepts only `reactome` for
> now; keep the flag so the value flows into the `source_db` column),
> `--no-download`, `--seed`, `--min-set-size` (default 10),
> `--max-set-size` (default 200), `--fi-curated-only` (default true),
> `--fi-min-score` (float, nullable), and `--allow-low-coverage`.
>
> **(a) Version-pinned acquisition.** Define a module-level registry mapping
> each required source file to a pinned URL, an expected sha256, and a size.
> Sources:
>
> - Open Targets Platform release **26.06** (released 2026-06-24, the latest
>   release at time of writing), specifically the `target` and
>   `association_by_datasource_indirect` parquet directories:
>   `ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.06/output/<table>/`.
>   **Pin the release number — never fetch from `.../platform/latest/`**,
>   which is a moving pointer and would silently change results between
>   runs. Record `26.06` in the manifest as `opentargets_version`.
> - Reactome: `ReactomePathways.gmt`, the Ensembl-to-pathway mapping, and
>   the functional-interaction file (`FIsInGene_*_with_annotations.txt`).
>   Pin Reactome to a numbered release directory too, never `current/`.
>
> Verify every versioned URL resolves before hardcoding it and record what
> you verified in a comment; the Reactome versioned paths in particular have
> not been confirmed and must not be assumed from the v1 README, which used
> `current/`. If a file is absent from `--data-dir`, download it and verify
> its sha256; under `--no-download`, print clear obtain instructions and exit
> non-zero instead. Never silently proceed on a checksum mismatch. Note in a
> comment that the OT parquet directories are multi-gigabyte and that
> `--no-download` exists for constrained environments.
>
> **(b) Identifier canonicalisation — the core responsibility of this
> stage.** The canonical internal identifier is the **Ensembl gene ID**.
> Build `genes.parquet` from the Open Targets `target` parquet, which is the
> ID authority here: using it means the pipeline's IDs agree with its
> primary evidence source by construction and no external mapping file is
> needed. Map every incoming Reactome symbol to a `gene_id` by trying the
> approved symbol first, then the synonym field, and record everything that
> fails to map. **Important:** verify whether the Reactome mapping file in
> use is gene-level or physical-entity-level before relying on it — v1 of
> this project used `Ensembl2Reactome_PE_All_Levels.txt` and leaked a
> physical-entity identifier (`GRB2-1`) into its published results. Use the
> gene-level file, or collapse PE identifiers to genes explicitly and say so
> in a comment.
>
> **(c) Disease subsetting.** Read the OT association parquets with
> `pyarrow` predicate pushdown filtered to `--disease`, and cache the result
> as a small `ot_disease_subset.parquet` so no later stage re-reads the full
> multi-gigabyte dump.
>
> **(d) Normalization.** Emit `gene_sets.parquet` and
> `interactions.parquet` conforming to the Stage 0 schemas. Parse direction
> and sign from the Reactome FI annotation columns; **inspect the real file
> header first and document the actual column names and the value
> vocabulary you found in a comment** rather than assuming a schema. Reactome
> FIs mix curated and predicted interactions, and predicted ones carry no
> reliable direction or sign, so honour `--fi-curated-only` /
> `--fi-min-score` here and record the filter in the manifest.
>
> **The interaction source is an open decision** — see "Open decision — what
> counts as pathway topology?" in README.md. The FI file assumed here is
> option B of five, chosen as a placeholder, not settled. Keep the FI parsing
> behind a single adapter function so swapping to a Reactome
> reaction-derived source (Pathway Commons SIF or BioPAX) is one function
> body, not a rewrite. Do not resolve this decision yourself.
>
> **(e) Reports.** Write `coverage_report.json` (per source: raw gene count,
> mapped count, percent mapped, dropped count, and up to 20 example dropped
> identifiers) and `scale_report.json` (gene-set size distribution before and
> after the size cap, count of sets retained, interaction counts broken down
> by sign and by `evidence_type`, and — critically — the projected Stage 3
> co-membership edge count as the sum of `C(n,2)` over retained sets). The
> projected edge count exists so an infeasible graph is caught now rather
> than on Day 2.
>
> **(f) Loud failures.** Exit non-zero with an actionable message when: the
> target cannot be resolved (offer fuzzy symbol suggestions); the disease has
> implausibly few `genetic_association` genes for GSEA to be meaningful; or
> any source maps below 90% of its genes, unless `--allow-low-coverage` is
> passed.
>
> Finally write `manifest.json`. Keep `ingest.py` a pure function of its
> inputs plus the pinned registry — no hidden global state.

**Test:**
> Add pytest tests for `ingest.py` using tiny fixture files under
> `tests/fixtures/` (a handful of fake genes, pathways and interactions).
> Cover exactly three things, and skip broader coverage until the pipeline
> runs end to end:
> 1. **Identifier mapping**: a clean HGNC symbol maps to the right Ensembl
>    ID; a deprecated symbol maps via the synonym field; an unmappable
>    symbol appears in the dropped list and is counted correctly in
>    `coverage_report.json`.
> 2. **Schema conformance**: every emitted artifact validates against its
>    Stage 0 schema, and a source file with a corrupt row fails loudly
>    rather than being silently skipped.
> 3. **Manifest determinism**: two runs over identical inputs produce
>    byte-identical manifests apart from `created_at`/`run_id`, and every
>    artifact sha256 recorded in the manifest matches the file on disk.

**Wire to next stage:**
> Make `manifest.json` plus the normalized parquet artifacts the sole input
> contract for `gsea_discovery.py` (Stage 2). No later stage may read a raw
> source file, resolve a file path from CLI args, or perform its own
> identifier mapping — all three belong to Stage 1 alone. Assert this by
> having each later stage take only `--manifest` and its own parameters.

### Stage 2 — Disease-pathway discovery (GSEA)

**Implement:**
> Write `gsea_discovery.py` that takes only the Stage 1 `--manifest` and
> reads the normalized artifacts it points to — never a raw source file.
> Retrieve disease-associated genes from `ot_disease_subset.parquet` using
> ONLY the `genetic_association` datatype (not the aggregated/overall score
> — this matters for avoiding circularity later), and run Gene Set
> Enrichment Analysis of those genes against the `gene_sets.parquet` sets
> (use an existing GSEA library such as `blitzgsea`, don't reimplement
> GSEA). Gene sets arrive already size-capped by Stage 1; before testing,
> additionally collapse near-duplicate sets (Jaccard > 0.7, keeping the
> smaller/more specific set) because `ReactomePathways.gmt` is an
> all-hierarchy-levels file in which ancestor and descendant sets are
> heavily redundant, which makes Benjamini-Hochberg anticonservative and
> inflates the significant-pathway count. Apply BH FDR correction across all
> tested pathways, and output a TSV of disease-relevant pathways with
> p-value, FDR, gene set, and `source_db`, filtered by user-specified
> `--pval-threshold` and `--fdr-threshold` CLI args (defaults 0.05 and 0.1).
> Also compute the set of pathways containing the `--target` gene from the
> manifest and include it in the output (a `contains_target` boolean
> column). Log the count of sets tested before and after collapsing.
>
> **Correcting for annotation bias is an open decision** — see "Open
> decision — correcting Stage 2 for annotation bias" in README.md. BH-FDR
> alone is option C of three, chosen as a placeholder. Well-studied genes
> carry both more Open Targets evidence and more Reactome annotation, so
> enrichment is inflated systematically and the significant-pathway count
> will likely be large (v1 returned 250). Structure the code so a
> per-pathway specificity score can be added as an extra column and filter
> without reshaping the output. Do not resolve this decision yourself.

**Test:**
> Add a pytest test for `gsea_discovery.py` using a small synthetic gene set
> fixture (e.g. 20 genes, 5 pathways, with a pathway you construct to be
> deliberately enriched) that verifies: the deliberately-enriched pathway is
> returned with p-value below threshold, FDR correction is actually being
> applied (test with a fixture where an individually-significant pathway
> becomes non-significant after correction), pathways not meeting threshold
> are excluded from output, and that two near-identical gene sets
> (Jaccard > 0.7) are collapsed to one before testing.

**Wire to next stage:**
> Modify `gsea_discovery.py` to write its pathway list as a stable,
> documented TSV schema, and add a `--benchmark-holdout-file` optional CLI
> arg: a list of genes (resolved to Ensembl IDs via Stage 1's `genes`
> table, not matched on symbol) to exclude from the disease-associated seed
> gene set before running GSEA. This will be used by Stage 7's benchmark
> validation to prevent circularity when known resistance-pair genes overlap
> with disease seed genes.

### Stage 3 — Pathway-based gene-gene graph construction

**Implement:**
> Write `build_graph.py` that takes Stage 2's disease-relevant pathway list
> and the target's own pathways, and constructs a **pathway-based gene-gene
> graph** — use this exact term in the module docstring, log messages, and
> variable names, not the generic "gene-gene graph". Build it as a
> `networkx.DiGraph`: nodes are genes keyed by Ensembl gene ID; edges come
> from two sources — (a) co-membership in at least one pathway from the
> union of (disease-relevant pathways) ∪ (target-containing pathways), added
> as a bidirectional pair of edges with `sign=0` (unknown effect); and (b)
> directional, signed edges read straight from Stage 1's
> `interactions.parquet` (already parsed, filtered and canonicalized — do
> not re-parse the Reactome FI file here). Tag every edge with `source_db`.
> Which curation those edges come from is an open decision (see README.md,
> "Open decision — what counts as pathway topology?"); this stage must work
> unchanged whichever option is chosen, since it only ever reads the
> normalized table.
>
> **Weight co-membership edges by `1 / (|pathway| - 1)`** (or another
> documented decreasing function of set size), summed over the pathways that
> induced them. Without this, a single 200-gene pathway contributes ~20,000
> edges and a handful of large Reactome sets determine the entire graph's
> topology, which is the failure mode behind v1's hub-dominated results.
> Compare the realized edge count against the projection in Stage 1's
> `scale_report.json` and log both.
>
> **Represent the graph as a `scipy.sparse` adjacency matrix plus a gene
> index for anything downstream that iterates over all edges** (Stage 5's
> diffusion in particular); keep `networkx` for construction, inspection and
> small test fixtures only. GraphML serialization of a million-edge graph is
> slow and large — serialize the sparse matrix plus the index instead, and
> say so in the docstring.
>
> Record the Reactome version and the random seed from the Stage 1 manifest
> as graph metadata.

**Test:**
> Add a pytest test with a small synthetic pathway list plus a small
> synthetic interactions table (a handful of directed, signed relations) and
> hand-verify the expected edge list; assert `build_graph.py` produces
> exactly those directed, signed edges with correct `source_db` tags, that
> co-membership edge weights match hand-computed `1 / (|pathway| - 1)`
> values for a fixture containing one small and one large pathway, and that
> a gene appearing in only one pathway alone (no shared pathway with any
> other gene) is an isolated node, not silently dropped.

**Wire to next stage:**
> Ensure `build_graph.py` outputs the target gene's node ID explicitly
> alongside the graph file (a small sidecar JSON with `target_node`,
> `graph_path`, `index_path` and `source_dbs`), since Stage 4's weighting
> and Stage 5's scoring methods (topology, RWR) all need to know which
> node is the source without re-deriving it.

---

## Day 2

### Stage 4 — Genetic-evidence weighting

**Implement:**
> Write `genetic_evidence_weights.py` that takes the Stage 3 graph (with
> its target node) and the Open Targets association data, and computes
> per-gene genetic-evidence scores restricted explicitly to non-pathway-
> derived datatypes (`genetic_association`, `known_drug` — exclude
> `affected_pathway` and any literature/pathway-derived datatype — note
> that OT 26.06 ships an `evidence_reactome` datasource feeding
> `affected_pathway`, so this exclusion is load-bearing, not theoretical:
> without it the pipeline would weight its graph with Reactome-derived
> evidence and then validate Reactome pathway structure against itself). Add a
> `--datatypes` CLI arg (comma-separated) defaulting to
> `genetic_association,known_drug` so the exclusion is explicit and
> overridable, with a code comment explaining why `affected_pathway` is
> excluded by default (orthogonality with Stages 2-3). Map these scores
> onto the graph as node weight attributes (`genetic_evidence_score`,
> default 0 for genes with no matching evidence), and derive edge weights
> by default via a `--edge-weight-mode` CLI arg (default `avg`, also
> supporting `product`, of the two endpoint node weights — both node and
> edge weighting are on by default, not optional). Output the updated
> graph in the same serialization format Stage 3 produces, plus a plain
> TSV of per-gene weights for inspection.

**Test:**
> Add a pytest test with a small fixture Open Targets association table
> containing a mix of datatypes for a few genes on a small synthetic graph;
> verify that node weights only reflect the allowed datatypes, that a gene
> with only `affected_pathway` evidence gets weight 0 by default, and that
> edge weights (when enabled) are computed correctly from endpoint node
> weights for a hand-picked edge.

**Wire to next stage:**
> Ensure `genetic_evidence_weights.py` writes the weighted graph in exactly
> the same format Stage 3 produces, so `score_candidates.py` (Stage 5) can
> load either the unweighted Stage 3 graph or this weighted graph
> interchangeably via the same `--graph` argument. Document the node/edge
> weight attribute names as a comment at the top of the script.

### Stage 5 — Candidate scoring (topology, RWR)

**Implement:**
> Write `score_candidates.py` that loads a graph (optionally weighted, from
> Stage 3 or Stage 4) and target node, and supports `--method` as a
> comma-separated list drawn from `{topology, rwr}`, **defaulting to
> `topology` only** — RWR is opt-in (e.g. `--method topology,rwr` or
> `--method rwr`), never run unless explicitly requested, since it's
> materially more expensive (a converged stationary distribution) than the
> deterministic topology method. Topology method: shortest-path distance to
> the target node and co-membership count, using edge weights if present.
> RWR method: personalized PageRank / random-walk-with-restart (via
> `networkx.pagerank` or an existing RWR library) personalized on the
> target node, using node/edge weights if present — do not hand-roll this
> algorithm. Use the fixed seed from the graph metadata for any stochastic
> step. Output a ranked TSV with a `gene` column and a `topology_score`
> column always present, plus an `rwr_score` column only when that method
> was also requested.
>
> (Backtrack-free walk / non-backtracking RWR via GBA-centrality was
> considered as a third method and dropped — see README.md plan update 5:
> it needs a compiled C extension and its own submodule, out of proportion
> to a hackathon prototype, and the core hypotheses don't require it.)

**Test:**
> Add a pytest test on a small synthetic graph containing both a short
> cycle and a longer chain where you can hand-compute expected
> shortest-path distances; assert that a default (no `--method`) run
> produces only `topology_score`, matching the hand-computed values
> exactly, with no `rwr_score` column at all; then assert that running with
> `--method topology,rwr` produces both columns, that `rwr_score` is a
> valid probability distribution (non-negative, sums to ~1), and that the
> target's own graph neighbours outrank far-away chain nodes.

**Wire to next stage:**
> Confirm `score_candidates.py`'s output TSV uses a `gene` column that
> matches the gene identifier convention used by the Open Targets data
> (document whether it's HGNC symbol or Ensembl ID, and add an ID-mapping
> step if they differ) so Stage 6 can join on it directly without silent
> mismatches.

### Stage 6 — Contextual annotation & filtering

**Implement:**
> Write `annotate_context.py` that takes the running candidate TSV and
> joins in, per gene, an Open Targets tractability bucket and a safety flag
> (from Open Targets safety data), computing a `composite_score` as a
> configurable weighted sum (`--weights` CLI arg, e.g.
> `topology:0.3,rwr:0.3,genetic_evidence:0.2,tractability:0.1,safety:0.1`,
> with sane documented defaults). Document that safety data has sparse
> coverage and genes with no safety annotation should not be penalized as
> if they were confirmed safe (flag them as `safety: unknown`, not
> `safety: safe`). Do not incorporate any tissue-expression data — this
> pipeline is intentionally pathway- and association-evidence-based only,
> with no expression filter of any kind.

**Test:**
> Add a pytest test with a small fixture tractability/safety table covering
> a few genes plus at least one gene with no safety annotation at all;
> verify the composite score formula matches a hand-computed value for a
> known input, and that the no-annotation gene is marked `unknown` rather
> than defaulting to a numeric safety score.

**Wire to next stage:**
> Confirm `annotate_context.py` writes the final per-candidate evidence
> trace (which pathways, which datatypes, which tractability bucket
> contributed) alongside the composite score, since Stage 8's report needs
> this trace for interpretability, not just the final number.

---

## Day 3

### Stage 7 — Benchmark validation

**Implement:**
> Write `benchmark_validate.py` that takes a small hand-curated TSV of known
> resistance/compensation gene pairs (columns: `original_target`,
> `alternative_target`, `disease_or_context`, `source_citation` — create
> this file **must be curated by a human, not generated here** — it is the
> single file the entire validation rests on, and every `source_citation`
> needs a real, checked PMID. Assume it already exists at
> `benchmark/resistance_pairs.tsv`; if it does not, stop and say so rather
> than inventing cases or citations). Run the full pipeline (Stages 1-6) for
> each known pair with `--benchmark-holdout-file` set to exclude the
> `alternative_target` gene from Stage 2's seed genes, memoizing the disease
> GSEA across pairs that share a disease so the benchmark does not re-run
> enrichment unnecessarily. Record the rank and percentile of
> `alternative_target` in the resulting candidate list. Print a clear
> disclaimer that with n<20 cases results must be reported descriptively
> (e.g. "X of N known pairs ranked in top 5%"), not as a significance test —
> do not compute or report a p-value against this benchmark.

**Test:**
> Add a pytest test using a tiny synthetic pipeline setup (small graph,
> small association table) with one deliberately "easy" known pair (short
> graph distance) and confirm `benchmark_validate.py` reports it near the
> top of the ranking with the correct percentile calculation.

**Wire to next stage:**
> Ensure `benchmark_validate.py` writes its summary as both a TSV and a
> short human-readable text block, so Stage 8's report can embed the
> validation summary directly without re-parsing raw output.

### Stage 8 — Output/report + web UI

**Implement:**
> Write a small pipeline orchestrator `run_pipeline.py` that chains Stages
> 1-6 given `--target` and `--disease`, and a minimal web app (e.g. Flask or
> Streamlit) with a form for target + disease, running `run_pipeline.py` and
> displaying the ranked candidate table (composite score, topology, RWR,
> genetic evidence, tractability, safety) plus, for each candidate on
> click/expand, the evidence trace (which shared pathways, which datatypes,
> which interactions contributed). Display genes by HGNC symbol while
> keeping Ensembl gene IDs as the internal key — resolve symbols via Stage
> 1's `genes` table at render time only. Include the Stage 7 benchmark
> summary as a static "validation" panel/page in the app, not recomputed per
> request.

**Test:**
> Write an integration test that runs `run_pipeline.py` end-to-end on the
> PTGS2 / rheumatic disease (EFO_0005755) example from the original
> project's README, and asserts the output is non-empty, PTGS2 itself
> appears in the pathway data, and the pipeline completes within a
> reasonable time budget (e.g. under 5 minutes) — this is a smoke test, not
> a correctness check on rankings.

---

## Final integration, documentation, demo

**Integration testing:**
> Write an end-to-end pytest that runs `ingest.py` → `gsea_discovery.py` →
> `build_graph.py` → `genetic_evidence_weights.py` → `score_candidates.py` →
> `annotate_context.py` as a single chained test on the small synthetic
> fixtures used in the individual stage tests, asserting that every
> hand-off validates against its `schemas.py` contract (no missing or
> renamed columns between stages, Ensembl gene IDs and `source_db`
> attribution preserved end to end) and that the final output contains a
> `composite_score` for every input gene. This is the highest-value test in
> the suite — it is what catches a schema break on Day 3 morning.

**Documentation:**
> Update `README.md` in this directory with a "How to run" section showing
> the exact CLI invocation chain for a full run (mirroring the original
> project's README style), and add a `--help` docstring to every stage
> script's argparse setup so the CLI is self-documenting.

**Demo/web app polish:**
> Add basic input validation and error messages to the web app from Stage 8
> (e.g. clear message if the EFO ID or target gene isn't found in the
> loaded data, instead of a stack trace), and a loading indicator while
> `run_pipeline.py` executes, since a full run may take up to a few minutes.
