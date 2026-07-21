# Implementation Prompts

Ready-to-use prompts for an agentic coding assistant, organized by pipeline
stage and matched to the Day 1/2/3 plan in README.md. Each prompt is
self-contained — hand it to a fresh agent session with no other context
beyond this repo.

---

## Day 1

### Stage 1 — Ingest

**Implement:**
> Write a standalone Python module `ingest.py` for a target-discovery
> pipeline. It should accept CLI args `--target` (gene symbol or Ensembl ID),
> `--disease` (EFO ID), `--data-dir` (path to a local cache directory), and
> `--pathway-db` (choice: `reactome`, default `reactome`). If the required
> Open Targets association/target parquet files and Reactome pathway files
> (ReactomePathways.gmt, Ensembl2Reactome_PE_All_Levels.txt) are not already
> present in `--data-dir`, print clear instructions for downloading them
> (reuse the download commands from `../README.md`'s "Data" section) rather
> than downloading automatically. Output: a config/manifest object (as JSON)
> recording the resolved file paths, pathway DB version/date, and a fixed
> random seed, to be consumed by later stages. Keep it a pure function of its
> inputs — no hidden global state.

**Test:**
> Add a small pytest test for `ingest.py` using a tiny fixture directory
> (create fixture files with a handful of fake genes/pathways under
> `tests/fixtures/`) that checks: (1) the manifest is produced with correct
> resolved paths when files exist, (2) a clear, non-crashing error message is
> printed when a required file is missing, listing exactly which file and how
> to get it.

**Wire to next stage:**
> Update `ingest.py` so the manifest JSON it produces is the sole input
> contract for `gsea_discovery.py` (Stage 2) — no other stage should read
> raw file paths directly from CLI args, only from the manifest. Document the
> manifest schema (field names and types) as a comment at the top of
> `ingest.py`.

### Stage 2 — Disease-pathway discovery (GSEA)

**Implement:**
> Write `gsea_discovery.py` that takes the Stage 1 manifest and a
> `--disease` EFO ID, retrieves disease-associated genes from the Open
> Targets association data using ONLY the `genetic_association` datatype
> (not the aggregated/overall score — this matters for avoiding circularity
> later), runs Gene Set Enrichment Analysis of those genes against the
> Reactome pathway gene sets (use an existing GSEA library such as
> `blitzgsea`, don't reimplement GSEA), applies Benjamini-Hochberg FDR
> correction across all tested pathways, and outputs a TSV of disease-
> relevant pathways with p-value, FDR, and gene set, filtered by
> user-specified `--pval-threshold` and `--fdr-threshold` CLI args (defaults
> 0.05 and 0.1). Also compute the set of pathways containing the
> `--target` gene from Stage 1's manifest and include it in the output
> (a `contains_target` boolean column).

**Test:**
> Add a pytest test for `gsea_discovery.py` using a small synthetic gene set
> fixture (e.g. 20 genes, 5 pathways, with a pathway you construct to be
> deliberately enriched) that verifies: the deliberately-enriched pathway is
> returned with p-value below threshold, FDR correction is actually being
> applied (test with a fixture where an individually-significant pathway
> becomes non-significant after correction), and pathways not meeting
> threshold are excluded from output.

**Wire to next stage:**
> Modify `gsea_discovery.py` to write its pathway list as a stable,
> documented TSV schema, and add a `--benchmark-holdout-file` optional CLI
> arg: a list of gene symbols to exclude from the disease-associated seed
> gene set before running GSEA. This will be used by Stage 7's benchmark
> validation to prevent circularity when known resistance-pair genes overlap
> with disease seed genes.

### Stage 3 — Pathway graph construction

**Implement:**
> Write `build_graph.py` that takes Stage 2's disease-relevant pathway list
> and the target's own pathways, and constructs a gene-gene graph (using
> `networkx`) where nodes are genes and edges connect genes that co-occur in
> at least one pathway from the union of (disease-relevant pathways) ∪
> (target-containing pathways). Add an optional `--functional-interactions`
> flag to also add edges from a Reactome functional-interaction file (same
> format as `FIsInGene_04142025_with_annotations.txt` in the original
> project) if provided. Record the pathway DB version/date (from the Stage 1
> manifest) and the random seed as graph metadata attributes. Serialize the
> graph to a `.graphml` or pickled networkx format.

**Test:**
> Add a pytest test with a small synthetic pathway list (3-4 pathways, each
> with 4-6 genes, some genes shared across pathways) and hand-verify the
> expected edge list; assert `build_graph.py` produces exactly those edges,
> and that a gene appearing in only one pathway alone (no shared pathway
> with any other gene) is an isolated node, not silently dropped.

**Wire to next stage:**
> Ensure `build_graph.py` outputs the target gene's node ID explicitly
> alongside the graph file (e.g. a small sidecar JSON with `target_node` and
> `graph_path`), since Stage 4's weighting and Stage 5's scoring methods
> (topology, RWR, BFW) all need to know which node is the source, without
> re-deriving it.

---

## Day 2

### Stage 4 — Genetic-evidence weighting

**Implement:**
> Write `genetic_evidence_weights.py` that takes the Stage 3 graph (with
> its target node) and the Open Targets association data, and computes
> per-gene genetic-evidence scores restricted explicitly to non-pathway-
> derived datatypes (`genetic_association`, `known_drug` — exclude
> `affected_pathway` and any literature/pathway-derived datatype). Add a
> `--datatypes` CLI arg (comma-separated) defaulting to
> `genetic_association,known_drug` so the exclusion is explicit and
> overridable, with a code comment explaining why `affected_pathway` is
> excluded by default (orthogonality with Stages 2-3). Map these scores
> onto the graph as node weight attributes (`genetic_evidence_score`,
> default 0 for genes with no matching evidence), and optionally derive
> edge weights via a `--edge-weight-mode` CLI arg (`none` default, or
> `avg`/`product` of the two endpoint node weights). Output the updated
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

### Stage 5 — Candidate scoring (topology, RWR, BFW)

**Implement:**
> Write `score_candidates.py` that loads a graph (optionally weighted, from
> Stage 3 or Stage 4) and target node, and supports `--method topology`,
> `--method rwr`, and `--method bfw` (default: run all three, output as
> separate columns). Topology method: shortest-path distance to the target
> node and co-membership count, using edge weights if present. RWR method:
> personalized PageRank / random-walk-with-restart (via `networkx.pagerank`
> or an existing RWR library) personalized on the target node, using
> node/edge weights if present — do not hand-roll this algorithm. BFW
> method: backtrack-free walk (non-backtracking random-walk-with-restart)
> — use the existing [GBA-centrality](https://github.com/jedrzejkubica/GBA-centrality)
> tool (add it as a dependency, e.g. git submodule or local pip-installable
> package) rather than reimplementing non-backtracking walk logic; check
> its README for the expected graph input format and adapt Stage 3/4's
> graph export if needed. Use the fixed seed from the graph metadata for
> any stochastic step. Output a ranked TSV: gene, topology_score, rwr_score,
> bfw_score.

**Test:**
> Add a pytest test on a small synthetic graph containing both a short
> cycle and a longer chain (backtrack-free walks differ most from RWR when
> short cycles are present) where you can hand-compute expected
> shortest-path distances; assert `score_candidates.py`'s topology output
> matches exactly, that `rwr_score` and `bfw_score` are both produced, and
> that they differ on the cyclic portion of the fixture in the expected
> direction (BFW should not artificially inflate a node reachable only by
> immediate back-and-forth on a single edge).

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
> `topology:0.2,rwr:0.2,bfw:0.2,genetic_evidence:0.2,tractability:0.1,safety:0.1`,
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
> this file yourself with 5-15 literature-documented cases, e.g.
> EGFR/MET, BRAF/CRAF, ESR1/PIK3CA-AKT1-MTOR pathway crosstalk; look up real
> citations, don't fabricate them) and, for each pair, runs the full
> pipeline (Stages 1-6) with `--benchmark-holdout-file` set to exclude the
> `alternative_target` gene from Stage 2's seed genes, then records the rank
> and percentile of `alternative_target` in the resulting candidate list.
> Output a summary table and print a clear disclaimer that with n<20 cases
> results should be reported descriptively (e.g. "X of N known pairs ranked
> in top 5%"), not as a significance test — do not compute or report a
> p-value against this benchmark.

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
> BFW, genetic evidence, tractability, safety) plus, for each
> candidate on click/expand, the evidence trace (which shared pathways,
> which datatypes). Include the Stage 7 benchmark summary as a static
> "validation" panel/page in the app, not recomputed per request.

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
> `build_graph.py` → `score_candidates.py` → `genetic_evidence.py` →
> `annotate_context.py` as a single chained test on the small synthetic
> fixtures used in the individual stage tests, asserting the schema
> contract holds at every hand-off (no missing/renamed columns between
> stages) and the final output contains a `composite_score` for every input
> gene.

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
