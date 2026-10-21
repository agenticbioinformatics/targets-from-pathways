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

There is no target gene: the disease's genes with an Open Targets
`genetic_association` score > 0 are the **seed genes** (Stage 2 writes them to
`disease_seed_genes.tsv`). GSEA of the seeds finds the disease-relevant
pathways (Stage 2), the graph is the Reactome functional interactions between
the genes of those pathways (Stage 3), and Stage 5 propagates from the seeds
to rank every non-seed graph gene.

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

`python3 -m pytest tests/test_run_directory.py` (9 tests) checks run-directory
integrity on copies of the example run: a copied run reads its own Stage 1
files even when the manifest lists the original's, an artifact changed after
Stage 1 is refused, a benchmark holdout gene's own evidence never feeds Stage
4 (and a re-run without a holdout clears it), Stage 4 leaves Stage 3's
unweighted graph unchanged without user weights (and applies and resets them
correctly with them), genes with a zero genetic-association score are
not seeds, and a disease with too few seeds for GSEA, or with no seed on a
disease-relevant pathway, stops Stage 3 with a message.

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
`disease_seed_genes.tsv` (`gene_id`, `score`: the genes with a
`genetic_association` score > 0 — the GSEA ranking signature — minus any
benchmark holdout) and `benchmark_holdout_genes.tsv` (`gene_id`: the
resolved holdout, header only without `--benchmark-holdout-file`). Stage 5
propagates from `disease_seed_genes.tsv` (Stage 3 only checks that some seeds
are in the graph); Stage 4 reads the holdout.

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
`benchmark_holdout_genes.tsv`), so it does not show as their Stage 6
`genetic_evidence_score`. Genes are listed one per line
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
longer a seed (17 seeds), and since the Stage 3 graph holds only the genes of
the GSEA pathways, G17–G21 leave the graph with `R-HSA-300`: G17 cannot be
ranked. A held-out gene can come back as a candidate only if the remaining
seeds still detect one of its pathways.

**Diseases with few genetic-association genes** (e.g. a Mendelian disease
run with Stage 1's `--allow-low-coverage`): GSEA cannot always run. Stage 2
then logs a warning saying why and writes an empty `disease_pathways.tsv`;
with no disease-relevant pathway there is no graph, so Stage 3 stops with a
message. This happens when:
- there are fewer than 5 seed genes (blitzgsea's minimum overlap per gene
  set);
- no gene set holds 5 of them;
- blitzgsea's calibration fails. It breaks on small signatures, including
  whenever the seed count equals one of its calibration set sizes (up to
  about 230 genes).

If no seed gene is on a disease-relevant pathway, Stage 3 also stops with a
message: there is nothing to propagate from. Spermatogenic failure 5
(`MONDO_0009461`) is an example: its one seed gene, AURKC, is in no Reactome
pathway (and one seed is too few for GSEA anyway).

`python3 -m pytest tests/test_gsea_discovery.py` runs the unit-level GSEA
tests (`tests/test_gsea_discovery.py`, 12 tests): GSEA skipped with a
warning and an empty result when there are too few seed genes, when
blitzgsea's calibration fails, or when no gene set can be tested; a deliberately-enriched
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
all genes of the Stage 2 disease pathways (the pathways GSEA found), and its
edges are the Reactome functional interactions between those genes; sharing
a pathway adds no edge. The graph is **unweighted**: every edge has weight 1
in `graph_weight.npz` (interaction confidence is not used as a weight, and
there is no genetic-evidence weighting). A Reactome interaction with no
direction (`-` in the FI file) is added in both directions; a directed one
only in its own direction. The seed genes are read only to check that some
are in the graph. Writes `graph_weight.npz`, `graph_sign.npz`,
`graph_gene_index.json`, and `graph_metadata.json` next to the manifest —
not GraphML; see the module docstring for why.

Stage 3 stops with a message if `disease_pathways.tsv` is empty (GSEA found
no pathway) or if no seed gene is on a disease pathway.

Against the checked-in example (regenerate with
`python3 example_data/build_gsea_example.py` if needed):

```
python3 pipeline/stage3_build_graph.py --manifest example_data/example_run/manifest.json
```

**Expected output** (log lines):

```
INFO build_graph: Graph: 19 nodes (genes on the 3 disease-relevant pathways), 10 edges (10 interaction edges added, 3 out of scope).
INFO build_graph: Seed genes in the graph: 16/18.
```

The two seeds left out (`G12`, `G13`) are only on `R-HSA-100` (dropped by
Stage 2's collapsing) and `R-HSA-400` (not significant). The genes of
`R-HSA-500` (too few seeds for GSEA) and `R-HSA-600` (no seed) are not
nodes either, so the three interaction rows touching only those genes
(`G12 -> G1`, `G2 -> N25`, `N29 -> N22`) are dropped as out of scope — proof
that interaction edges are scoped to the genes of the GSEA pathways, not
added genome-wide. `graph_metadata.json` records `n_disease_pathways`,
`n_nodes`, `n_edges`, `n_seed_genes_in_graph` and `n_seed_genes_total`.

`python3 -m pytest tests/test_build_graph.py` runs the unit-level graph
tests (4 tests): the edges are exactly the in-scope interactions (with
sign, confidence, `evidence_type` and `source_db`; no co-membership edges),
the nodes are exactly the genes of the GSEA pathways (a seed-gene pathway
that GSEA did not detect stays out, and a gene without interactions
survives as an isolated node), the serialized weight matrix is 1 on every
edge, and an undirected interaction comes out symmetric.

## Running Stage 4 (`pipeline/stage4_genetic_evidence_weights.py`)

**Setup:** same as Stage 3, plus `numpy` (already required transitively).

Stage 4 does two independent things: it **optionally weights the graph with
your own node and/or edge weights**, and it scores each graph gene's genetic
evidence for Stage 6's annotation.

**Run command**, after Stage 3 has produced the graph artifacts:

```
python3 pipeline/stage4_genetic_evidence_weights.py --manifest out/manifest.json \
    [--node-weights node_weights.tsv] [--edge-weights edge_weights.tsv]
```

Reads `ot_disease_subset.parquet` via the manifest and Stage 3's graph
artifacts alongside it by default (override with `--graph-dir`/`--out-dir`).
Writes **exactly Stage 3's four filenames** (`graph_weight.npz` is
replaced; `graph_sign.npz` and `graph_gene_index.json` are copied through;
`graph_metadata.json` is extended), so Stage 5 reads a Stage 3 or a Stage 4
directory the same way, plus `gene_weights.tsv`.

### Your own node and edge weights (optional)

Stage 3's graph is unweighted: every edge has weight 1. If you have your own
prior knowledge — e.g. expression in a relevant tissue, a druggability
score, or interaction confidences from another database — you can turn it
into weights that Stage 5's propagation (BFWalk and RWR) then follows: a
heavier edge carries more of the seeds' signal. Both files are optional and
can be given alone or together; without either, the graph stays unweighted.

**`--node-weights`**: a tab-separated file with a header and two columns:

```
gene_id	weight
ENSG00000105974	2.0
ENSG00000141510	0.5
```

Each graph edge `u → v` gets the **mean** of `u`'s and `v`'s weights
(`--node-weight-mode product` multiplies them instead, which penalises an
edge more strongly when either gene has a low weight). Genes not in the file
weigh `--default-node-weight` (default 1).

**`--edge-weights`**: a tab-separated file with a header and three columns:

```
gene_a	gene_b	weight
ENSG00000105974	ENSG00000141510	0.8
```

A row weights the graph edge(s) between `gene_a` and `gene_b` in **whichever
directions the graph has** them (an undirected Reactome interaction is two
directed edges, and both get the weight), so the order of the two genes does
not matter. Listing the same pair twice, in either order, is an error. Edges
not in the file weigh `--default-edge-weight` (default 1).

**With both files**, an edge's weight is the product of its edge weight and
its combined node weight. For example, with node weights `A = 0.5`, `B = 1`
and an edge weight of 0.4 for the pair (A, B), the edge `A → B` weighs
`0.4 × (0.5 + 1) / 2 = 0.3`.

Rules that apply to both files:
- Genes are **Ensembl gene IDs** (`ENSG` + 11 digits), as everywhere else in
  the pipeline; symbols are rejected.
- Every weight must be a finite number **> 0**. A weight of 0 would cut an
  edge out of propagation; to remove an edge, leave its interaction out of
  the graph instead.
- Weights are relative, not bounded: RWR uses them as they are, and BFWalk
  divides them all by the largest one (it needs weights in `]0, 1]`). Mind
  the defaults: if your weights are all below 1, every gene or edge missing
  from the file (weight 1) outweighs the ones you listed — set
  `--default-node-weight` / `--default-edge-weight` to match your scale.
- Stage 4 **never adds or removes an edge**. Rows for genes that are not
  graph nodes, and edge rows for gene pairs that are not graph edges, are
  skipped, and the log says how many (e.g. `Edge weights: 1 gene pairs
  weighted, 3 not graph edges (skipped).`).
- The weights are always applied to Stage 3's structure, starting from 1 on
  every edge. Re-running Stage 4 with other files, or with none, replaces the
  previous weighting rather than compounding it; running it without files
  restores the unweighted graph.

`graph_metadata.json` records what was used under `stage4_user_weights`
(`weighted`, the two file paths, the mode and defaults, and how many rows
were used or skipped), and the Stage 7 report header says whether the graph
was weighted. `run_pipeline.py` passes `--node-weights` / `--edge-weights`
through to Stage 4; for `--node-weight-mode` and the defaults, run Stage 4
on its own and then Stages 5-6.

### Genetic-evidence gene scores

Independently of the weights, Stage 4 scores each graph gene with its
strongest Open Targets `genetic_association` evidence for the disease, for
Stage 6's `genetic_evidence_score` column only — **these scores are never
used as weights**. Only `genetic_association` evidence counts; every other
Open Targets datatype is ignored, and there is no option to add one:
- `affected_pathway` would be circular: OT 26.06 maps its own `reactome`
  evidence datasource to it (see `pipeline/stage1_ingest.py`'s
  `DATASOURCE_TO_DATATYPE`), and the graph is built from Reactome.
- `known_drug` would favour genes that are already drug targets for the
  disease.

The scores go into `graph_metadata.json` as a `genetic_evidence_score`
gene->score mapping, which Stage 6 reads, and into a human-facing
`gene_weights.tsv`. Genes listed in Stage 2's `benchmark_holdout_genes.tsv`
(next to the manifest; override with `--holdout-genes`) score 0.0 whatever
their Open Targets rows say, and `graph_metadata.json` records them under
`stage4_genetic_evidence.holdout_genes`. With Stage 2's seeds every candidate
scores 0 (a gene with a score > 0 is a seed); the column is informative only
when Stage 5 is given other seeds with `--seed-genes`.

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

Every non-seed gene scores 0: seeds are exactly the genes with a
`genetic_association` score > 0. `ENSG...022`'s `known_drug` row is ignored.
The log also says `No --node-weights / --edge-weights: the graph stays
unweighted.` To try the weights on the example, give `ENSG00000000022` a
low node weight:

```
printf 'gene_id\tweight\nENSG00000000022\t0.1\n' > /tmp/node_weights.tsv
python3 pipeline/run_pipeline.py --from-manifest example_data/example_run/manifest.json \
    --node-weights /tmp/node_weights.tsv
```

`ENSG00000000022`'s edges (all with seeds) now weigh `(1 + 0.1) / 2 = 0.55`,
and its BFWalk score drops from 0.5 to about 0.27, below `ENSG00000000024`
(0.5). Run the command again without `--node-weights` to restore the
checked-in result.

`python3 -m pytest tests/test_genetic_evidence_weights.py` runs the
unit-level tests (13 tests): only `genetic_association` rows count (the
strongest one per gene), and a gene with only `affected_pathway` or
`known_drug` evidence scores 0; without user weights every edge is reset to
1; node weights combine by mean or product, with the default for missing
genes; edge weights match either gene order, skip non-edges and never add
an edge; node and edge weights multiply; a pair listed twice, and invalid
files (weight ≤ 0, a symbol instead of an Ensembl ID, a wrong column, a
duplicate gene), exit non-zero. `tests/test_run_directory.py` also runs
Stage 4 on the example with both files and checks that a re-run without
them restores the unweighted graph.

## Running Stage 5 (`pipeline/stage5_score_candidates.py`)

**Setup:** `pip install networkx scipy numpy pandas pandera pydantic` (or `pip install -r requirements.txt`).

**Run command**, against a Stage 3 or Stage 4 graph directory (the two
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

`--method` is a comma-separated, non-empty subset of `{bfwalk, rwr}`,
**default `bfwalk`**; anything else (including an empty list) exits
non-zero. Each requested method writes its own column (`bfwalk_score`,
`rwr_score`). With both (`--method bfwalk,rwr`) the output is ranked by
`bfwalk_score`. Both start from the seed set:

- `rwr` = personalised PageRank restarting uniformly on the seeds.
  `--restart-prob` (default 0.5; PageRank `alpha = 1 - restart_prob`).
  Every edge of Stage 3's unweighted graph has weight 1; Stage 4's user
  weights, if given, are used as they are.

Hub bias, measured on inflammatory bowel disease (`--disease MONDO_0005265`, Open Targets 26.06 +
Reactome 97, 9,127 candidates from 502 seeds) as the rank correlation of each score with
the candidate's number of graph neighbours and of Reactome pathways. These numbers come from an earlier
graph (genetic-evidence edge weights, co-membership edges, and the pathways of every seed gene) and have
not been re-measured on the current unweighted GSEA-pathway graph:

| Method | Spearman ρ with node degree | Spearman ρ with no. of Reactome pathways |
|---|---|---|
| `rwr` | 0.74 | 0.62 |
| `bfwalk` | 0.13 | 0.13 |

RWR's top IBD candidates mix existing drug targets (PPARG, TNF, ALOX5) with
ubiquitin hubs (RPS27A, UBA52, UBB, UBC). BFWalk's top 100 have a median
degree of 48, against 174 for all candidates, and its top hits are
small metabolic modules (PDK2, PDK3, MPC1). There is no topology score
any more: its (seed neighbours / distance) formula had ρ = 0.85 with node
degree and ranked ubiquitin and proteasome genes first.

**Which to use:** use `bfwalk` (the default) to find candidates beyond well-connected hub genes, and
add `rwr` when you prefer well-connected, well-studied genes close to the seeds, e.g. to recover
established targets.

`bfwalk` runs [BFWalk](https://github.com/jedrzejkubica/BFWalk)
(backtrack-free network propagation with in-degree normalisation) with the
seed genes as its seeds, on the same directed graph. BFWalk needs edge weights in
`]0, 1]`, so all weights are divided by the largest one (a no-op on the
unweighted graph); zero-weight edges from a hand-made graph would be dropped. Knobs: `--bfwalk-alpha` (attenuation coefficient, default 0.5)
and `--bfwalk-threads` (default 0 = all cores). BFWalk is not on PyPI yet:
`requirements.txt` installs it from TestPyPI (see the top-level README); it
will be available via `pip install bfwalk` shortly. Until it is installed,
run Stage 5 with `--method rwr`. Stages 6 and 7 carry whichever of
`bfwalk_score` / `rwr_score` was produced through, and keep its ranking.

Against the checked-in example (after Stages 1–4 have run):

```
python3 pipeline/stage5_score_candidates.py --graph-dir example_data/example_run
```

**Expected output** (`example_data/example_run/candidate_scores.tsv`, all
3 candidates; columns per `schemas.CandidateScoresSchema`, ranked by
`bfwalk_score`):

```
gene             bfwalk_score
ENSG00000000022  0.5
ENSG00000000024  0.5
ENSG00000000023  0.25
```

No candidate has genetic evidence of its own, so the ranking comes from
network position alone: `N22` and `N24` have an edge from a seed, `N23` is
reached only through `N24`. The 16 seeds in the graph are not in the table.

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
tests (23 tests) on a small two-seed cycle+chain graph:
- RWR restarting on the seeds, reconstructed as a full stationary
  distribution (sums to 1), with an unreachable node scoring 0.
- BFWalk following edge direction and ignoring zero-weight edges.
- The output-column contract:
  - the default `run()` writes `bfwalk_score` only, and `--method rwr`
    writes `rwr_score` only;
  - both methods write both columns, ranked by `bfwalk_score` whatever
    order they are given in;
  - seed rows are excluded;
  - an empty, unknown or `topology` `--method` exits non-zero.
- `--seed-genes` overriding the default seed file, including a seed list
  without a `score` column.
- Seeds outside the graph dropped; no seed in the graph, or a missing seed
  file, exiting non-zero.
- An all-seed graph writing an empty ranking.
- `--restart-prob` outside ]0, 1[ exiting non-zero.

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

Per candidate it adds the columns below. They are for reading only: there
is no composite score, and candidates keep Stage 5's ranking (by
`bfwalk_score`, else `rwr_score`).

- **`genetic_evidence_score`** — Stage 4's score for the gene (0.0 for a
  benchmark holdout gene); left out when Stage 4 did not run.
- **`tractability`** — `clinical` (an Open Targets `value==true` flag for a
  clinical stage), `discovery` (some other true flag), or `unknown` (no
  true flag, or the gene is absent from the `target` parquet).
- **`safety`** — `has_liabilities` (one or more `safetyLiabilities`
  entries) or `unknown`. **There is no `safe` value**: Open Targets records
  known liabilities or says nothing, and its safety coverage is sparse, so
  a gene with no annotation is `unknown`, never presented as confirmed
  safe. `n_safety_liabilities` carries the count.

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

**Expected output** (`example_data/example_run/candidates_annotated.tsv`; columns per `schemas.AnnotatedCandidatesSchema`, in Stage 5's
order by `bfwalk_score`):

```
gene             bfwalk_score  genetic_evidence_score  tractability  safety           n_safety_liabilities
ENSG00000000022  0.5           0.0                     clinical      has_liabilities  1
ENSG00000000024  0.5           0.0                     discovery     unknown          0
ENSG00000000023  0.25          0.0                     unknown       unknown          0
```

The order is exactly `candidate_scores.tsv`'s. `ENSG...022`'s safety
liability is shown, but does not move it. `genetic_evidence_score` is 0 for every candidate: candidates
are non-seed genes, and Stage 4 counts only `genetic_association`
evidence, which is what makes a gene a seed. Genes absent from the
synthetic `ot_target_subset.parquet` are `unknown`/`unknown`.

`python3 -m pytest tests/test_annotate_context.py` runs the unit-level
tests (7 tests):
- tractability bucketing and the safety flag/count;
- a scores table with neither `rwr_score` nor `bfwalk_score` rejected by
  the schema;
- a gene with an empty `safetyLiabilities` list and a gene absent from the
  `target` parquet both coming out `unknown`, never `safe`, in
  `build_annotations` and end-to-end;
- `run()` restoring Stage 5's ranking (by `bfwalk_score`, else
  `rwr_score`) from an unsorted candidates file, with or without Stage 4;
- the annotations never changing that order.

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
  Targets versions, and the kind of graph (Stage 3's unweighted functional
  interactions on the GSEA pathways, or with Stage 4's user weights);
- a **sortable** ranked table, in Stage 5's order — `symbol`, `gene`
  (Ensembl), `rwr`, `bfwalk`, `tractability`, `safety` (click a
  header to sort; `rwr`/`bfwalk` columns appear only if those methods
  ran). Stage 6's `genetic_evidence_score` is not shown: it is 0 for every
  candidate;
- an **evidence-trace side panel**: clicking the symbol of one of the top 100
  candidates (`--n-traces`; the table still lists every candidate) opens a
  panel on the right of the table with that gene's trace, and highlights its
  row; clicking another symbol switches the panel to that gene, and `×` or
  Esc closes it (on a narrow screen the panel covers the page). The trace
  shows its ranking
  score, the disease-relevant (GSEA) pathways it shares with seed
  genes (with the number of seeds in each; and the disease-relevant pathways
  without seeds it is a member of), its Open
  Targets datatypes + scores (strongest first), and the signed interaction edges touching it (edges with a seed
  are badged "with seed"; interactions to non-graph genes are omitted,
  matching Stage 3's scoping).

Genes are displayed by HGNC symbol; the Ensembl gene ID stays the internal
key (its own column, and the id of each trace in the panel). Symbols are
resolved from `genes.parquet` at render time only.

Only `manifest.json` and `candidates_annotated.tsv` are required; every
other artifact is optional. A partial run directory still produces a
report — missing trace inputs (`gene_sets.parquet`, `interactions.parquet`,
`ot_disease_subset.parquet`, …) show that section as unavailable,
and no `genes.parquet` falls the symbol column back to Ensembl ids — a
header line lists whatever was missing.

**Expected output** (`example_data/example_run/report.html`, top of the
ranked table — same ordering as Stage 6's `candidates_annotated.tsv`):

```
#  symbol  gene             bfwalk  tractability  safety
1  GENE22  ENSG00000000022  0.5000  clinical      has_liabilities
2  GENE24  ENSG00000000024  0.5000  discovery     unknown
3  GENE23  ENSG00000000023  0.2500  unknown       unknown
```

Clicking `GENE22` opens the panel with: pathways shared with seed genes *Child Broad
Signaling* (10 seed genes); its `known_drug` Open Targets row (0.500,
chembl); and the interactions `GENE03 → GENE22`, `GENE05 → GENE22` and
`GENE07 — GENE22`, each badged "with seed" (`N29 → N22` is omitted: N29 is
not in the graph).

`python3 -m pytest tests/test_stage7_report.py tests/test_run_pipeline.py`
runs the tests (10): `evidence_trace`'s disease-pathway / graph-scoping /
datatype logic on a hand-built run, `build_report` producing a
self-contained symbol-labelled HTML with one trace per candidate, all in the
side panel (closed by default, after the table), and no genetic-evidence column,
traces (and symbol links / row highlighting) only for the top `--n-traces`, a partial run
directory (only `manifest.json` + `candidates_annotated.tsv`) still
rendering with the optional columns/trace sections marked unavailable, a
`run_pipeline --from-manifest --report` smoke test that Stages 2–6 chain
and produce every expected artifact, and a stage-failure producing a
helpful "which stage / which command / exit code" message rather than a
traceback.
