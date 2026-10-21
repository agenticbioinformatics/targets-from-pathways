# Targets from Pathways - User guide

Per-stage reference: what each stage does, what it reads and writes, and its options. For
installation and a quick overview see the [main README](../README.md).

Throughout, `data/` is the source cache (`--data-dir`) and `out/` is the run directory. Every stage
reads from and writes into the run directory, next to Stage 1's `manifest.json`. Run the scripts from
the repository root (`python pipeline/<script> ...`); every script has `--help`.

**Conventions**
- **Genes are Ensembl gene IDs** (`ENSG` + 11 digits) in every file. Stage 1 is the only stage that maps
  identifiers; the report shows symbols for reading only.
- **Every table is validated** against its schema in `pipeline/stage0_schemas.py` when it is read and
  written, so a malformed file stops the stage with a message naming the column.
- **Seed genes** are the disease's genes with an Open Targets `genetic_association` score > 0. They
  drive GSEA (Stage 2) and the propagation (Stage 5), and are never candidates themselves.


## Whole pipeline (`pipeline/run_pipeline.py`)

```
python pipeline/run_pipeline.py --disease MONDO_0005265 --data-dir data/ --out-dir out/ --report
```

Runs Stages 1-6 one after the other, each as its own subprocess with its normal command line. If a
stage fails, `run_pipeline.py` says which stage, with its exit code and exact command, and stops; the
stage's own error message is just above. `--report` also writes `out/report.html` (Stage 7).

| Option | Meaning |
|---|---|
| `--disease`, `--data-dir`, `--out-dir` | Passed to Stage 1. |
| `--from-manifest out/manifest.json` | Skip Stage 1 and rerun Stages 2-6 on an existing run directory. |
| `--method` | Stage 5's method: `bfwalk` (default), `rwr` or `bfwalk,rwr`. |
| `--node-weights`, `--edge-weights` | Stage 4's optional weight files (see Stage 4). |
| `--report` | Also write the HTML report. |
| anything else | Passed on to Stage 1, e.g. `--allow-low-coverage` or `--min-set-size 5`. |

**Installation check.** Every run, `--help` included, first imports every stage. If a package is
missing it stops and lists what is missing, so `python pipeline/run_pipeline.py --help` tells you
whether the installation is complete. A missing BFWalk only gives a warning, because
`--method rwr` works without it.

To change Stage 4's other options or Stage 5's parameters, run those stages on their own afterwards,
then Stage 6 and Stage 7.


## Stage 1 - ingest (`pipeline/stage1_ingest.py`)

```
python pipeline/stage1_ingest.py --disease MONDO_0005265 --data-dir data/ --out-dir out/
```

Downloads the pinned source files into `--data-dir` on the first run and checks each against a
stored sha256 (later runs reuse the cache):
- Open Targets Platform 26.06: `association_by_datasource_direct`, `target`, `disease`;
- Reactome 97: `ReactomePathways.gmt`, `Ensembl2Reactome.txt`, and the functional-interaction file
  `FIsInGene_04142025_with_annotations.txt`.

It then maps every gene identifier to an Ensembl gene ID, using the Open Targets `target` table as the
authority, and keeps the Open Targets evidence for the requested disease.

**Writes:** `genes.parquet` (Ensembl ID, symbol, synonyms), `gene_sets.parquet` (Reactome pathways),
`interactions.parquet` (Reactome functional interactions, with direction and sign),
`ot_disease_subset.parquet` (the disease's Open Targets evidence, one row per gene and datasource),
`coverage_report.json` (how many identifiers were mapped), `scale_report.json` (pathway sizes and
interaction counts) and `manifest.json` (versions, checksums and parameters). Later stages find each
file next to `manifest.json` and refuse one whose checksum no longer matches the manifest.

| Option | Default | Meaning |
|---|---|---|
| `--disease` | required | EFO, MONDO or another ontology ID that Open Targets knows. |
| `--no-download` | off | Never use the network; `--data-dir` must already hold the files. |
| `--min-set-size`, `--max-set-size` | 10, 200 | Keep only pathways with this many genes. |
| `--fi-curated-only` | on | Keep only curated functional interactions (predicted ones carry no reliable direction or sign). |
| `--fi-min-score` | none | Drop functional interactions scoring below this. |
| `--allow-low-coverage` | off | Do not stop when fewer than 90% of identifiers map, or the disease has fewer than 5 seed genes. |
| `--seed` | 0 | Random seed, recorded in the manifest and used by GSEA. |


## Stage 2 - disease pathways (`pipeline/stage2_gsea_discovery.py`)

```
python pipeline/stage2_gsea_discovery.py --manifest out/manifest.json
```

Ranks the seed genes by their strongest `genetic_association` score and runs preranked GSEA
([blitzgsea](https://github.com/MaayanLab/blitzgsea)) against the Reactome pathways. Only genetic
evidence is used: Open Targets' overall score includes Reactome-derived evidence, which would make
testing against Reactome pathways circular. Near-duplicate pathways (Jaccard > 0.7, common in
Reactome's hierarchy) are collapsed to the smallest one first, and p-values are BH-FDR corrected.

**Writes:** `disease_pathways.tsv` (the significant pathways: `set_id`, `gene_set`, `source_db`,
`pval`, `fdr`), `disease_seed_genes.tsv` (`gene_id`, `score`) and `benchmark_holdout_genes.tsv`.

| Option | Default | Meaning |
|---|---|---|
| `--pval-threshold` | 0.05 | Nominal p-value cutoff. |
| `--fdr-threshold` | 0.1 | FDR cutoff. |
| `--benchmark-holdout-file` | none | Genes to treat as not associated with the disease (see below). |

**Benchmark holdout.** `--benchmark-holdout-file` lists genes, one per line (symbol, synonym or Ensembl
ID; `#` starts a comment). They are removed from the seeds and the GSEA ranking, and Stage 4 ignores
their own evidence, so a held-out gene can only come back as a ranked candidate. This is what
`validation/leave_one_out.py` uses. A held-out gene can be ranked only if the remaining seeds still
detect one of its pathways; otherwise it is not in the graph.

**Diseases with few seed genes.** GSEA needs at least 5 seed genes in a pathway, and blitzgsea's
calibration can fail on small rankings. In those cases Stage 2 logs a warning and writes an empty
`disease_pathways.tsv`, and Stage 3 then stops with a message: without a disease pathway there is no
graph.


## Stage 3 - graph (`pipeline/stage3_build_graph.py`)

```
python pipeline/stage3_build_graph.py --manifest out/manifest.json
```

Builds a directed, signed, unweighted graph:
- **nodes:** every gene on the disease pathways in `disease_pathways.tsv`;
- **edges:** the Reactome functional interactions between two of those genes. Sharing a pathway adds
  no edge, and interactions with genes outside these pathways are left out. An interaction without a
  direction becomes an edge in both directions.

Every edge has weight 1 (interaction confidence is not used as a weight). The sign (+1 activating,
-1 inhibiting, 0 unknown) is stored separately. The seed genes are read only to check that some are
in the graph: Stage 3 stops with a message if `disease_pathways.tsv` is empty or if no seed gene is on
a disease pathway.

**Writes:** `graph_weight.npz` and `graph_sign.npz` (sparse matrices, row = source gene, column =
target gene), `graph_gene_index.json` (the gene of each row and column) and `graph_metadata.json`
(counts and source versions).

| Option | Default | Meaning |
|---|---|---|
| `--disease-pathways` | next to the manifest | Another pathway table. |
| `--seed-genes` | next to the manifest | Another seed list. |
| `--out-dir` | the manifest's directory | Where to write the graph. |


## Stage 4 - optional weights (`pipeline/stage4_weights.py`)

```
python pipeline/stage4_weights.py --manifest out/manifest.json \
    [--node-weights node_weights.tsv] [--edge-weights edge_weights.tsv]
```

Stage 4 does two independent things.

**1. Your own weights (optional).** If you have prior knowledge, such as expression in a relevant
tissue, a druggability score or interaction confidences from another database, you can turn it into
weights that Stage 5's propagation follows: a heavier edge carries more of the seeds' signal. Both
files are optional and can be used alone or together; without either, the graph stays unweighted.

- `--node-weights`: a TSV with the columns `gene_id` and `weight`. Each edge gets the **mean** of its
  two genes' weights (`--node-weight-mode product` multiplies them instead, which penalises an edge
  more when either gene has a low weight).
- `--edge-weights`: a TSV with the columns `gene_a`, `gene_b` and `weight`. A row weights the edge
  between the two genes in whichever directions the graph has it, so the order of the genes does not
  matter; listing a pair twice is an error.
- **Both files:** an edge's weight is its edge weight times its combined node weight. For example,
  with node weights A = 0.5 and B = 1 and an edge weight of 0.4 for (A, B), the edge A → B weighs
  0.4 × (0.5 + 1) / 2 = 0.3.

```
gene_id          weight              gene_a           gene_b           weight
ENSG00000105974  2.0                 ENSG00000105974  ENSG00000141510  0.8
ENSG00000141510  0.5
```

Rules:
- Genes are Ensembl gene IDs; weights are finite numbers > 0 (a 0 would silently cut the edge out).
- Genes and edges missing from a file weigh `--default-node-weight` / `--default-edge-weight` (both 1).
  Weights are relative, so set the defaults to match your scale: if all your weights are below 1, the
  genes and edges you did not list (weight 1) would outweigh the ones you did.
- Stage 4 never adds or removes an edge. Rows for genes or gene pairs that are not in the graph are
  skipped, and the log says how many.
- The weights are always applied to Stage 3's graph with every edge reset to 1 first, so rerunning
  Stage 4 with other files, or with none, replaces the previous weighting instead of stacking on it.

`graph_metadata.json` records what was used under `stage4_user_weights`, and the report header says
whether the graph is weighted.

**2. Genetic-evidence scores, for Stage 6.** Each graph gene gets its strongest Open Targets
`genetic_association` score. Other evidence types are not counted: `affected_pathway` contains
Reactome evidence (circular, as in Stage 2), and `known_drug` would favour genes that are already drug
targets for the disease. These scores are **never used as weights**. With Stage 2's seeds every
candidate scores 0 (a gene scoring > 0 is a seed), so the score is informative only when Stage 5 is
given other seeds. Benchmark holdout genes score 0.

**Writes:** the four graph files again (`graph_weight.npz` replaced, `graph_metadata.json` extended
with `genetic_evidence_score` and `stage4_user_weights`), and `genetic_evidence_scores.tsv`.

| Option | Default | Meaning |
|---|---|---|
| `--node-weights`, `--edge-weights` | none | The weight files described above. |
| `--node-weight-mode` | `mean` | `mean` or `product` of an edge's two node weights. |
| `--default-node-weight`, `--default-edge-weight` | 1 | Weight of a gene or edge missing from its file. |
| `--graph-dir`, `--out-dir` | the manifest's directory | Where to read and write the graph. |
| `--holdout-genes` | next to the manifest | Stage 2's `benchmark_holdout_genes.tsv`. |


## Stage 5 - candidate scoring (`pipeline/stage5_score_candidates.py`)

```
python pipeline/stage5_score_candidates.py --graph-dir out/
```

Every graph gene that is not a seed is a candidate. Candidates are scored by how strongly the
seeds' signal reaches them over the directed graph:
- **`bfwalk`** (default): [BFWalk](https://github.com/jedrzejkubica/BFWalk), backtrack-free network
  propagation from the seeds. It needs edge weights in ]0, 1], so all weights are divided by the
  largest one.
- **`rwr`**: random walk with restart (personalised PageRank), restarting uniformly on the seeds.

BFWalk is less biased towards highly connected hub genes than RWR: use it to find candidates beyond
well-studied hubs, and add RWR when you prefer well-connected genes close to the seeds, e.g. to recover
established targets. With both methods (`--method bfwalk,rwr`) each gets its own column and the table
is ranked by `bfwalk_score`. Seeds that are not graph nodes are dropped with a log line.

**Writes:** `candidate_scores.tsv` (`gene`, `bfwalk_score` and/or `rwr_score`), ranked.

| Option | Default | Meaning |
|---|---|---|
| `--method` | `bfwalk` | `bfwalk`, `rwr` or `bfwalk,rwr`. |
| `--seed-genes` | `<graph-dir>/disease_seed_genes.tsv` | Your own seeds: any TSV with a `gene_id` column. |
| `--bfwalk-alpha` | 0.5 | BFWalk attenuation coefficient, between 0 and 1. |
| `--bfwalk-threads` | 0 (all cores) | BFWalk threads. |
| `--restart-prob` | 0.5 | RWR restart probability, between 0 and 1. |
| `--out` | `<graph-dir>/candidate_scores.tsv` | Output file. |


## Stage 6 - annotation (`pipeline/stage6_annotate_context.py`)

```
python pipeline/stage6_annotate_context.py --manifest out/manifest.json
```

Adds columns for reading; it filters nothing and keeps Stage 5's order:
- `genetic_evidence_score`: from Stage 4 (left out if Stage 4 did not run);
- `tractability`: `clinical` (an Open Targets clinical-stage flag), `discovery` (any other positive
  flag) or `unknown`;
- `safety`: `has_liabilities` or `unknown`, with `n_safety_liabilities`. There is no "safe" value:
  Open Targets' safety data is sparse, so no annotation is not evidence of safety.

**Writes:** `candidates_annotated.tsv`.

| Option | Default | Meaning |
|---|---|---|
| `--candidates` | next to the manifest | Stage 5's `candidate_scores.tsv`. |
| `--graph-dir` | the manifest's directory | Where Stage 4's `graph_metadata.json` is. |
| `--target-parquet` | from the manifest | The Open Targets `target` table. |
| `--out` | next to the candidates | Output file. |


## Stage 7 - report (`pipeline/stage7_report.py`)

```
python pipeline/stage7_report.py --run-dir out/
```

Writes one self-contained HTML file (no server or internet needed) with:
- a header: disease, numbers of seeds and candidates, data versions, and whether the graph is weighted;
- a sortable table of every candidate: symbol, Ensembl ID, scores, tractability, safety;
- an **evidence trace** for each of the top `--n-traces` candidates (default 100): click a symbol to open
  it in a panel on the right (× or Esc closes it). It shows the disease pathways the gene shares with
  seeds (and how many seeds each has), its Open Targets evidence, and its interactions in the graph,
  with those to a seed marked.

Only `manifest.json` and `candidates_annotated.tsv` are required. If other files are missing, the
report still renders and says which parts are unavailable.

| Option | Default | Meaning |
|---|---|---|
| `--run-dir` | required | The run directory. |
| `--out` | `<run-dir>/report.html` | Output file. |
| `--n-traces` | 100 | Number of top candidates with an evidence trace. |


## Tests

```
python -m pytest
```

The tests need no network and no downloads. Stage 1 runs on tiny Open Targets and Reactome files in
`tests/fixtures/`. Stages 2-7 run on a synthetic Stage 1 output written by
`tests/fixtures/synthetic_run.py`: 30 genes, 18 seeds and 7 hand-designed pathways (an ancestor/child
near-duplicate pair, a pathway with too few seeds for GSEA, a seedless one and a noise control).
`tests/test_end_to_end.py` checks every hand-off between stages against its schema.
