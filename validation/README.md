# Validation

Does the pipeline recover a disease's known genes, and does it do so better than simply sorting the Open
Targets Platform data by score? We build one graph per disease from a genome-wide disease-pathway GSEA
catalogue, hold out each of its genetically associated genes in turn and check how highly the pipeline ranks
it. We then rank the same held-out genes by their other Open Targets evidence and compare the two.

## Data

- **Disease-pathway GSEA catalogue** (Open Targets 26.03), in `validation/Disease_graph_catalogue/`. Not in
  git; available on request from Open Targets (EMBL-EBI). GSEA of every Open Targets disease's genes, ranked by
  their `genetic_association` score, against GO BP, Reactome and MSigDB Hallmark gene sets. Used here:
  - `genetic_disease_pathway_catalogue.parquet`: one row per disease and gene set (`es`, `pval`, `fdr`, ...);
  - `merged_gene_sets.gmt`: the gene sets (gene symbols);
  - `genetic_gene_scores_input.parquet`: each disease's genes and genetic-association scores (`diseaseId`,
    `gene`, `score`).
- **Open Targets Platform 26.06** `target`, `disease` and `association_by_datasource_direct`, and the **Reactome
  functional interactions** `FIsInGene_04142025_with_annotations.txt` (unzipped): Stage 1 of the pipeline
  downloads them into `data/`.

Genes are mapped from symbols to Ensembl gene IDs with the 26.06 `target` table, as Stage 1 does.

## Disease graph catalogue

```
python validation/build_graph_catalogue.py \
    --catalogue-dir validation/Disease_graph_catalogue \
    --target-dir data/opentargets/26.06/target \
    --fi-file data/reactome/FIsInGene_04142025_with_annotations.txt
```

One graph per disease, built exactly as Stage 3 builds it:
- **Pathways:** the disease's Reactome gene sets with `es > 0` and FDR < 0.05 (`--fdr`). The FDR is recomputed
  over Reactome only, because the catalogue corrected over all 8,217 gene sets. Diseases without such a pathway
  are left out (2,555 of the 3,766 diseases with GSEA results are kept).
- **Graph:** nodes are all genes on those pathways; edges are the curated Reactome functional interactions
  between them, directed and signed, all of weight 1.
- **Seeds:** the disease's genes with a genetic-association score > 0. Only seeds that are graph nodes can be
  used for propagation.

`--diseases EFO_0003767,MONDO_0004975` or `--n-diseases 10` builds only some diseases. Disease IDs are those
of Open Targets 26.03.

```
validation/Disease_graph_catalogue/disease_graph_catalogue/      (--out-dir changes the parent)
├── graph_catalogue_diseases.parquet      one row per disease: n_pathways, n_nodes, n_edges, n_isolated_nodes,
│                                         n_seeds, n_seeds_in_graph, n_seed_symbols_unmapped
└── EFO_0003767_disease_graph/            the same files as a pipeline run directory
    ├── graph_weight.npz, graph_sign.npz  the graph (sparse matrices)
    ├── graph_gene_index.json             the Ensembl gene of each matrix row / column
    ├── graph_metadata.json               the summary row, input paths and FDR cutoff
    ├── disease_seed_genes.tsv            gene_id, score
    └── disease_pathways.tsv              set_id, gene_set, source_db, pval, fdr
```

So Stage 5 can score any disease directly: `python pipeline/stage5_score_candidates.py --graph-dir <disease directory>`.

## Leave-one-out

```
python validation/leave_one_out.py \
    --catalogue validation/Disease_graph_catalogue/disease_graph_catalogue \
    --diseases EFO_0003767,MONDO_0004975 --jobs 8
```

For each disease, each seed in its graph is held out in turn: the other graph genes that are not seeds are the
candidates, the held-out gene among them, and Stage 5's BFWalk and RWR score them from the remaining seeds
(all seeds weighted equally). The held-out gene's rank among the candidates (ties count half) gives:
- **AUROC:** the fraction of the other candidates it outscores, averaged over the held-out genes;
- **recall@k:** the fraction of held-out genes ranked within the top k (`--k`).

A random ranking would score AUROC 0.500 and recall@k = k / number of candidates; on the ROC plots it is the
diagonal. A disease with fewer than 2 seeds in its graph is skipped. The graph is not rebuilt without the
held-out gene, although its pathways were selected by a GSEA that included it, so the results are somewhat
optimistic.

Without `--diseases` or `--n-diseases N`, every disease of the catalogue is run. `--max-seeds N` holds out a
random sample of N seeds per disease (`--seed`). A disease already scored from the same input files (sha256)
and parameters is not scored again; `--recompute` forces it.

```
validation/leave_one_out/                  (--out-dir)
├── leave_one_out_summary.tsv              one row per disease and method: AUROC, recall@k, random expectation
├── leave_one_out.log
├── roc_curve_all_diseases.png             all held-out genes of all diseases pooled
├── roc_curves_per_disease.png             one curve per disease and method (BFWalk blue, RWR orange shades)
├── roc_curves_per_disease_bfwalk.png      one curve per disease, BFWalk
├── roc_curves_per_disease_rwr.png         one curve per disease, RWR
└── EFO_0003767_leave_one_out/
    ├── scores_left_out_details.tsv        one row per held-out gene: scores, ranks, candidates, genetic evidence
    ├── metrics.tsv, leave_one_out.log
    ├── roc_curve.png
    ├── cache_key.json                     what the scores were computed from
    ├── full/candidate_scores.tsv          all seeds, nothing held out
    └── ENSG00000.../                      one directory per held-out gene
        ├── candidate_scores.tsv           every candidate's bfwalk_score and rwr_score
        ├── result.json                    the held-out gene's scores and ranks
        └── seed.log
```

The per-disease plots also draw the pooled curve of all diseases.

## This pipeline vs Open Targets Platform scores

```
python validation/compare_scores.py \
    --association-dir data/opentargets/26.06/association_by_datasource_direct \
    --disease-file data/opentargets/26.06/disease/disease.parquet
```

Run after `leave_one_out.py`, on its output (`--leave-one-out-dir`; `--diseases` / `--n-diseases` as above).
Each held-out gene is also ranked among the candidates of its own leave-one-out run by its highest Open
Targets score for the disease over every datasource except the genetic-association ones, which made it a
seed. This is "sorting the Open Targets data by score" without the answer. Per disease and method, a paired
Wilcoxon test on the held-out genes' AUROC answers whether the pipeline is `better`, `worse` or shows
`no_difference` (p < 0.05). The log also gives the Spearman correlations between the genetic evidence, Open
Targets, BFWalk and RWR scores.

The catalogue's 26.03 disease IDs are looked up in the 26.06 `disease` table, by ID or else by cross-reference
(26.06 keeps many older EFO IDs only there, e.g. `EFO_0003767` is `MONDO_0005265`). A disease not found is left
out and logged. The seeds come from 26.03 indirect genetic scores, the baseline from 26.06 direct associations.

```
validation/leave_one_out/
├── compare_scores_summary.tsv                    one row per disease and method: AUROCs, wins, p-value, verdict
├── compare_scores.log                            how many diseases get each answer
├── roc_curve_vs_open_targets_all_diseases.png
└── EFO_0003767_leave_one_out/
    ├── compare_scores_details.tsv                the held-out genes with their Open Targets score and rank
    ├── compare_scores.log
    └── roc_curve_vs_open_targets.png
```

Tests: `tests/test_graph_catalogue.py` and `tests/test_leave_one_out.py`, on small hand-built inputs.
