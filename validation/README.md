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

For each disease, each seed in its graph is held out in turn: the graph genes that are not remaining seeds are
the candidates (the held-out gene among them), and Stage 5's BFWalk and RWR score them from the remaining seeds,
all weighted equally. Everything below is computed from each held-out gene's rank among its candidates (tied
genes share the average place). A disease with fewer than 2 seeds in its graph is skipped. The graph is not
rebuilt without the held-out gene, although its pathways were selected by a GSEA that included it, so the
results are somewhat optimistic.

Without `--diseases` or `--n-diseases N`, every disease of the catalogue is run. `--max-seeds N` holds out a
random sample of N seeds per disease (`--seed`). A disease already scored from the same input files (sha256)
and parameters is not scored again; `--recompute` forces it.

### Combining BFWalk and RWR

BFWalk and RWR propagate from the same seeds in different ways, so each may recover genes the other misses.
Five rankings merge a candidate's two ranks into one; they use ranks, not scores, because the two scores are on
different scales (BFWalk's are about 100 times larger).
- **RRF**, reciprocal rank fusion (Cormack et al. 2009): the sum over the two methods of 1 / (60 + rank). It is a
  standard way to fuse rankings from different systems: a gene near the top of both lists wins, but one very
  high rank still counts.
- **Borda** count: in each method, one point per candidate ranked below the gene, summed; this orders genes
  exactly like the sum of their ranks. The simplest consensus: it asks whether the two methods' average opinion
  beats either one.
- **RRA**, robust rank aggregation (Kolde et al. 2012): compares the gene's two normalised ranks (rank / number of
  candidates) with uniformly random ones; the score is the smallest probability that random ranks are as good
  as its better rank, or as both. Built for noisy gene lists, it puts first the genes that either or both methods
  rank unusually high.
- **union top-k**: the genes in either method's top k; a gene joins at k = its better rank. It asks whether
  following up both methods' top hits finds more known genes per gene followed up.
- **intersection top-k**: the genes in both methods' top k; a gene joins at k = its worse rank. It asks whether
  agreement between the two methods is a stronger signal than either method alone.

### Values

Each value is given for BFWalk, RWR, each combination, and a uniformly random ranking of the same candidates
(row `random`).
- **AUROC**: the fraction of the other candidates the held-out gene outranks, averaged over the held-out genes
  (random: 0.5). It does not depend on the number of candidates, so it compares diseases directly, but it weighs
  the whole list, not only the top that would be followed up.
- **AUPRC**: the area under the precision-recall curve of the held-out genes pooled, every run cut at the same
  depth k (recall: the fraction of held-out genes in the top k; precision: those genes / all candidates in the
  top-k lists), computed as average precision, without interpolation. With one positive among thousands of
  candidates it is dominated by the top of the list, which is where a combination has to help.
- **AUPRC fold enrichment**: AUPRC / prevalence, the prevalence being held-out genes / candidates (about 1 /
  number of candidates). A random ranking's AUPRC is about the prevalence, which differs between diseases, so
  AUPRC is compared as is only within a disease and across diseases as this ratio. With few held-out genes even
  a random ranking scores above 1 (about 1.2 with 50 genes, 2.5 with 3): read it against the `random` row.
- **precision@k and its fold enrichment** at the follow-up depth (`--follow-up-k`, default 20): the held-out
  genes in the top-k lists / all candidates in them, as is and / prevalence. It is the hit rate among the genes
  that would actually be followed up; with one held-out gene per run it is at most 1 / k, so the fold enrichment
  is easier to read.
- **recall@k** (`--k`): the fraction of held-out genes ranked within the top k. It says how many of the known
  genes a list of k would find.

### Is a combination better?

`combined_score_tests.tsv` has one row per combination and per-gene measure: the measure's mean under the
combination, BFWalk, RWR and a random ranking (one random sample, `--seed`), and one-sided Wilcoxon
signed-rank p-values that the combination's is greater. The two measures are the gene's AUROC (the whole list)
and its average precision vs random (the top of the list): 1 / rank divided by a random ranking's expected
1 / rank, so 1 = random. The test is paired, each held-out gene compared with itself under the two rankings,
which removes the large differences in how hard genes are to recover. Rows are given for all diseases pooled
(`all`) and per disease; the p-values are not corrected for multiple testing.

### Plots

`roc_pr_curve.png` (one disease) and `roc_pr_curve_combined_all_diseases.png` (all diseases pooled) show the ROC
curves with their AUROC, and the precision-recall curves with their AUPRC and precision@k; the pooled plot
gives precision, AUPRC and precision@k as fold enrichment ("x"), since diseases differ in prevalence. Precision
is on a log scale. The dot on each curve is precision at the follow-up depth; it lies below the steps, which
show precision where each held-out gene is recovered. `roc_curve_all_diseases.png` shows only the BFWalk and RWR
ROC curves.

```
validation/leave_one_out/                    (--out-dir)
├── leave_one_out_summary.tsv                one row per disease and ranking (bfwalk, rwr, rrf, borda, rra,
│                                            union, intersection, random): the values above
├── combined_score_tests.tsv                 each combination vs BFWalk, RWR, random: means, p-values
├── roc_pr_curve_combined_all_diseases.png   all rankings, all diseases pooled
├── roc_curve_all_diseases.png               BFWalk and RWR, all diseases pooled
├── leave_one_out.log
└── EFO_0003767_leave_one_out/
    ├── metrics.tsv                          one row per ranking: the values above
    ├── roc_pr_curve.png                     all rankings
    ├── scores_left_out_details.tsv          one row per held-out gene: scores, ranks under every ranking,
    │                                        candidates, genetic evidence
    ├── leave_one_out.log
    ├── cache_key.json                       what the scores were computed from
    ├── full/candidate_scores.tsv            all seeds, nothing held out
    └── ENSG00000.../                        one directory per held-out gene
        ├── candidate_scores.tsv             every candidate's bfwalk_score and rwr_score
        ├── result.json                      the held-out gene's scores and ranks
        └── seed.log
```

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
Wilcoxon test on the held-out genes' AUROC (two-sided) answers whether the pipeline is `better`, `worse` or shows
`no_difference` (p < 0.05). The log also gives the Spearman correlations between the genetic evidence, Open
Targets, BFWalk and RWR scores.

The catalogue's 26.03 disease IDs are looked up in the 26.06 `disease` table, by ID or else by cross-reference
(26.06 keeps many older EFO IDs only there, e.g. `EFO_0003767` is `MONDO_0005265`). A disease not found is left
out and logged. The seeds come from 26.03 indirect genetic scores, the baseline from 26.06 direct associations.

```
validation/leave_one_out/
├── compare_scores_summary.tsv                    one row per disease and method: AUROC and AUPRC fold
│                                                 enrichment of the pipeline and Open Targets, wins, p-value, verdict
├── compare_scores.log                            how many diseases get each answer
├── roc_curve_vs_open_targets_all_diseases.png
└── EFO_0003767_leave_one_out/
    ├── compare_scores_details.tsv                the held-out genes with their Open Targets score and rank
    ├── compare_scores.log
    └── roc_curve_vs_open_targets.png
```

Tests: `tests/test_graph_catalogue.py` and `tests/test_leave_one_out.py`, on small hand-built inputs.
