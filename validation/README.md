# Validation

Does the pipeline recover disease's known genes, and does it do so better than simply sorting the Open Targets data by score?

## Disease graph catalogue

```
python validation/build_graph_catalogue.py \
    --catalogue-dir validation/Disease_graph_catalogue \
    --target-dir data/opentargets/26.06/target \
    --fi-file data/reactome/FIsInGene_04142025_with_annotations.txt
```

All three paths are required: the GSEA catalogue, the Open Targets `target` table and the unzipped Reactome
functional-interaction file (Stage 1 downloads the last two into `data/`). Add `--n-diseases 10` to build
only the first 10 included diseases for a test run, or `--diseases EFO_0003767,MONDO_0004975` for diseases
of interest.

One graph per Open Targets disease, built from the genome-wide disease-pathway GSEA catalogue in
`validation/Disease_graph_catalogue/` (Open Targets 26.03; see its `README.md` for how the GSEA
was run). Only Reactome pathways are used for now.

How it is built:
1. **Pathways.** The Reactome rows of `genetic_disease_pathway_catalogue.parquet`. The catalogue's FDR was
   corrected over all 8,217 gene sets (GO BP, Reactome, Hallmark), so the Benjamini-Hochberg FDR is
   recomputed per disease over its Reactome pathways only (both ES signs; rows with a null p-value are
   not tested). A disease's pathways are those with `es > 0` and recomputed FDR < 0.05 (`--fdr`).
2. **Diseases.** Of the 22,771 diseases with a genetic score, 19,005 were never tested by GSEA (fewer than
   25 scored genes in the gene-set library) and 3,766 were. Of those, 2,555 have at least one Reactome
   pathway with `es > 0` and FDR < 0.05 and are included; 1,211 have none and are left out.
3. **Graph.** Exactly Stage 3's graph (`pipeline/stage3_build_graph.py`): nodes are all genes on the
   disease's pathways (gene lists from `merged_gene_sets.gmt`), edges are the curated Reactome functional
   interactions between them (`data/reactome/FIsInGene_04142025_with_annotations.txt`), directed and
   signed (+1, -1, 0), unweighted. An interaction without a direction is an edge in both directions.
4. **Seeds.** Each gene's genetic-association score for the disease from `genetic_gene_scores_input.parquet`
   (Open Targets `genetic_association`, indirect, in (0, 1]); genes without a score get 0. These scores are
   the seeds for propagation.
5. **Gene IDs.** Symbols (GMT, FI file, scores) are mapped to Ensembl gene IDs with the Open Targets 26.06
   `target` table (`data/opentargets/26.06/target/`), exactly as Stage 1 does: approved symbol first, then a
   unique synonym or obsolete symbol. 96.5% of the GMT's Reactome symbols map (most of the rest are not
   genes, e.g. rRNAs and viral proteins) and 99.4% of the FI file's. For 81% of the pathways the mapped gene
   count equals the catalogue's `geneset_size` (the catalogue mapped onto 26.03 symbols); 95% are within 11%.

Written to `validation/Disease_graph_catalogue/disease_graph_catalogue/` (`--out-dir` changes the parent),
one directory per included disease:

```
disease_graph_catalogue/
├── graph_catalogue_diseases.parquet      one summary row per included disease
├── EFO_0003767_disease_graph/
│   ├── graph_weight.npz                  Stage 3's graph: weight 1 per directed edge
│   ├── graph_sign.npz                    edge signs (+1, -1, 0)
│   ├── graph_gene_index.json             the Ensembl gene of each matrix row / column
│   ├── graph_metadata.json               the summary row plus the input paths and FDR cutoff
│   ├── disease_seed_genes.tsv            gene_id, score: the seeds, strongest first
│   └── disease_pathways.tsv              set_id (R-HSA ID), gene_set, source_db, pval, fdr (recomputed)
└── ...
```

Each disease directory has the same files as a pipeline run directory, so Stage 5 scores a disease directly:

```
python pipeline/stage5_score_candidates.py --graph-dir validation/Disease_graph_catalogue/disease_graph_catalogue/EFO_0003767_disease_graph
```

`graph_catalogue_diseases.parquet` columns: `diseaseId`, `n_pathways`, `n_nodes`, `n_edges`,
`n_isolated_nodes`, `n_seeds`, `n_seeds_in_graph`, `n_seed_symbols_unmapped`. `n_seeds` counts all the
disease's scored genes and `n_seeds_in_graph` those on its pathways: seeds off the pathways cannot take part
in propagation (Stage 5 drops them with a log line). Disease IDs are Open Targets' own (e.g. inflammatory bowel disease
is `EFO_0003767`, not `MONDO_0005265`).

The log ends with a summary (median graph size, and how many diseases have no seed in their graph or no
edge). For example, Alzheimer disease (`MONDO_0004975`) has 63 pathways, 5,173 genes, 135,967 edges and
507 of its 1,700 seeds in the graph.

Other options: `--fdr` (default 0.05) and `--out-dir` (default: `--catalogue-dir`). Tested by
`tests/test_graph_catalogue.py` on tiny synthetic inputs in the same formats.

## Leave-one-out

```
python validation/leave_one_out.py \
    --catalogue validation/Disease_graph_catalogue/disease_graph_catalogue \
    --diseases EFO_0003767,MONDO_0004975 --jobs 8
```

Runs on the disease graph catalogue above. Choose diseases with `--diseases` (comma-separated IDs) or
`--n-diseases N` (the first N, sorted by ID); without either, every disease of the catalogue is run.

For each disease, every seed that is in its graph is held out in turn:
1. the seed is removed from the seeds;
2. every graph gene that is not one of the remaining seeds is a candidate, the held-out gene among them;
3. the candidates are scored with Stage 5's BFWalk and RWR (`pipeline/stage5_score_candidates.py`, uniform
   seeds), and the held-out gene's rank among them is recorded (mid-rank: ties count half, 1 = best).

Over the held-out seeds of a disease it reports, for BFWalk, RWR and a random ranking: **AUROC** (per gene, the
fraction of the other candidates it outscores, averaged) and **recall@k** (the fraction ranked within the top k),
with a permutation p-value against `--n-random` random rankings. The ROC curve plots, for each fraction of
candidates taken from the top, the fraction of held-out seeds recovered; its area is the AUROC.

Seeds that are not in the graph (on none of the disease's pathways) cannot be held out and are left out.
A disease with fewer than 2 seeds in its graph is skipped. The graph is not rebuilt for each held-out seed:
its pathways were selected by a GSEA that included that seed's score, so the result measures recovery given
the disease's pathways, and is somewhat optimistic.

Written to `--out-dir` (default `validation/leave_one_out/`):

```
leave_one_out/
├── leave_one_out_summary.tsv          one row per disease and method: AUROC, recall@k, p-values
├── leave_one_out.log
├── roc_curve_all_diseases.png         all held-out seeds of all diseases pooled
└── EFO_0003767_leave_one_out/
    ├── scores_left_out_details.tsv    one row per held-out seed: scores, ranks, candidates, genetic evidence
    ├── metrics.tsv, leave_one_out.log
    ├── roc_curve.png                  BFWalk, RWR and random
    ├── full/candidate_scores.tsv      all seeds, nothing held out
    └── ENSG00000.../                  one directory per held-out seed
        ├── candidate_scores.tsv       every candidate's bfwalk_score and rwr_score
        ├── result.json                the held-out gene's scores and ranks
        └── seed.log
```

Other options: `--max-seeds N` (a random sample of N seeds per disease), `--jobs` (seeds scored in
parallel), `--k`, `--n-random`, `--seed`, `--restart-prob`, `--bfwalk-alpha`.

## Pipeline vs Open Targets

```
python validation/compare_scores.py \
    --association-dir data/opentargets/26.06/association_by_datasource_direct
```

Run after `leave_one_out.py` (same `--diseases` / `--n-diseases` options; reads `--leave-one-out-dir`,
default `validation/leave_one_out/`). For each disease it ranks every held-out seed among the candidates of its
own leave-one-out run by its Open Targets score for the disease over every datatype except
`genetic_association` (the evidence that made it a seed): "sorting the Open Targets data by score" without
the answer. It compares that ranking with BFWalk and RWR (AUROC, recall@k, a paired Wilcoxon test on the
per-seed AUROC), correlates genetic evidence, Open Targets, BFWalk and RWR scores, and answers per method
"is the pipeline better than sorting the Open Targets data by score?" (`better`, `worse` or `no_difference`
at p < 0.05).

The baseline uses the local Open Targets 26.06 **direct** associations, while the seeds come from the
catalogue's 26.03 **indirect** genetic scores; a disease missing from the 26.06 data gets a baseline of 0 for
every gene (logged).

Writes, in each `<diseaseId>_leave_one_out/`: `compare_scores.log`, `compare_scores_details.tsv` and
`roc_curve_vs_open_targets.png`; and in `--leave-one-out-dir`: `compare_scores_summary.tsv` (one row per disease
and method, with the verdict), `compare_scores.log` (how many diseases each answer) and
`roc_curve_vs_open_targets_all_diseases.png`.

Both scripts are tested by `tests/test_leave_one_out.py` on a hand-built catalogue.
