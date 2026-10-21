# Validation

Does the pipeline recover disease's known genes, and does it do so better than simply sorting the Open Targets data by score?

## Disease graph catalogue

```
python validation/build_graph_catalogue.py \
    --catalogue-dir validation/Disease_pathways_graph_catalogue \
    --target-dir data/opentargets/26.06/target \
    --fi-file data/reactome/FIsInGene_04142025_with_annotations.txt
```

All three paths are required: the GSEA catalogue, the Open Targets `target` table and the unzipped Reactome
functional-interaction file (Stage 1 downloads the last two into `data/`). Add `--n-diseases 10` to build
only the first 10 included diseases for a test run, or `--diseases EFO_0003767,MONDO_0004975` for diseases
of interest.

One graph per Open Targets disease, built from the genome-wide disease-pathway GSEA catalogue in
`validation/Disease_pathways_graph_catalogue/` (Open Targets 26.03; see its `README.md` for how the GSEA
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

Written to `validation/Disease_pathways_graph_catalogue/`:

| File | One row per | Columns |
|---|---|---|
| `graph_catalogue_diseases.parquet` | included disease | `diseaseId`, `n_pathways`, `reactome_ids` (list), `n_nodes`, `n_edges`, `n_isolated_nodes`, `n_seeds`, `n_seeds_in_graph`, `n_seed_symbols_unmapped` |
| `graph_catalogue_nodes.parquet` | (disease, graph gene) | `diseaseId`, `gene_id`, `seed_score` (0 for a non-seed) |
| `graph_catalogue_edges.parquet` | (disease, directed edge) | `diseaseId`, `gene_a`, `gene_b`, `sign` |

`n_seeds` counts all the disease's scored genes, `n_seeds_in_graph` those on its pathways; seeds off the
pathways cannot take part in propagation. Disease IDs are Open Targets' own (e.g. inflammatory bowel disease
is `EFO_0003767`, not `MONDO_0005265`).

The log ends with a summary (median graph size, and how many diseases have no seed in their graph or no
edge). For example, Alzheimer disease (`MONDO_0004975`) has 63 pathways, 5,173 genes, 135,967 edges and
507 of its 1,700 seeds in the graph.

Other options: `--fdr` (default 0.05) and `--out-dir` (default: `--catalogue-dir`). Tested by
`tests/test_graph_catalogue.py` on tiny synthetic inputs in the same formats.

## Leave-one-out

```
python validation/leave_one_out.py --disease MONDO_0005265 --from-manifest out/manifest.json --jobs 8
```

For each seed gene (one with a genetic-association score > 0 for the disease), the pipeline is rerun as
if that gene were not associated with the disease: it is removed from the seed set, its score and rank among the candidates are recorded.
A held-out gene that is on no pathway GSEA still detects without it drops out of the graph; it is recorded as not ranked, and the metrics are computed over the genes that were ranked.

- `validation/scores_left_out.tsv`: `gene`, `score` (BFWalk by default; `--method rwr` for RWR)
- `validation/leave_one_out.log`: AUROC and recall@k for BFWalk and RWR, random ranking

Each gene is a full Stage 2-5 run; add `--max-genes 50` to try a random sub-sample
first. Without an existing `out/`, use `--data-dir data/` instead of `--from-manifest` to run Stage 1.

## Pipeline vs Open Targets

```
python validation/compare_scores.py
```

Ranks the same held-out genes by their Open Targets score for the disease, leaving out the genetic
evidence that was hidden from the pipeline, and compares that ranking with BFWalk and RWR. It also
correlates genetic evidence, Open Targets, BFWalk and RWR scores. `validation/compare_scores.log` ends
with a YES/NO answer per method to "is the pipeline better than sorting the Open Targets data by score?".
