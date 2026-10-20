# Targets from Pathways

> The scientific idea was introduced by Jędrzej Kubica, Polina Rusina, Siddharth Sethi and Elvis Poku-Adusei at the Open Targets Hackathon on October 21-22, 2025. The technology was created by agentic AI provided by Claude Sonnet 4.6. The original project is available under the MIT license at [https://github.com/jedrzejkubica/targets-from-pathways](https://github.com/jedrzejkubica/targets-from-pathways)


## Quick start

**Install:** `pip install -r requirements.txt` (Python 3.9+, CPU only).

> **Note:** [BFWalk](https://github.com/jedrzejkubica/BFWalk) (optional Stage 5 scoring method) is still in development and is not on PyPI yet — `requirements.txt` installs it from TestPyPI. It will be available via `pip install bfwalk` shortly. To install it on its own: `pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ bfwalk`. Installing it compiles C code, so it needs gcc, zlib headers (e.g. `zlib-devel`) and OpenMP.

**Demo** — a tiny synthetic Stage 1 output (30 genes, 18 seeds) is provided:

```
python3 example_data/build_gsea_example.py
```

```
python3 pipeline/run_pipeline.py --from-manifest example_data/example_run/manifest.json --report
```

Expected: the log ends with `Scoring 7 candidates from 18 seeds` and
`Report written (7 candidates)`; `example_data/example_run/candidates_annotated.tsv`
ranks the 7 non-seed graph genes, led by `ENSG00000000022` (composite 0.71).
The 5 genes of the seedless pathway (`ENSG...029`–`033`) are not in the graph.

**Tests:** `python3 -m pytest` — ~120 offline tests to check if the pipeline is set up correctly.


## Overview

Given a disease, this pipeline proposes candidate drug targets based on
pathway topology, using public data: Open Targets (disease–target genetic
associations, tractability, safety) and Reactome (pathway gene sets and
functional interactions). The disease's genes with Open Targets
`genetic_association` evidence are the **seed genes**. The pipeline discovers
disease-relevant pathways (GSEA of the seed genes), builds a pathway-based
gene–gene graph from all genes on those pathways plus all genes on any
pathway that contains a seed gene, applies genetic evidence to the graph as
edge weights, then **propagates the signal from the seeds** (BFWalk by
default, or random walk with restart) to score every non-seed gene, and
annotates the ranking with tractability and safety information. The output
is an interactive HTML report with a ranked and annotated list of candidate
targets for the disease.

Known limits (WIP): single-database coverage and curation bias, sparse Open
Targets safety annotation, no systematic validation of the rankings yet, and
no data types beyond pathways and Open Targets associations. For diseases
with many genetic-association genes, the seed-pathway rule can pull a large
share of Reactome into the graph (see `docs/README.md`, Stage 3). It is
meant to complement genetic and clinical evidence for target prioritisation.


### Pipeline, stage by stage

Each stage is a standalone script (`python3 pipeline/<script> --help`);
each reads the previous stage's output from the run directory.

| # | Script | What it does |
|---|---|---|
| 1 | `stage1_ingest.py` | Resolves the disease, then downloads (or reads from `--data-dir`) the version-pinned Open Targets + Reactome sources. Canonicalises every gene identifier to an Ensembl gene ID and writes the normalized `genes` / `gene_sets` / `interactions` / `ot_disease_subset` tables plus `manifest.json`. |
| 2 | `stage2_gsea_discovery.py` | Runs FDR-corrected GSEA of the disease's genetic-association genes against Reactome gene sets, over size-capped and redundancy-collapsed sets. Writes `disease_pathways.tsv`, the disease-relevant pathway list, and `disease_seed_genes.tsv`, the seed genes (the same genetic-association genes). |
| 3 | `stage3_build_graph.py` | Builds a directed, signed gene–gene graph over all genes of the disease-relevant pathways plus all genes of every pathway containing a seed gene, with size-down-weighted co-membership edges plus Reactome functional-interaction edges. Writes the sparse `graph_weight.npz` / `graph_sign.npz` / `graph_gene_index.json` / `graph_metadata.json`. |
| 4 | `stage4_genetic_evidence_weights.py` | Scores each graph gene from non-pathway-derived Open Targets datatypes and re-derives the graph's edge weights from those scores, preserving Stage 3's sparsity. Overwrites `graph_weight.npz`, extends `graph_metadata.json`, and writes `gene_weights.tsv`. |
| 5 | `stage5_score_candidates.py` | Propagates from the seed genes over the graph and ranks every non-seed graph gene — [BFWalk](https://github.com/jedrzejkubica/BFWalk) by default, random-walk-with-restart opt-in (`--method bfwalk`, `rwr` or `bfwalk,rwr`; at least one). Writes `candidate_scores.tsv`. |
| 6 | `stage6_annotate_context.py` | Joins an Open Targets tractability bucket and a safety flag onto each candidate and computes a configurable weighted `composite_score` (`--weights`). Writes `candidates_annotated.tsv`; uses no tissue-expression data. |
| 7 | `run_pipeline.py` + `stage7_report.py` | `run_pipeline.py` chains Stages 1–6 for a disease (subprocess per stage, `--from-manifest` to skip Stage 1). `stage7_report.py` turns the finished run directory into one self-contained interactive HTML report — sortable ranked table, click-to-expand evidence trace for the top 100 candidates (pathways and interactions shared with seed genes). |

**Run pipeline**, one disease → an interactive HTML report:

```
python3 pipeline/run_pipeline.py \
    --disease MONDO_0005265 \
    --data-dir data/ --out-dir out/ --report
```

On the first run Stage 1 downloads the version-pinned Open Targets 26.06 +
Reactome 97 sources into `--data-dir` (a few GB; cached afterwards). The
report is written to `out/report.html`.


**Manual stage-by-stage run** (or use `run_pipeline.py` to automatically run all stages):

#### Stage 1: Ingesting data
```
python3 pipeline/stage1_ingest.py --disease MONDO_0005265 --data-dir data/ --out-dir out/
```

#### Stage 2: Running GSEA
```
python3 pipeline/stage2_gsea_discovery.py           --manifest out/manifest.json
```

#### Stage 3: Building a graph
```
python3 pipeline/stage3_build_graph.py              --manifest out/manifest.json
```

#### Stage 4: Calculating graph weights
```
python3 pipeline/stage4_genetic_evidence_weights.py --manifest out/manifest.json
```

#### Stage 5: Scoring candidates
```
python3 pipeline/stage5_score_candidates.py         --graph-dir out/
```
Seeds are read from `out/disease_seed_genes.tsv` (Stage 2) by default; pass `--seed-genes my_seeds.tsv` (a TSV with a `gene_id` column) to use your own. Scores with BFWalk by default; `--method rwr` scores with random walk with restart instead, and `--method bfwalk,rwr` with both (BFWalk then ranks the output).

Hub bias, measured on inflammatory bowel disease (`--disease MONDO_0005265`, Open Targets 26.06 +
Reactome 97, Stage 4 graph, 9,127 candidates from 502 seeds) as the rank correlation of each score with
the candidate's number of graph neighbours and of Reactome pathways:

| Method | Spearman ρ with node degree | Spearman ρ with no. of Reactome pathways |
|---|---|---|
| `rwr` | 0.74 | 0.62 |
| `bfwalk` | 0.13 | 0.13 |

**Which to use:** use `bfwalk` (the default) to find candidates beyond well-connected hub genes, and
add `rwr` when you prefer well-connected, well-studied genes close to the seeds, e.g. to recover
established targets.

#### Stage 6: Annotating candidates
```
python3 pipeline/stage6_annotate_context.py         --manifest out/manifest.json
```

#### Stage 7: HTML reporting
```
python3 pipeline/stage7_report.py                   --run-dir out/
```

### More detail

- **[`docs/README.md`](docs/README.md)** — per-stage reference: inputs, all
  flags, expected output, and each stage's unit tests.
- **[`development/README.md`](development/README.md)** — how this pipeline
  was built (implementation prompts, session log, a transferable skill).
- **`validation/`** — gitignored; a place for local validation data. Nothing
  in the pipeline reads it. To check whether known disease genes are
  recovered, hold them out of the seeds with Stage 2's
  `--benchmark-holdout-file` and see where they rank among the candidates
  (worked example in `docs/README.md`, Stage 2).
