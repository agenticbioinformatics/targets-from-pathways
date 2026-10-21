# Targets from Pathways

> The scientific idea was introduced by Jędrzej Kubica, Polina Rusina, Siddharth Sethi and Elvis Poku-Adusei at the Open Targets Hackathon on October 21-22, 2025. The technology was created by agentic AI provided by Claude Sonnet 4.6 and Opus 5.5. The original project is available under the MIT license at [https://github.com/jedrzejkubica/targets-from-pathways](https://github.com/jedrzejkubica/targets-from-pathways)


## Quick start

**Install:** `pip install -r requirements.txt`

> **Note:** [BFWalk](https://github.com/jedrzejkubica/BFWalk) (for Stage 5) is still in development and is not on PyPI yet - `requirements.txt` installs it from TestPyPI. It will be available via `pip install bfwalk` shortly. To install it on its own: `pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ bfwalk`. For installation details see [BFWalk](https://github.com/jedrzejkubica/BFWalk).

**Demo** - a tiny example input (30 genes, 18 seeds) is provided:

```
python example_data/build_gsea_example.py
```

```
python pipeline/run_pipeline.py --from-manifest example_data/example_run/manifest.json --report
```

Expected: the log ends with `Scoring 3 candidates from 16 seeds` and
`Report written (3 candidates)`; `example_data/example_run/candidates_annotated.tsv`
ranks the 3 non-seed genes on the GSEA pathways: `ENSG00000000022` and `ENSG00000000024`
(one step from a seed) above `ENSG00000000023` (two steps).
Genes on pathways GSEA did not detect (`ENSG...025`-`033`) are not in the graph.

Run `python -m pytest` to check if the pipeline is set up correctly (~135 tests).


## Overview

Given a disease, this pipeline proposes candidate drug targets based on
pathway topology, using public data: Open Targets (disease-target genetic
associations, tractability, safety) and Reactome (pathway gene sets and
functional interactions). The disease's genes with an Open Targets
`genetic_association` score > 0 are the **seed genes**. The pipeline discovers
disease-relevant pathways (GSEA of the seed genes), builds an unweighted
gene-gene graph of the Reactome functional interactions between all genes
on those pathways, then **propagates the signal from the seeds** (BFWalk by
default, or random walk with restart) to score every non-seed gene, and
annotates the ranking with tractability and safety information. The output
is an interactive HTML report with a ranked and annotated list of candidate
targets for the disease.


### Pipeline, stage by stage

Each stage is a standalone script (`python pipeline/<script> --help`), each reads the previous stage's output.

| # | Stage | Description |
|---|---|---|
| 1 | `stage1_ingest.py` | Downloads (or reads from `--data-dir`) the version-pinned Open Targets + Reactome sources. Maps every gene name to an Ensembl ID and writes `genes` / `gene_sets` / `interactions` / `ot_disease_subset` tables plus `manifest.json`. |
| 2 | `stage2_gsea_discovery.py` | Runs FDR-corrected GSEA of the disease-associated genes against Reactome gene sets. Writes `disease_pathways.tsv` and `disease_seed_genes.tsv`. |
| 3 | `stage3_build_graph.py` | Builds a directed, signed, unweighted gene-gene graph: its nodes are all genes of the disease-relevant pathways found by GSEA, its edges the Reactome functional interactions between them. Writes `graph_weight.npz` / `graph_sign.npz` / `graph_gene_index.json` / `graph_metadata.json`. |
| 4 | `stage4_genetic_evidence_weights.py` | Optionally weights the graph with your own node and/or edge weights (`--node-weights`, `--edge-weights`; without them the graph stays unweighted), and scores each graph gene from its Open Targets `genetic_association` evidence for Stage 6's annotation. Rewrites the graph files, extends `graph_metadata.json`, and writes `gene_weights.tsv`. |
| 5 | `stage5_score_candidates.py` | Propagates from the seeds over the graph and ranks every non-seed gene - [BFWalk](https://github.com/jedrzejkubica/BFWalk) by default, or Random Walk with Restart (RWR). Use `--method bfwalk`, `rwr` or `bfwalk,rwr`. Writes `candidate_scores.tsv`. |
| 6 | `stage6_annotate_context.py` | Joins genetic evidence, and tractability and safety from Open Targets, onto candidates for reading; candidates keep Stage 5's ranking. Writes `candidates_annotated.tsv`. |
| 7 | `run_pipeline.py` + `stage7_report.py` | `run_pipeline.py` chains Stages 1-6 (subprocess per stage, use `--from-manifest` to skip Stage 1). `stage7_report.py` turns the finished run directory into an interactive HTML report for the top 100 candidates. |

**Run pipeline**, one disease to an interactive HTML report:

```
python pipeline/run_pipeline.py \
    --disease MONDO_0005265 \
    --data-dir data/ \
    --out-dir out/ \
    --report
```

On the first run Stage 1 downloads the version-pinned Open Targets 26.06 +
Reactome 97 sources into `--data-dir` (a few GB; cached afterwards). The
report is written to `out/report.html`.


**Run pipeline stage-by-stage**:

#### Stage 1: Ingesting data
```
python pipeline/stage1_ingest.py --disease MONDO_0005265 --data-dir data/ --out-dir out/
```

#### Stage 2: Running GSEA
```
python pipeline/stage2_gsea_discovery.py --manifest out/manifest.json
```

#### Stage 3: Building a graph
```
python pipeline/stage3_build_graph.py --manifest out/manifest.json
```

#### Stage 4: Adding your own weights (optional)
```
python pipeline/stage4_genetic_evidence_weights.py --manifest out/manifest.json \
    --node-weights node_weights.tsv --edge-weights edge_weights.tsv
```
Both files are optional TSVs with Ensembl gene IDs and positive weights (`gene_id`, `weight` and `gene_a`, `gene_b`, `weight`) that weight the graph for Stage 5's propagation (also accepted by `run_pipeline.py`); see [docs/README.md](docs/README.md#running-stage-4-pipelinestage4_genetic_evidence_weightspy) for the details.

#### Stage 5: Scoring candidates
```
python pipeline/stage5_score_candidates.py --graph-dir out/
```
Seeds are read from `out/disease_seed_genes.tsv` (Stage 2) by default; use `--seed-genes seeds.tsv` (a TSV with a `gene_id` column) to use your own seeds. Scores with BFWalk by default; `--method rwr` scores with random walk with restart instead, and `--method bfwalk,rwr` scores with both methods.

> **Note:** Hub bias: correlation of each score with the candidates degrees are **0.13 for BFWalk** and **0.74 for RWR** (measured on inflammatory bowel disease `--disease MONDO_0005265`, Open Targets 26.06 + Reactome 97, on an earlier genetic-evidence-weighted graph; not yet re-measured on the current unweighted graph). We suggest to use BFWalk to find candidates beyond well-connected hub genes, and add RWR when you prefer well-connected, well-studied genes.

#### Stage 6: Annotating candidates
```
python pipeline/stage6_annotate_context.py --manifest out/manifest.json
```

#### Stage 7: HTML reporting
```
python pipeline/stage7_report.py --run-dir out/
```

### More detail

- **[`docs/README.md`](docs/README.md)** - per-stage reference file: inputs, all
  flags, expected output, and each stage's unit tests.
- **[`development/README.md`](development/README.md)** - how this pipeline
  was built (implementation prompts, session log, a transferable skill).
