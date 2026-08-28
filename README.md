# Targets from Pathways

> The scientific idea was introduced by Jędrzej Kubica, Polina Rusina, Siddharth Sethi and Elvis Poku-Adusei at the Open Targets Hackathon on October 21-22, 2025. The technology was created by agentic AI provided by Claude Sonnet 4.6. The original project is available under the MIT license at [https://github.com/jedrzejkubica/targets-from-pathways](https://github.com/jedrzejkubica/targets-from-pathways)


## Quick start

**Install:** `pip install -r requirements.txt` (Python 3.9+, CPU only).

**Demo** — a tiny synthetic Stage 1 output is provided:

```
python3 example_data/build_gsea_example.py
```

```
python3 pipeline/run_pipeline.py --from-manifest example_data/example_run/manifest.json --report
```

**Tests:** `python3 -m pytest` — ~100 offline tests to check if the pipeline is set up correctly.


## Overview

When a drug's molecular target stops working for a patient, clinicians need
alternative targets for treatments. This pipeline provides an orthogonal line
of evidence for that search based on pathway topology alone, using public data: Open Targets
(disease–target genetic associations, tractability, safety) and Reactome
(pathway gene sets and functional interactions). Given a
target and a disease, it discovers disease-relevant pathways, builds a
pathway-based gene–gene graph over the target's and disease's pathway
neighbourhoods, applied genetic evidence into the graph
as node and edge weights, then scores every other gene by topological
and random-walk proximity to the target, and annotates the ranking with
tractability and safety information. The output is provided in an interactive HTML
report with a ranked and annotated list of candidate targets.

Known limits (WIP): single-database coverage
and curation bias, sparse Open Targets safety annotation, a benchmark too
small (5–15 cases), and no other data types beyon pathways. It is meant to complement genetic and clinical
evidence for target prioritisation.


### Pipeline, stage by stage

Each stage is a standalone script (`python3 pipeline/<script> --help`);
each reads the previous stage's output from the run directory.

| # | Script | What it does |
|---|---|---|
| 1 | `stage1_ingest.py` | Resolves the target gene and disease, then downloads (or reads from `--data-dir`) the version-pinned Open Targets + Reactome sources. Canonicalises every gene identifier to an Ensembl gene ID and writes the normalized `genes` / `gene_sets` / `interactions` / `ot_disease_subset` tables plus `manifest.json`. |
| 2 | `stage2_gsea_discovery.py` | Runs FDR-corrected GSEA of the disease's genetic-association genes against Reactome gene sets, over size-capped and redundancy-collapsed sets. Writes `disease_pathways.tsv`, the disease-relevant pathway list. |
| 3 | `stage3_build_graph.py` | Builds a directed, signed gene–gene graph over the union of disease-relevant and target-containing pathways, with size-down-weighted co-membership edges plus Reactome functional-interaction edges. Writes the sparse `graph_weight.npz` / `graph_sign.npz` / `graph_gene_index.json` / `graph_metadata.json`. |
| 4 | `stage4_genetic_evidence_weights.py` | Scores each graph gene from non-pathway-derived Open Targets datatypes and re-derives the graph's edge weights from those scores, preserving Stage 3's sparsity. Overwrites `graph_weight.npz`, extends `graph_metadata.json`, and writes `gene_weights.tsv`. |
| 5 | `stage5_score_candidates.py` | Ranks every non-target graph gene by proximity to the target — a deterministic topology score always, random-walk-with-restart opt-in via `--method topology,rwr`. Writes `candidate_scores.tsv`. |
| 6 | `stage6_annotate_context.py` | Joins an Open Targets tractability bucket and a safety flag onto each candidate and computes a configurable weighted `composite_score` (`--weights`). Writes `candidates_annotated.tsv`; uses no tissue-expression data. |
| 7 | `run_pipeline.py` + `stage7_report.py` | `run_pipeline.py` chains Stages 1–6 for a (target, disease) (subprocess per stage, `--from-manifest` to skip Stage 1). `stage7_report.py` turns the finished run directory into one self-contained interactive HTML report — sortable ranked table, click-to-expand evidence trace per candidate, static benchmark panel. |

**Run pipeline**, one target + disease → an interactive HTML report:

```
python3 pipeline/run_pipeline.py \
    --target NOD2 --disease MONDO_0005265 \
    --data-dir data/ --out-dir out/ --report
```

On the first run Stage 1 downloads the version-pinned Open Targets 26.06 +
Reactome 97 sources into `--data-dir` (a few GB; cached afterwards). The
report is written to `out/report.html`.


**Manual stage-by-stage run** (or use `run_pipeline.py` to automatically run all stages):

#### Stage 1: Ingesting data
```
python3 pipeline/stage1_ingest.py --target NOD2 --disease MONDO_0005265 --data-dir data/ --out-dir out/
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
python3 pipeline/stage5_score_candidates.py         --graph-dir out/ --method topology,rwr
```

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
- **[`benchmarking/README.md`](benchmarking/README.md)** — the
  benchmark-validation module (a separate driver over the pipeline; WIP).
