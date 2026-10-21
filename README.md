# Targets from Pathways

> The scientific idea was introduced by Jędrzej Kubica, Polina Rusina, Siddharth Sethi and Elvis Poku-Adusei at the Open Targets Hackathon on October 21-22, 2025. The technology was created by agentic AI provided by Claude Sonnet 4.6 and Opus 5.5. The original project is available under the MIT license at [https://github.com/jedrzejkubica/targets-from-pathways](https://github.com/jedrzejkubica/targets-from-pathways)


## Overview

Given a disease, the pipeline proposes candidate drug targets from pathway topology. The disease's
**seed genes** are its genes with an Open Targets `genetic_association` score > 0. GSEA of the seeds
finds the disease-relevant Reactome pathways, and the pipeline builds an unweighted graph of the
Reactome functional interactions between all genes on those pathways. It then propagates the
seeds' signal over the graph (BFWalk by default, or random walk with restart) to rank every
non-seed gene, and annotates the ranking with Open Targets tractability and safety. The result is
an interactive HTML report.


## Installation

```
pip install -r requirements.txt
```

[BFWalk](https://github.com/jedrzejkubica/BFWalk), Stage 5's default method, is not on PyPI yet:
`requirements.txt` installs it from TestPyPI, and building it needs gcc, zlib-devel and OpenMP.
Without it, run Stage 5 with `--method rwr`.

Check the installation:

```
python pipeline/run_pipeline.py --help
```

It imports every stage first, so if it prints the help, every package is installed; otherwise it lists
what is missing. `python -m pytest` runs the test suite (no network or downloads needed).


## Usage

**Whole pipeline**, one disease to `out/report.html`:

```
python pipeline/run_pipeline.py --disease MONDO_0005265 --data-dir data/ --out-dir out/ --report
```

**Stage by stage** - each stage reads the previous stage's files in `out/` (`--help` lists all options):

| Stage | What it does | Command |
|---|---|---|
| 1 | Downloads the pinned Open Targets + Reactome data and maps every gene to an Ensembl ID. | `python pipeline/stage1_ingest.py --disease MONDO_0005265 --data-dir data/ --out-dir out/` |
| 2 | Runs GSEA of the seed genes against the Reactome pathways. | `python pipeline/stage2_gsea_discovery.py --manifest out/manifest.json` |
| 3 | Builds the graph: the functional interactions between the genes of the GSEA pathways. | `python pipeline/stage3_build_graph.py --manifest out/manifest.json` |
| 4 | Optionally weights the graph with your own node and/or edge weights (`--node-weights`, `--edge-weights`; both optional TSV files). | `python pipeline/stage4_weights.py --manifest out/manifest.json` |
| 5 | Propagates from the seeds and ranks every other graph gene (`--method bfwalk`, `rwr` or `bfwalk,rwr`). | `python pipeline/stage5_score_candidates.py --graph-dir out/` |
| 6 | Adds Open Targets tractability and safety to the candidates. | `python pipeline/stage6_annotate_context.py --manifest out/manifest.json` |
| 7 | Writes the HTML report; click a gene to see its evidence. | `python pipeline/stage7_report.py --run-dir out/` |

BFWalk is less biased towards highly connected hub genes than RWR: use it to find candidates beyond
well-studied hubs, and add RWR when you prefer well-connected genes close to the seeds.

See [`docs/README.md`](docs/README.md) for each stage's inputs, outputs and options.


## Data

Stage 1 downloads the pinned versions into `--data-dir` on the first run (about 1.5 GB), checks
their sha256, and reuses them afterwards:

- [Open Targets Platform](https://platform.opentargets.org/) 26.06: `association_by_datasource_direct`
  (disease-gene evidence), `target` (gene identifiers, tractability, safety) and `disease`.
- [Reactome](https://reactome.org/) release 97: `ReactomePathways.gmt` (pathway gene sets) and
  `Ensembl2Reactome.txt`; [Reactome functional interactions](https://reactome.org/tools/reactome-fiviz)
  `FIsInGene_04142025_with_annotations.txt`.

Every stage writes into the run directory (`out/`): parquet tables from Stage 1, TSV tables
(`disease_pathways.tsv`, `candidate_scores.tsv`, `candidates_annotated.tsv`, ...), the graph as sparse
matrices (`graph_weight.npz`, `graph_sign.npz`) with `graph_gene_index.json`, and `report.html`.
Genes are always Ensembl gene IDs.

Development history: [`development/README.md`](development/README.md). Validation scripts:
[`validation/README.md`](validation/README.md).
