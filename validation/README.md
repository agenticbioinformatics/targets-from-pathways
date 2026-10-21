# Validation

Does the pipeline recover disease's known genes, and does it do so better than simply sorting the Open Targets data by score?

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
