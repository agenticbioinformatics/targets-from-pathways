"""pytest tests for run-directory integrity across stages, on fresh synthetic runs (tests/fixtures/synthetic_run.py).

- ``find_artifact`` (stage0_schemas.py): a stage reads the Stage 1 files next to its --manifest, even when
  the manifest lists a path that exists elsewhere (a copied or moved run directory), and refuses a file
  whose sha256 no longer matches the manifest.
- Small diseases (Stages 2 + 3): too few seed genes for GSEA, or seed genes on no disease pathway, stop
  Stage 3 with a message instead of a traceback.
- Benchmark holdout (Stages 2 + 4): a held-out gene is not a seed, and its own Open Targets evidence does
  not feed Stage 4's scores - otherwise it would be ranked by the very association the benchmark is meant
  to recover.
- Stage 4 leaves Stage 3's unweighted graph unchanged without user weights, applies --node-weights /
  --edge-weights when given, and a re-run without them restores the unweighted graph.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest
import scipy.sparse as sp

import stage2_gsea_discovery as gsea
import stage3_build_graph as bg
import stage4_weights as wts
from fixtures.synthetic_run import write_synthetic_stage1
from stage0_schemas import find_artifact, load_manifest, sha256_of_file

HOLDOUT_FILE = Path(__file__).parent / "fixtures" / "benchmark_holdout_synthetic.txt"
HELD_OUT = "ENSG00000000017"  # OLDG17 in benchmark_holdout_synthetic.txt
# G5: unlike G17 it stays in the graph when held out (its GSEA pathways remain significant without it).
HELD_OUT_IN_GRAPH = "ENSG00000000005"


def _copy_run(tmp_path, listed_dir=None):
    """Write a synthetic Stage 1 run to tmp_path/run. listed_dir: rewrite the manifest's artifact paths
    to point there instead (as in a run directory copied from elsewhere)."""
    run_dir = write_synthetic_stage1(tmp_path / "run").parent
    if listed_dir is not None:
        manifest = json.loads((run_dir / "manifest.json").read_text())
        for artifact in manifest["output_artifacts"]:
            artifact["path"] = str(listed_dir / artifact["path"].split("/")[-1])
        (run_dir / "manifest.json").write_text(json.dumps(manifest))
    return run_dir


def _replace_artifact(run_dir, name, df):
    """Rewrite a Stage 1 artifact and record its new sha256, as if Stage 1 had written it."""
    path = run_dir / name
    df.to_parquet(path, index=False)
    manifest = json.loads((run_dir / "manifest.json").read_text())
    for artifact in manifest["output_artifacts"]:
        if artifact["path"].split("/")[-1] == name:
            artifact["sha256"] = sha256_of_file(path)
    (run_dir / "manifest.json").write_text(json.dumps(manifest))


# --------------------------------------------------------------------------
# find_artifact
# --------------------------------------------------------------------------


def test_copied_run_reads_its_own_files_not_the_listed_ones(tmp_path):
    # The manifest lists the original run's files, which exist - the copy next to the manifest must win.
    original = write_synthetic_stage1(tmp_path / "original").parent
    run_dir = _copy_run(tmp_path, listed_dir=original)
    manifest_path = run_dir / "manifest.json"
    assert find_artifact(load_manifest(manifest_path), manifest_path, "genes.parquet") == run_dir / "genes.parquet"


def test_listed_path_is_the_fallback(tmp_path):
    original = write_synthetic_stage1(tmp_path / "original").parent
    run_dir = _copy_run(tmp_path, listed_dir=original)
    (run_dir / "genes.parquet").unlink()
    manifest_path = run_dir / "manifest.json"
    assert find_artifact(load_manifest(manifest_path), manifest_path, "genes.parquet") == original / "genes.parquet"


def test_artifact_changed_after_stage1_is_refused(tmp_path):
    run_dir = _copy_run(tmp_path)
    ot = pd.read_parquet(run_dir / "ot_disease_subset.parquet")
    ot.iloc[:-1].to_parquet(run_dir / "ot_disease_subset.parquet", index=False)
    manifest_path = run_dir / "manifest.json"
    with pytest.raises(SystemExit):
        find_artifact(load_manifest(manifest_path), manifest_path, "ot_disease_subset.parquet")


# --------------------------------------------------------------------------
# benchmark holdout
# --------------------------------------------------------------------------


def _run_stages_2_to_4(run_dir, holdout_file=HOLDOUT_FILE):
    manifest = str(run_dir / "manifest.json")
    stage2_argv = ["--manifest", manifest]
    if holdout_file is not None:
        stage2_argv += ["--benchmark-holdout-file", str(holdout_file)]
    gsea.run(gsea.parse_args(stage2_argv))
    bg.run(bg.parse_args(["--manifest", manifest]))
    wts.run(wts.parse_args(["--manifest", manifest]))


def test_holdout_gene_evidence_does_not_feed_stage4(tmp_path):
    run_dir = _copy_run(tmp_path)
    # Give the held-out gene strong genetic evidence: a leak would now be large.
    ot = pd.read_parquet(run_dir / "ot_disease_subset.parquet")
    ot.loc[ot["gene_id"] == HELD_OUT_IN_GRAPH, "score"] = 0.9
    _replace_artifact(run_dir, "ot_disease_subset.parquet", ot)
    holdout_file = tmp_path / "holdout.txt"
    holdout_file.write_text(HELD_OUT_IN_GRAPH + "\n")

    _run_stages_2_to_4(run_dir, holdout_file)

    assert HELD_OUT_IN_GRAPH not in set(pd.read_csv(run_dir / "disease_seed_genes.tsv", sep="\t")["gene_id"])
    assert list(pd.read_csv(run_dir / "benchmark_holdout_genes.tsv", sep="\t")["gene_id"]) == [HELD_OUT_IN_GRAPH]

    metadata = json.loads((run_dir / "graph_metadata.json").read_text())
    assert metadata["genetic_evidence_score"][HELD_OUT_IN_GRAPH] == 0.0
    assert metadata["stage4_genetic_evidence"]["holdout_genes"] == [HELD_OUT_IN_GRAPH]



def test_rerun_without_holdout_clears_the_holdout_file(tmp_path):
    run_dir = _copy_run(tmp_path)
    _run_stages_2_to_4(run_dir)
    _run_stages_2_to_4(run_dir, holdout_file=None)
    assert pd.read_csv(run_dir / "benchmark_holdout_genes.tsv", sep="\t").empty
    metadata = json.loads((run_dir / "graph_metadata.json").read_text())
    assert metadata["genetic_evidence_score"][HELD_OUT] > 0.0


def test_stage4_user_weights_and_rerun_without_them(tmp_path):
    run_dir = _copy_run(tmp_path)
    manifest = str(run_dir / "manifest.json")
    gsea.run(gsea.parse_args(["--manifest", manifest]))
    bg.run(bg.parse_args(["--manifest", manifest]))
    stage3_weight = sp.load_npz(run_dir / "graph_weight.npz")
    assert stage3_weight.nnz > 0 and set(stage3_weight.data) == {1.0}

    # Without weight files the graph is unchanged.
    wts.run(wts.parse_args(["--manifest", manifest]))
    assert abs(stage3_weight - sp.load_npz(run_dir / "graph_weight.npz")).max() == 0
    assert json.loads((run_dir / "graph_metadata.json").read_text())["stage4_user_weights"]["weighted"] is False

    # G6 -> N24 and N24 -> N23 are graph edges; the node weight of N23 halves N24 -> N23 (mean of 1 and 0).
    index = json.loads((run_dir / "graph_gene_index.json").read_text())
    g6, n23, n24 = "ENSG00000000006", "ENSG00000000023", "ENSG00000000024"
    (tmp_path / "edges.tsv").write_text(f"gene_a\tgene_b\tweight\n{n24}\t{g6}\t0.25\n")
    (tmp_path / "nodes.tsv").write_text(f"gene_id\tweight\n{n23}\t0.5\n")
    wts.run(wts.parse_args(["--manifest", manifest, "--edge-weights", str(tmp_path / "edges.tsv"),
                            "--node-weights", str(tmp_path / "nodes.tsv")]))
    weighted = sp.load_npz(run_dir / "graph_weight.npz").toarray()
    assert weighted[index.index(g6), index.index(n24)] == pytest.approx(0.25)
    assert weighted[index.index(n24), index.index(n23)] == pytest.approx(0.75)
    assert (weighted > 0).sum() == stage3_weight.nnz  # no edge added or removed
    stats = json.loads((run_dir / "graph_metadata.json").read_text())["stage4_user_weights"]
    assert stats["weighted"] is True and stats["n_edge_weights_used"] == 1 and stats["n_node_weights_used"] == 1

    # A re-run without the files restores the unweighted graph: weights never compound.
    wts.run(wts.parse_args(["--manifest", manifest]))
    assert abs(stage3_weight - sp.load_npz(run_dir / "graph_weight.npz")).max() == 0


# --------------------------------------------------------------------------
# diseases with few genetic-association genes
# --------------------------------------------------------------------------


def test_disease_too_small_for_gsea_stops_stage3_with_a_message(tmp_path, caplog):
    # Three seed genes: GSEA cannot test any pathway, so disease_pathways.tsv is empty and there is no graph.
    run_dir = _copy_run(tmp_path)
    ot = pd.read_parquet(run_dir / "ot_disease_subset.parquet")
    genetic = ot["datatype_id"] == "genetic_association"
    keep = set(ot.loc[genetic, "gene_id"].sort_values().head(3))
    _replace_artifact(run_dir, "ot_disease_subset.parquet", ot[~genetic | ot["gene_id"].isin(keep)])

    manifest = str(run_dir / "manifest.json")
    gsea.run(gsea.parse_args(["--manifest", manifest]))
    disease_pathways = pd.read_csv(run_dir / "disease_pathways.tsv", sep="\t")
    assert disease_pathways.empty and list(disease_pathways.columns) == ["set_id", "gene_set", "source_db", "pval", "fdr"]
    with pytest.raises(SystemExit) as excinfo:
        bg.run(bg.parse_args(["--manifest", manifest]))
    assert excinfo.value.code != 0
    assert "lists no disease-relevant pathway" in caplog.text


def test_genes_with_zero_genetic_score_are_not_seeds(tmp_path):
    run_dir = _copy_run(tmp_path)
    ot = pd.read_parquet(run_dir / "ot_disease_subset.parquet")
    genetic = ot["datatype_id"] == "genetic_association"
    zeroed = ot.loc[genetic, "gene_id"].sort_values().iloc[0]
    ot.loc[genetic & (ot["gene_id"] == zeroed), "score"] = 0.0
    _replace_artifact(run_dir, "ot_disease_subset.parquet", ot)

    gsea.run(gsea.parse_args(["--manifest", str(run_dir / "manifest.json")]))
    seeds = set(pd.read_csv(run_dir / "disease_seed_genes.tsv", sep="\t")["gene_id"])
    assert zeroed not in seeds
    assert seeds == set(ot.loc[genetic & (ot["score"] > 0), "gene_id"])


def test_no_seed_gene_on_a_disease_pathway_stops_stage3_with_a_message(tmp_path, caplog):
    # A disease whose only seed gene is in no Reactome pathway (e.g. AURKC for spermatogenic failure 5) gave an
    # empty graph and a networkx traceback.
    run_dir = _copy_run(tmp_path)
    gsea.run(gsea.parse_args(["--manifest", str(run_dir / "manifest.json")]))
    seeds = run_dir / "no_pathway_seed.tsv"
    seeds.write_text("gene_id\nENSG99999999999\n")
    with pytest.raises(SystemExit) as excinfo:
        bg.run(bg.parse_args(["--manifest", str(run_dir / "manifest.json"), "--seed-genes", str(seeds)]))
    assert excinfo.value.code != 0
    assert "is on a disease-relevant pathway" in caplog.text
