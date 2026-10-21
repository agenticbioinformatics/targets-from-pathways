"""pytest tests for validation/build_graph_catalogue.py on tiny synthetic inputs in the catalogue's formats.

Genes come from the Stage 1 fixtures (tests/fixtures/opentargets/target.parquet: TP53, PTGS2 with synonym
COX2, BRCA1; the FI fixture also has an unmappable FAKEGENE123 and a predicted row).

- D1: Reactome "Pathway One" {TP53, COX2 -> PTGS2, BRCA1} significant; "Pathway Two" not (FDR 0.2); a
  negative-ES Reactome row counts in BH's m but is never selected; a GO row is ignored.
- D2: nothing significant -> left out of the catalogue.
- D3: "Pathway Two" {TP53, BRCA1} significant, but its only seed (PTGS2) is not on it.
"""

from __future__ import annotations

import importlib.util
import shutil
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

FIXTURES = Path(__file__).parent / "fixtures"
SCRIPT = Path(__file__).parent.parent / "validation" / "build_graph_catalogue.py"
spec = importlib.util.spec_from_file_location("build_graph_catalogue", SCRIPT)
gc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gc)

TP53, PTGS2, BRCA1 = "ENSG00000141510", "ENSG00000073756", "ENSG00000012048"
P1, P2, P3 = "Pathway One [Reactome]", "Pathway Two [Reactome]", "Pathway Three [Reactome]"


def _row(disease, pathway, source, es, pval, size=3):
    return {"diseaseId": disease, "pathway": pathway, "source": source, "es": es, "nes": es * 2, "pval": pval,
            "geneset_size": size}


@pytest.fixture(scope="module")
def catalogue_run(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("catalogue")
    cat_dir = tmp / "catalogue"
    cat_dir.mkdir()
    pd.DataFrame([
        _row("D1", P1, "Reactome", 0.5, 0.001),
        _row("D1", P2, "Reactome", 0.3, 0.2),
        _row("D1", P3, "Reactome", -0.4, 0.001),
        _row("D1", "Some term (GO:1) [GO BP]", "GO BP", 0.9, 1e-6),
        _row("D2", P1, "Reactome", 0.4, 0.5),
        _row("D3", P2, "Reactome", 0.5, 0.0001),
        _row("D3", P1, "Reactome", -0.1, None),  # null p-value: not tested
    ]).to_parquet(cat_dir / "genetic_disease_pathway_catalogue.parquet", index=False)
    (cat_dir / "merged_gene_sets.gmt").write_text(
        f"{P1}\tR-HSA-1\tTP53\tCOX2\tBRCA1\n"
        f"{P2}\tR-HSA-2\tTP53\tBRCA1\n"
        "Some term (GO:1) [GO BP]\tGO:1\tTP53\tPTGS2\n"
    )
    pd.DataFrame([
        {"diseaseId": "D1", "gene": "TP53", "score": 0.8},
        {"diseaseId": "D1", "gene": "BRCA1", "score": 0.2},
        {"diseaseId": "D1", "gene": "FAKEGENE123", "score": 0.5},
        {"diseaseId": "D2", "gene": "TP53", "score": 0.9},
        {"diseaseId": "D3", "gene": "PTGS2", "score": 0.9},
    ]).to_parquet(cat_dir / "genetic_gene_scores_input.parquet", index=False)

    target_dir = tmp / "target"
    target_dir.mkdir()
    shutil.copy(FIXTURES / "opentargets" / "target.parquet", target_dir / "part-0.parquet")
    with zipfile.ZipFile(FIXTURES / "reactome" / "FIsInGene_test_with_annotations.txt.zip") as archive:
        archive.extractall(tmp)
    fi_file = tmp / "FIsInGene_test_with_annotations.txt"

    inputs = ["--catalogue-dir", str(cat_dir), "--target-dir", str(target_dir), "--fi-file", str(fi_file)]

    def build(name, *options):
        out_dir = tmp / name
        gc.run(gc.parse_args(inputs + ["--out-dir", str(out_dir), *options]))
        return {
            "edges": pd.read_parquet(out_dir / "graph_catalogue_edges.parquet"),
            "nodes": pd.read_parquet(out_dir / "graph_catalogue_nodes.parquet"),
            "diseases": pd.read_parquet(out_dir / "graph_catalogue_diseases.parquet").set_index("diseaseId"),
        }

    result = build("all")
    result["build"] = build
    return result


def test_benjamini_hochberg_matches_statsmodels():
    from statsmodels.stats.multitest import multipletests

    pvals = [0.01, 0.04, 0.03, 0.2, 0.001, 0.5]
    # Hand check on the first four: p * m / rank = 0.04, 0.0533, 0.06, 0.2, then the step-up minimum.
    assert gc.benjamini_hochberg(pvals[:4]) == pytest.approx([0.04, 0.16 / 3, 0.16 / 3, 0.2])
    assert gc.benjamini_hochberg(pvals) == pytest.approx(multipletests(pvals, method="fdr_bh")[1])


def test_only_diseases_with_a_significant_reactome_pathway_are_included(catalogue_run):
    diseases = catalogue_run["diseases"]
    assert sorted(diseases.index) == ["D1", "D3"]  # D2's only pathway has FDR 0.5
    assert diseases.loc["D1", "n_pathways"] == 1    # Pathway Two (FDR 0.2) and the negative-ES row are not selected
    assert list(diseases.loc["D1", "reactome_ids"]) == ["R-HSA-1"]


def test_graph_is_stage3s_curated_interactions_between_pathway_genes(catalogue_run):
    edges = catalogue_run["edges"]
    d1 = set(map(tuple, edges[edges["diseaseId"] == "D1"][["gene_a", "gene_b", "sign"]].values.tolist()))
    # TP53 -> COX2 (PTGS2) activates; TP53 - BRCA1 undirected (both directions); COX2 -| BRCA1 inhibits.
    # The predicted and FAKEGENE123 rows are dropped.
    assert d1 == {(TP53, PTGS2, 1), (TP53, BRCA1, 0), (BRCA1, TP53, 0), (PTGS2, BRCA1, -1)}
    d3 = set(map(tuple, edges[edges["diseaseId"] == "D3"][["gene_a", "gene_b"]].values.tolist()))
    assert d3 == {(TP53, BRCA1), (BRCA1, TP53)}


def test_seed_scores_and_unmapped_symbols(catalogue_run):
    nodes = catalogue_run["nodes"]
    d1 = dict(zip(nodes[nodes["diseaseId"] == "D1"]["gene_id"], nodes[nodes["diseaseId"] == "D1"]["seed_score"]))
    assert d1 == {TP53: 0.8, BRCA1: 0.2, PTGS2: 0.0}
    diseases = catalogue_run["diseases"]
    assert diseases.loc["D1", "n_seeds_in_graph"] == 2
    assert diseases.loc["D1", "n_seed_symbols_unmapped"] == 1  # FAKEGENE123
    # D3's only seed (PTGS2) is not on its pathway: the disease is kept but flagged by n_seeds_in_graph.
    assert diseases.loc["D3", "n_seeds"] == 1 and diseases.loc["D3", "n_seeds_in_graph"] == 0
    assert np.all(nodes["seed_score"] >= 0)


def test_n_diseases_and_diseases_options(catalogue_run):
    build = catalogue_run["build"]
    assert list(build("first", "--n-diseases", "1")["diseases"].index) == ["D1"]
    # Comma-separated IDs; D2 has no significant pathway and NOPE is not in the catalogue.
    chosen = build("chosen", "--diseases", "D3, D2,NOPE")
    assert list(chosen["diseases"].index) == ["D3"]
    assert set(chosen["edges"]["diseaseId"]) == {"D3"}
