"""Refit a catalogue disease's GSEA without a held-out gene, for leave_one_out.py --refit-gsea.

A catalogue graph is built from the pathways GSEA found with every seed in the ranking, the held-out seed
included. Holding a seed out of the propagation only, on that graph, leaks it: the gene helped choose the
pathways, and so the nodes, the graph is made of. The refit removes the gene from the GSEA ranking too, selects
the pathways again exactly as build_graph_catalogue.py does (Reactome sets, BH over Reactome, es > 0, FDR below
the cutoff) and rebuilds the graph from them. A gene whose pathways are no longer selected is not in the refitted
graph: leave_one_out.py ranks it with score 0, as an unreachable gene.

GSEA runs as in DisPathCatalogue v2 (pipeline/gsea_core.py): the background is every gene of the catalogue's
merged library, with symbols normalised to Open Targets approved symbols as the catalogue did, and only the
Reactome sets are tested (a set's p-value does not depend on the other sets tested). The permutation null of the
full ranking is reused for every refit (gsea_core null_cache): removing one gene from a ranking of thousands
barely changes it, and sampling it again would make each refit minutes instead of a second.
"""

import logging

import numpy as np
import pandas as pd
import pyarrow.compute as pc
import pyarrow.dataset as ds

import gsea_core
from stage1_ingest import resolve_symbol

logger = logging.getLogger("leave_one_out.gsea_refit")

REACTOME_SUFFIX = " [Reactome]"


def symbol_normaliser(target_paths):
    """{alias: approved symbol} from Open Targets' target table: synonyms and obsolete symbols that are not
    themselves an approved symbol and point to a single gene (ambiguous aliases are left as they are), as
    DisPathCatalogue normalised its library."""
    table = ds.dataset(sorted(map(str, target_paths)), format="parquet").to_table(
        columns=["approvedSymbol", "symbolSynonyms", "obsoleteSymbols"]).to_pandas()
    table = table.dropna(subset=["approvedSymbol"])
    approved = set(table["approvedSymbol"])
    aliases = {}
    for column in ("symbolSynonyms", "obsoleteSymbols"):
        for symbol, entries in zip(table["approvedSymbol"], table[column]):
            if entries is None:
                continue
            for entry in entries:
                label = entry["label"]
                if label not in approved:
                    aliases.setdefault(label, set()).add(symbol)
    return approved, {alias: next(iter(symbols)) for alias, symbols in aliases.items() if len(symbols) == 1}


def read_library(gmt_path, approved, aliases):
    """The merged GMT as {set name: normalised symbols, first occurrence order}."""
    library = {}
    for line in gmt_path.read_text().splitlines():
        fields = line.split("\t")
        seen = {}
        for symbol in fields[2:]:
            if symbol:
                seen.setdefault(symbol if symbol in approved else aliases.get(symbol, symbol), None)
        library[fields[0]] = list(seen)
    return library


def read_disease_scores(scores_path, disease):
    """{symbol: genetic_association score} of one disease (genetic_gene_scores_input.parquet)."""
    table = ds.dataset(str(scores_path)).to_table(filter=pc.field("diseaseId") == disease,
                                                columns=["gene", "score"]).to_pandas()
    return dict(zip(table["gene"], table["score"].astype(float)))


class CatalogueGsea:
    """The library shared by every disease: the background order and the Reactome sets to test."""

    def __init__(self, gmt_path, target_paths, order_seed=0):
        approved, aliases = symbol_normaliser(target_paths)
        library = read_library(gmt_path, approved, aliases)
        self.background = gsea_core.ordered_background(library, order_seed)
        self.reactome = {name: genes for name, genes in library.items() if name.endswith(REACTOME_SUFFIX)}
        logger.info("GSEA library: %d sets (%d Reactome), %d background genes.", len(library), len(self.reactome),
                    len(self.background))


class DiseaseRefit:
    """One disease's GSEA, refittable without given genes. The full ranking's null is sampled once."""

    def __init__(self, library, scores, fdr_cutoff, seed=0):
        self.library, self.scores, self.fdr_cutoff, self.seed = library, scores, fdr_cutoff, seed
        self.null_cache = {}
        self.full_selection = self.select(scores)

    def select(self, scores):
        """The pathways selected from scores: Reactome sets with es > 0 and BH FDR (over the Reactome sets
        tested) below the cutoff, as their GMT names."""
        result = gsea_core.run_gsea(scores, self.library.reactome, background=self.library.background, seed=self.seed,
                                    null_cache=self.null_cache)
        return frozenset(gsea_core.significant_sets(result, self.fdr_cutoff)["set_id"])

    def selection_without(self, symbols):
        if not set(symbols) & set(self.scores):
            return self.full_selection
        return self.select({gene: score for gene, score in self.scores.items() if gene not in symbols})


def symbols_by_gene(scores, index):
    """{Ensembl gene ID: the disease's scored symbols mapping to it} (the catalogue's seed mapping)."""
    by_gene = {}
    for symbol in scores:
        gene_id = resolve_symbol(symbol, index)[0]
        if gene_id is not None:
            by_gene.setdefault(gene_id, set()).add(symbol)
    return by_gene


def jaccard(a, b):
    a, b = set(a), set(b)
    return float(len(a & b) / len(a | b)) if a | b else np.nan


def selection_summary(full, refitted):
    """How a refit changed the pathway selection."""
    return {"n_pathways": len(refitted), "n_pathways_lost": len(full - refitted),
            "n_pathways_gained": len(refitted - full), "pathway_jaccard": jaccard(full, refitted)}


def catalogue_selection(disease_dir, gene_sets):
    """The catalogue graph's pathways as GMT names, mapped through gene_sets' reactome_id."""
    pathways = pd.read_csv(disease_dir / "disease_pathways.tsv", sep="\t")
    names = gene_sets.drop_duplicates("reactome_id").set_index("reactome_id")["set_id"]
    return frozenset(names.reindex(pathways["set_id"]).dropna())
