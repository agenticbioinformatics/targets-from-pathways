"""Preranked GSEA with a gene-set permutation null: the DisPathCatalogue v2 method, in numpy.

Stage 2 and the validation scripts (build_graph_catalogue.py, leave_one_out.py) call this module, so the
pipeline and the catalogue its validation graphs come from use one algorithm.

The method (DisPathCatalogue docs/GSEA_PARAMETERS.md, v2):
- Background: every gene in the library. A gene without a score for the disease scores 0.
- Ties: most genes tie at 0, so the background is put in a seeded random order (numpy default_rng(order_seed))
  before a stable sort by score. Alphabetical order would bias ES.
- ES: blitzGSEA's enrichment score on the mean-centred signature: running sum of |score| over hits, minus a
  constant step over misses, taking the extreme value. Gene sets keep only background genes; sets with fewer
  than min_size or more than max_size genes are not tested.
- p-values: per gene-set size k, the ES of random k-gene sets on the same signature (blitzGSEA's O(k) sparse
  sampler). 1,000 draws to start, x10 while a positive-ES set of that size has fewer than 10 null draws as
  extreme, up to min(100,000, 2e7 / k). p = (1 + #null >= es) / (1 + #null > 0) for es > 0, mirrored for es < 0;
  nes = es / mean of the same-sign null.
- FDR: Benjamini-Hochberg over all tested sets, both ES signs.

blitzGSEA's own gamma-fit p-values are not used: on a zero-inflated genetic signature they are miscalibrated.
The ES and the sampler are ports of jdhayhurst/blitzgsea v0.0.1 (Apache-2.0), which the catalogue ran with;
with the same inputs this module reproduces the catalogue's es, pval and fdr (tests/test_gsea_core.py).
"""

import numpy as np
import pandas as pd

MIN_SIZE = 5
MAX_SIZE = 4000
PERM_MIN = 1_000
PERM_MAX = 100_000
PERM_BUDGET = 20_000_000
MIN_EXCEED = 10
RESULT_COLUMNS = ["set_id", "es", "nes", "pval", "n_perm", "fdr", "geneset_size"]


def ordered_background(library, order_seed=0):
    """The background (all library genes) in the seeded random order the catalogue used."""
    background = sorted({gene for genes in library.values() for gene in genes})
    return list(np.random.default_rng(order_seed).permutation(background))


def build_signature(scores, background):
    """Genes sorted by score, descending; ties keep the background order. Returns (genes, values)."""
    values = np.array([float(scores.get(gene, 0.0)) for gene in background])
    order = np.argsort(-values, kind="stable")
    return [background[i] for i in order], values[order]


def enrichment_score(abs_signature, hits):
    """blitzGSEA's ES for one gene set: hits are indices into the sorted signature."""
    hit = np.zeros(len(abs_signature))
    hit[hits] = 1
    norm_hit = 1.0 / abs_signature[hits].sum()
    norm_miss = 1.0 / (len(abs_signature) - len(hits))
    running = np.cumsum(hit * abs_signature * norm_hit - (1 - hit) * norm_miss)
    return running[np.argmax(np.abs(running))]


def null_es(abs_signature, k, n, seed):
    """n null ES values for random k-gene sets: port of blitzgsea.get_peak_size_adv (float32, same RNG stream)."""
    rng = np.random.default_rng(seed)
    N = len(abs_signature)
    norm_miss = np.float32(1.0 / (N - k))
    abs_f32 = abs_signature.astype(np.float32)
    k_idx = np.arange(k, dtype=np.float32)
    batch_size = max(1, min(n, int(50_000_000 // (k * 24))))
    chunks = []
    remaining = n
    while remaining > 0:
        batch = min(batch_size, remaining)
        remaining -= batch
        # Sorted sampling without replacement from exponential spacings (Devroye 1986).
        cum = np.cumsum(rng.exponential(1.0, size=(batch, k + 1)).astype(np.float32), axis=1)
        hs = np.clip(np.floor(cum[:, :k] / cum[:, -1:] * N).astype(np.int32), 0, N - 1)
        for _ in range(8):
            bad = (np.diff(hs, axis=1) == 0).any(axis=1)
            if not bad.any():
                break
            c2 = np.cumsum(rng.exponential(1.0, size=(int(bad.sum()), k + 1)).astype(np.float32), axis=1)
            hs[bad] = np.clip(np.floor((c2[:, :k] / c2[:, -1:]) * N).astype(np.int32), 0, N - 1)
        ah = abs_f32[hs]
        ca = np.cumsum(ah, axis=1)
        nh = np.float32(1.0) / ca[:, -1:]
        at_hit = ca * nh - (hs.astype(np.float32) - k_idx) * norm_miss
        before_hit = at_hit - ah * nh
        rows = np.arange(batch)
        i_at, i_bf = np.abs(at_hit).argmax(axis=1), np.abs(before_hit).argmax(axis=1)
        v_at, v_bf = np.abs(at_hit)[rows, i_at], np.abs(before_hit)[rows, i_bf]
        pos_at = hs[rows, i_at].astype(np.int32)
        pos_bf = np.maximum(0, hs[rows, i_bf].astype(np.int32) - 1)
        at_wins = (v_at > v_bf) | ((v_at == v_bf) & (pos_at <= pos_bf))
        es = np.where(at_wins, at_hit[rows, i_at], before_hit[rows, i_bf])
        chunks.append(es[~np.isnan(es)])
    return np.concatenate(chunks).astype(np.float64)


def permutation_stats(es, sizes, abs_signature, seed=0, null_cache=None):
    """Permutation p-values, NES and draw counts for every tested set (see the module docstring).

    null_cache: an optional dict {set size: sorted null ES}. A size already in it is not sampled again (and not
    grown); a size sampled here is added. Leave-one-out uses this to reuse the full signature's null for a
    signature with one gene removed, where the null barely changes and sampling it again dominates the run time.
    """
    pval = np.full(es.shape, np.nan)
    nes = np.full(es.shape, np.nan)
    n_perm = np.zeros(es.shape, dtype=np.int64)
    for k in np.unique(sizes):
        idx = np.flatnonzero(sizes == k)
        obs = es[idx]
        cached = None if null_cache is None else null_cache.get(int(k))
        cap = max(PERM_MIN, min(PERM_MAX, PERM_BUDGET // int(k)))
        chunks, done, target = [], 0, PERM_MIN
        while cached is None:
            # Same seed scheme as the catalogue, so its p-values are reproduced exactly.
            chunks.append(np.asarray(null_es(abs_signature, int(k), target - done, seed + int(k) * 7919 + done)))
            done = target
            null = np.sort(np.concatenate(chunks))
            exceed = np.where(obs > 0, len(null) - np.searchsorted(null, obs, side="left"),
                              np.searchsorted(null, obs, side="right"))
            pos_exceed = exceed[obs > 0]
            if done >= cap or not len(pos_exceed) or pos_exceed.min() >= MIN_EXCEED:
                break
            target = min(cap, target * 10)
        if cached is not None:
            null = cached
            exceed = np.where(obs > 0, len(null) - np.searchsorted(null, obs, side="left"),
                              np.searchsorted(null, obs, side="right"))
        elif null_cache is not None:
            null_cache[int(k)] = null
        n_pos, n_neg = int((null > 0).sum()), int((null < 0).sum())
        mean_pos = null[null > 0].mean() if n_pos else np.nan
        mean_neg = -null[null < 0].mean() if n_neg else np.nan
        denom = np.where(obs > 0, n_pos, n_neg)
        with np.errstate(invalid="ignore", divide="ignore"):
            pval[idx] = np.where(denom > 0, (1 + exceed) / (1 + denom), np.nan)
            nes[idx] = obs / np.where(obs > 0, mean_pos, mean_neg)
        n_perm[idx] = len(null)
    return pval, nes, n_perm


def bh_fdr(pvals):
    """Benjamini-Hochberg q-values over the non-NaN p-values (NaN stays NaN)."""
    q = np.full(pvals.shape, np.nan)
    ok = ~np.isnan(pvals)
    p = pvals[ok]
    if not len(p):
        return q
    order = np.argsort(p, kind="stable")
    ranked = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty(len(p))
    out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q


def run_gsea(scores, library, background=None, seed=0, order_seed=0, min_size=MIN_SIZE, max_size=MAX_SIZE,
             null_cache=None):
    """GSEA of one disease. scores: {gene: score} (genes absent score 0); library: {set_id: genes}.

    background: the gene order to break ties with; defaults to ordered_background(library, order_seed). Pass it
    in when running many signatures against one library. null_cache: see permutation_stats. Returns one row per
    tested set (both ES signs), in library order, with RESULT_COLUMNS.
    """
    if background is None:
        background = ordered_background(library, order_seed)
    genes, values = build_signature(scores, background)
    abs_signature = np.abs(values - values.mean())
    position = {gene: i for i, gene in enumerate(genes)}
    set_ids, es, sizes = [], [], []
    for set_id, members in library.items():
        hits = sorted({position[g] for g in members if g in position})
        # A set holding the whole background has no misses: its ES and null are undefined.
        if not min_size <= len(hits) <= max_size or len(hits) >= len(abs_signature):
            continue
        if abs_signature[hits].sum() == 0:
            continue
        set_ids.append(set_id)
        es.append(enrichment_score(abs_signature, hits))
        sizes.append(len(hits))
    if not set_ids:
        return pd.DataFrame({c: pd.Series(dtype="float64") for c in RESULT_COLUMNS}).astype({"set_id": "object"})
    es, sizes = np.array(es), np.array(sizes)
    pval, nes, n_perm = permutation_stats(es, sizes, abs_signature, seed=seed, null_cache=null_cache)
    return pd.DataFrame({"set_id": set_ids, "es": es, "nes": nes, "pval": pval, "n_perm": n_perm,
                         "fdr": bh_fdr(pval), "geneset_size": sizes})


def significant_sets(result, fdr_cutoff):
    """Sets enriched at the top of the ranking: es > 0 and fdr below the cutoff, most significant first."""
    keep = result[(result["es"] > 0) & (result["fdr"] < fdr_cutoff)]
    return keep.sort_values(["fdr", "pval", "set_id"]).reset_index(drop=True)
