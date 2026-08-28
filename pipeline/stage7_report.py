"""Stage 7 — Output/report.

Reads a run directory holding the Stage 1-6 artifacts and writes ONE
self-contained, interactive HTML file (no server, no external assets, no
JS/CSS dependencies). Open it in a browser or hand it to a teammate.

Run directory contents used (all produced by earlier stages):

    manifest.json               Stage 1  — target/disease/versions
    genes.parquet               Stage 1  — Ensembl -> HGNC symbol (display only)
    gene_sets.parquet           Stage 1  — pathway membership (evidence trace)
    interactions.parquet        Stage 1  — signed edges (evidence trace)
    ot_disease_subset.parquet   Stage 1  — OT datatypes/scores (evidence trace)
    disease_pathways.tsv        Stage 2  — the disease-relevant pathway list
    graph_gene_index.json       Stage 3  — which genes are in the graph
    graph_metadata.json         Stage 3/4 — seed, which datatypes Stage 4 admitted
    candidates_annotated.tsv    Stage 6  — the ranked table + composite trace

Genes are shown by **HGNC symbol**; the Ensembl gene ID stays the key
(its own column, and the anchor id for each candidate's detail block).
Symbols are resolved from ``genes.parquet`` here, at render time only.

The benchmark validation panel is *static* — it embeds
``benchmarking/benchmark_summary.txt`` (or, if absent,
``benchmarking/benchmark_summary.example.txt``, clearly labelled). This
stage never runs the benchmark; see ``benchmarking/``.
"""

from __future__ import annotations

import argparse
import html
import json
import logging
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from stage0_schemas import AnnotatedCandidatesSchema

PIPELINE_DIR = Path(__file__).resolve().parent
REPO_DIR = PIPELINE_DIR.parent

logger = logging.getLogger("stage7_report")

_SIGN_LABEL = {1: "activating (+)", -1: "inhibiting (−)", 0: "unsigned (0)"}
_DEFAULT_ADMITTED_DATATYPES = ["genetic_association", "known_drug"]


def _fail(msg: str) -> None:
    logger.error(msg)
    sys.exit(1)


# ==========================================================================
# Load
# ==========================================================================


def load_run(run_dir: Path) -> dict:
    run_dir = Path(run_dir)

    def _need(name: str) -> Path:
        p = run_dir / name
        if not p.exists():
            _fail(f"{name} missing from run directory {run_dir} — has the pipeline finished?")
        return p

    manifest = json.loads(_need("manifest.json").read_text())
    graph_meta = json.loads(_need("graph_metadata.json").read_text())
    graph_index = set(json.loads(_need("graph_gene_index.json").read_text()))

    annotated = pd.read_csv(_need("candidates_annotated.tsv"), sep="\t")
    AnnotatedCandidatesSchema.validate(annotated)

    genes = pd.read_parquet(_need("genes.parquet"))
    symbols = dict(zip(genes["gene_id"], genes["symbol"]))

    admitted = list(
        graph_meta.get("stage4_genetic_evidence", {}).get("datatypes", _DEFAULT_ADMITTED_DATATYPES)
    )

    return {
        "run_dir": run_dir,
        "manifest": manifest,
        "graph_meta": graph_meta,
        "graph_index": graph_index,
        "annotated": annotated,
        "symbols": symbols,
        "admitted_datatypes": set(admitted),
        "gene_sets": pd.read_parquet(_need("gene_sets.parquet")),
        "interactions": pd.read_parquet(_need("interactions.parquet")),
        "ot": pd.read_parquet(_need("ot_disease_subset.parquet")),
        "disease_pathways": pd.read_csv(_need("disease_pathways.tsv"), sep="\t"),
        "stage4_ran": "stage4_genetic_evidence" in graph_meta,
    }


def _sym(gene_id: str, symbols: dict) -> str:
    return symbols.get(gene_id, gene_id)


# ==========================================================================
# Evidence trace (assembled at render time from Stage 1-3 artifacts)
# ==========================================================================


def evidence_trace(gene_id: str, run: dict) -> dict:
    target_id = run["manifest"]["resolved_target"]["gene_id"]
    gs = run["gene_sets"]
    symbols = run["symbols"]

    set_name = dict(zip(gs["set_id"], gs["set_name"]))
    target_sets = set(gs.loc[gs["gene_id"] == target_id, "set_id"])
    disease_sets = set(run["disease_pathways"]["set_id"])
    union_sets = disease_sets | target_sets  # the Stage 3 pathway union
    cand_sets = set(gs.loc[gs["gene_id"] == gene_id, "set_id"]) & union_sets

    shared = sorted((set_name.get(s, s) for s in cand_sets & target_sets), key=str.lower)
    member_only = sorted((set_name.get(s, s) for s in cand_sets - target_sets), key=str.lower)

    ot = run["ot"]
    rows = ot[ot["gene_id"] == gene_id]
    datatypes = []
    for dt, grp in rows.groupby("datatype_id"):
        datatypes.append(
            {
                "datatype": str(dt),
                "score": float(grp["score"].max()),
                "datasources": sorted(str(x) for x in grp["datasource_id"].unique()),
                # "feeds" only means something when Stage 4 actually ran
                "feeds": run["stage4_ran"] and str(dt) in run["admitted_datatypes"],
            }
        )
    datatypes.sort(key=lambda d: (not d["feeds"], -d["score"]))

    it = run["interactions"]
    touching = it[(it["gene_a"] == gene_id) | (it["gene_b"] == gene_id)]
    interactions = []
    for r in touching.itertuples(index=False):
        other = r.gene_b if r.gene_a == gene_id else r.gene_a
        if other not in run["graph_index"]:
            continue  # dropped by Stage 3's pathway-pool scoping — not in the graph
        interactions.append(
            {
                "a": _sym(r.gene_a, symbols),
                "b": _sym(r.gene_b, symbols),
                "directed": bool(r.directed),
                "sign": _SIGN_LABEL.get(int(r.sign), str(r.sign)),
                "evidence_type": str(r.evidence_type),
                "confidence": None if pd.isna(r.confidence) else float(r.confidence),
                "to_target": other == target_id,
            }
        )

    return {"shared_pathways": shared, "member_pathways": member_only,
            "datatypes": datatypes, "interactions": interactions}


# ==========================================================================
# Benchmark summary (static)
# ==========================================================================


def benchmark_summary(benchmark_dir: Path | None = None) -> tuple[str | None, str]:
    benchmark_dir = Path(benchmark_dir) if benchmark_dir else (REPO_DIR / "benchmarking")
    real = benchmark_dir / "benchmark_summary.txt"
    example = benchmark_dir / "benchmark_summary.example.txt"
    if real.exists():
        return real.read_text().strip(), "benchmark_summary.txt"
    if example.exists():
        return example.read_text().strip(), "benchmark_summary.example.txt (EXAMPLE — not a real validation run)"
    return None, "none"


# ==========================================================================
# HTML
# ==========================================================================

_CSS = """
:root { color-scheme: light dark; }
body { font: 15px/1.5 system-ui, -apple-system, Segoe UI, Roboto, sans-serif; margin: 2rem auto;
       max-width: 60rem; padding: 0 1rem; }
h1 { font-size: 1.4rem; margin-bottom: .2rem; }
.meta { color: #666; font-size: .88rem; margin-top: 0; }
.lead { color: #444; font-size: .92rem; }
table { border-collapse: collapse; width: 100%; margin: 1rem 0; font-size: .9rem; }
th, td { text-align: right; padding: .35rem .55rem; border-bottom: 1px solid #ddd; }
th.txt, td.txt { text-align: left; }
th { cursor: pointer; user-select: none; border-bottom: 2px solid #999; white-space: nowrap; }
th:hover { color: #06c; }
tbody tr:nth-child(odd) { background: rgba(127,127,127,.06); }
.ens { color: #888; font-family: ui-monospace, monospace; font-size: .82rem; }
a { color: #06c; text-decoration: none; } a:hover { text-decoration: underline; }
details { border: 1px solid #ddd; border-radius: 6px; margin: .5rem 0; padding: .3rem .7rem; }
details[open] { background: rgba(127,127,127,.05); }
summary { cursor: pointer; }
.trace .k { font-weight: 600; margin: .6rem 0 .1rem; font-size: .86rem; text-transform: uppercase;
            letter-spacing: .03em; color: #777; }
.trace ul { margin: .1rem 0 .4rem 1.2rem; padding: 0; }
.trace code { font-size: .82rem; background: rgba(127,127,127,.12); padding: .05rem .3rem; border-radius: 3px; }
.badge { font-size: .72rem; background: #06c; color: #fff; padding: .05rem .35rem; border-radius: 3px;
         margin-left: .3rem; }
.badge.tgt { background: #c60; }
.src { color: #888; font-size: .82rem; }
.none { color: #999; font-style: italic; }
pre { background: rgba(127,127,127,.1); padding: .8rem; border-radius: 6px; overflow-x: auto;
      font-size: .84rem; }
.warn { color: #b30; font-size: .86rem; }
""".strip()

_JS = """
function sortTable(i, numeric){
  var tb = document.getElementById('ranktbody');
  var rows = Array.prototype.slice.call(tb.rows);
  var key = 'c' + i;
  var asc = tb.getAttribute('data-key') === key ? tb.getAttribute('data-dir') !== 'asc' : true;
  rows.sort(function(a, b){
    var x = a.cells[i].getAttribute('data-v'); if (x === null) x = a.cells[i].textContent;
    var y = b.cells[i].getAttribute('data-v'); if (y === null) y = b.cells[i].textContent;
    if (numeric){ x = parseFloat(x); y = parseFloat(y); }
    return (x < y ? -1 : x > y ? 1 : 0) * (asc ? 1 : -1);
  });
  rows.forEach(function(r){ tb.appendChild(r); });
  tb.setAttribute('data-key', key); tb.setAttribute('data-dir', asc ? 'asc' : 'desc');
}
function openGene(id){
  var d = document.getElementById(id);
  if (d){ d.open = true; d.scrollIntoView({behavior: 'smooth', block: 'start'}); }
}
""".strip()


def _e(x) -> str:
    return html.escape(str(x))


def _trace_html(gene_id: str, sym: str, row, tr: dict) -> str:
    """``row`` is a namedtuple from ``df.itertuples()``."""
    parts = [
        f'<details id="g-{_e(gene_id)}">',
        f"<summary><b>{_e(sym)}</b> <span class=ens>{_e(gene_id)}</span> "
        f"&mdash; composite <b>{row.composite_score:.3f}</b></summary>",
        '<div class="trace">',
        '<p class="k">Composite breakdown</p>',
        f"<p><code>{_e(row.composite_breakdown)}</code><br>"
        f"weights <code>{_e(row.composite_weights)}</code></p>",
    ]

    parts.append('<p class="k">Shared pathways with the target</p>')
    if tr["shared_pathways"]:
        parts.append("<ul>" + "".join(f"<li>{_e(n)}</li>" for n in tr["shared_pathways"]) + "</ul>")
    else:
        parts.append('<p class="none">none in the Stage 3 pathway union</p>')

    if tr["member_pathways"]:
        parts.append('<p class="k">Also a member of (disease-relevant; target is not)</p>')
        parts.append("<ul>" + "".join(f"<li>{_e(n)}</li>" for n in tr["member_pathways"]) + "</ul>")

    parts.append('<p class="k">Open Targets evidence datatypes</p>')
    if tr["datatypes"]:
        lis = []
        for d in tr["datatypes"]:
            badge = '<span class="badge">feeds genetic-evidence score</span>' if d["feeds"] else ""
            src = f' <span class="src">({", ".join(_e(s) for s in d["datasources"])})</span>'
            lis.append(f"<li>{_e(d['datatype'])} &mdash; {d['score']:.3f}{badge}{src}</li>")
        parts.append("<ul>" + "".join(lis) + "</ul>")
    else:
        parts.append('<p class="none">no Open Targets association rows for this gene in this disease</p>')

    parts.append('<p class="k">Interactions in the graph</p>')
    if tr["interactions"]:
        lis = []
        for it in tr["interactions"]:
            arrow = "&rarr;" if it["directed"] else "&mdash;"
            conf = "" if it["confidence"] is None else f", confidence {it['confidence']:.2f}"
            tgt = '<span class="badge tgt">to target</span>' if it["to_target"] else ""
            lis.append(
                f"<li>{_e(it['a'])} {arrow} {_e(it['b'])} &mdash; {_e(it['sign'])}, "
                f"{_e(it['evidence_type'])}{conf}{tgt}</li>"
            )
        parts.append("<ul>" + "".join(lis) + "</ul>")
    else:
        parts.append('<p class="none">no signed interaction edges touch this gene</p>')

    parts.append("</div></details>")
    return "".join(parts)


def render_html(run: dict, bench_text: str | None, bench_label: str) -> str:
    m = run["manifest"]
    df = run["annotated"].reset_index(drop=True)
    symbols = run["symbols"]
    tgt = m["resolved_target"]
    dis = m["disease"]
    versions = ", ".join(f"{s['db']} {s['version']}" for s in m.get("sources", []))
    graph_kind = "Stage 4 (genetic-evidence-weighted)" if run["stage4_ran"] else "Stage 3 (structural)"

    have_rwr = "rwr_score" in df.columns
    have_gen = "genetic_evidence_score" in df.columns

    # ---- ranked table ----
    # (name, numeric-sort, left-align-text)
    headers = [
        ("#", True, False), ("symbol", False, True), ("gene", False, True),
        ("composite", True, False), ("topology", True, False),
    ]
    if have_rwr:
        headers.append(("rwr", True, False))
    if have_gen:
        headers.append(("genetic ev.", True, False))
    headers += [("tractability", False, True), ("safety", False, True)]

    thead = "".join(
        f'<th class="txt" onclick="sortTable({i},false)">{_e(name)}</th>'
        if txt else f'<th onclick="sortTable({i},{str(num).lower()})">{_e(name)}</th>'
        for i, (name, num, txt) in enumerate(headers)
    )

    body_rows = []
    for rank, row in enumerate(df.itertuples(index=False), start=1):
        gid = row.gene
        sym = _sym(gid, symbols)
        cells = [
            f'<td data-v="{rank}">{rank}</td>',
            f'<td class="txt"><a href="#" onclick="openGene(\'g-{_e(gid)}\');return false;">{_e(sym)}</a></td>',
            f'<td class="txt ens">{_e(gid)}</td>',
            f'<td data-v="{row.composite_score}">{row.composite_score:.3f}</td>',
            f'<td data-v="{row.topology_score}">{row.topology_score:.3f}</td>',
        ]
        if have_rwr:
            cells.append(f'<td data-v="{row.rwr_score}">{row.rwr_score:.4f}</td>')
        if have_gen:
            cells.append(f'<td data-v="{row.genetic_evidence_score}">{row.genetic_evidence_score:.3f}</td>')
        cells.append(f'<td class="txt">{_e(row.tractability)}</td>')
        cells.append(f'<td class="txt">{_e(row.safety)}</td>')
        body_rows.append("<tr>" + "".join(cells) + "</tr>")

    # ---- per-candidate detail blocks (rank order) ----
    details = [
        _trace_html(row.gene, _sym(row.gene, symbols), row, evidence_trace(row.gene, run))
        for row in df.itertuples(index=False)
    ]

    # ---- validation panel ----
    if bench_text is None:
        bench_block = (
            '<p class="none">No benchmark summary found. Run '
            "<code>benchmarking/benchmark_validate.py</code> and re-generate this report.</p>"
        )
    else:
        warn = ""
        if "EXAMPLE" in bench_label:
            warn = f'<p class="warn">Source: {_e(bench_label)}</p>'
        else:
            warn = f'<p class="meta">Source: {_e(bench_label)}</p>'
        bench_block = warn + f"<pre>{_e(bench_text)}</pre>"

    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    tbody_html = "\n".join(body_rows)
    details_html = "\n".join(details)

    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Target Rescue via Pathways &mdash; {_e(_sym(tgt['gene_id'], symbols))} / {_e(dis['name'])}</title>
<style>{_CSS}</style>
</head><body>
<h1>Target Rescue via Pathways &mdash; candidate report</h1>
<p class="meta">
Target <b>{_e(_sym(tgt['gene_id'], symbols))}</b> <span class=ens>{_e(tgt['gene_id'])}</span>
(input <code>{_e(tgt['input'])}</code>) &middot;
Disease <b>{_e(dis['name'])}</b> <span class=ens>{_e(dis['efo_id'])}</span> &middot;
{len(df)} candidates &middot; seed {_e(m.get('seed'))} &middot; {_e(versions)} &middot;
graph: {_e(graph_kind)} &middot; generated {generated}
</p>
<p class="lead">
<code>composite_score</code> is a weighted average of the normalised components; each
candidate's <i>breakdown</i> below shows every component's weighted contribution
(the terms sum to the score). Click a header to sort; click a gene symbol to open
its evidence trace.
</p>

<table id="ranktable">
<thead><tr>{thead}</tr></thead>
<tbody id="ranktbody">
{tbody_html}
</tbody>
</table>

<h2>Evidence trace</h2>
{details_html}

<h2>Benchmark validation</h2>
<p class="meta">Static &mdash; embedded from <code>benchmarking/</code>, not recomputed for this report.</p>
{bench_block}

<script>{_JS}</script>
</body></html>
"""


# ==========================================================================
# Orchestration
# ==========================================================================


def build_report(run_dir: Path, out_path: Path, *, benchmark_dir: Path | None = None) -> Path:
    run = load_run(run_dir)
    text, label = benchmark_summary(benchmark_dir)
    html_str = render_html(run, text, label)
    out_path = Path(out_path)
    out_path.write_text(html_str)
    logger.info(
        "Wrote %s (%d candidates, benchmark: %s).",
        out_path, len(run["annotated"]), label,
    )
    return out_path


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="stage7_report.py",
        description="Render a self-contained interactive HTML report from a finished pipeline run.",
    )
    p.add_argument(
        "--run-dir",
        required=True,
        type=Path,
        help="Directory with the Stage 1-6 artifacts (manifest.json, candidates_annotated.tsv, ...).",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output HTML path. Default: <run-dir>/report.html.",
    )
    p.add_argument(
        "--benchmark-dir",
        type=Path,
        default=None,
        help="Directory holding benchmark_summary.txt / .example.txt. Default: <repo>/benchmarking.",
    )
    return p


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    args = build_arg_parser().parse_args(argv)
    out = args.out or (args.run_dir / "report.html")
    build_report(args.run_dir, out, benchmark_dir=args.benchmark_dir)


if __name__ == "__main__":
    main()
