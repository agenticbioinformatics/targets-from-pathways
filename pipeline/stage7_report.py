"""Stage 7 - HTML report.

Turns a run directory into one self-contained interactive HTML file (no server, no external assets):
a sortable ranked table of every candidate and an evidence trace for the top candidates (--n-traces).
Only manifest.json and candidates_annotated.tsv are required; a missing optional file only hides the
report sections that need it.
"""

import argparse
import html
import json
import logging
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from stage0_schemas import AnnotatedCandidatesSchema, fail, ranking_column

PIPELINE_DIR = Path(__file__).resolve().parent
REPO_DIR = PIPELINE_DIR.parent

logger = logging.getLogger("stage7_report")

SIGN_LABELS = {1: "activating (+)", -1: "inhibiting (−)", 0: "unsigned (0)"}
# Evidence traces are rendered for the top candidates only: one per candidate made a report of thousands
# of candidates tens of MB and slow to build.
DEFAULT_N_TRACES = 100

CSS = """
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

JAVASCRIPT = """
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


def load_run(run_dir):
    run_dir = Path(run_dir)
    for required in ["manifest.json", "candidates_annotated.tsv"]:
        if not (run_dir / required).exists():
            fail(f"{required} missing from {run_dir} - this is a required Stage 1/6 output; has the pipeline run?")
    manifest = json.loads((run_dir / "manifest.json").read_text())
    annotated = pd.read_csv(run_dir / "candidates_annotated.tsv", sep="\t")
    AnnotatedCandidatesSchema.validate(annotated)

    missing = []
    optional = {}
    for name in ["graph_metadata.json", "graph_gene_index.json", "genes.parquet", "gene_sets.parquet",
                 "interactions.parquet", "ot_disease_subset.parquet", "disease_pathways.tsv",
                 "disease_seed_genes.tsv"]:
        path = run_dir / name
        if not path.exists():
            missing.append(name)
            optional[name] = None
        elif name.endswith(".json"):
            optional[name] = json.loads(path.read_text())
        elif name.endswith(".parquet"):
            optional[name] = pd.read_parquet(path)
        else:
            optional[name] = pd.read_csv(path, sep="\t")

    graph_metadata = optional["graph_metadata.json"] or {}
    graph_index = None
    if optional["graph_gene_index.json"] is not None:
        graph_index = set(optional["graph_gene_index.json"])
    seeds = set()
    if optional["disease_seed_genes.tsv"] is not None:
        seeds = set(optional["disease_seed_genes.tsv"]["gene_id"])
    symbols = {}
    if optional["genes.parquet"] is not None:
        symbols = dict(zip(optional["genes.parquet"]["gene_id"], optional["genes.parquet"]["symbol"]))

    return {
        "run_dir": run_dir,
        "manifest": manifest,
        "graph_meta": graph_metadata,
        "graph_index": graph_index,
        "annotated": annotated,
        "symbols": symbols,
        "gene_sets": optional["gene_sets.parquet"],
        "interactions": optional["interactions.parquet"],
        "ot": optional["ot_disease_subset.parquet"],
        "disease_pathways": optional["disease_pathways.tsv"],
        "seeds": seeds,
        "missing": missing,
    }


def symbol_of(gene_id, symbols):
    return symbols.get(gene_id, gene_id)


def evidence_trace(gene_id, run):
    """Collect why a candidate is near the seed genes: pathways shared with seeds, Open Targets datatypes and
    interactions. Stage 5 keeps only scores, so this is rebuilt here from the Stage 1-3 files."""
    seeds = run["seeds"]
    symbols = run["symbols"]
    gene_sets = run["gene_sets"]
    ot = run["ot"]
    interactions_df = run["interactions"]
    graph_index = run["graph_index"]

    shared_pathways = []
    member_only_pathways = []
    if gene_sets is not None:
        set_names = dict(zip(gene_sets["set_id"], gene_sets["set_name"]))
        seed_rows = gene_sets[gene_sets["gene_id"].isin(seeds)]
        n_seeds_per_set = seed_rows.groupby("set_id")["gene_id"].nunique().to_dict()
        disease_sets = set()
        if run["disease_pathways"] is not None:
            disease_sets = set(run["disease_pathways"]["set_id"])
        # Only pathways that Stage 3 used to build the graph: the disease pathways found by GSEA.
        candidate_sets = set(gene_sets.loc[gene_sets["gene_id"] == gene_id, "set_id"]) & disease_sets
        for set_id in candidate_sets:
            name = set_names.get(set_id, set_id)
            if set_id in n_seeds_per_set:
                shared_pathways.append((name, n_seeds_per_set[set_id]))
            else:
                member_only_pathways.append(name)
        # Pathways with the most seed genes first.
        shared_pathways.sort(key=lambda entry: (-entry[1], entry[0].lower()))
        member_only_pathways.sort(key=str.lower)

    datatypes = []
    if ot is not None:
        for datatype, rows in ot[ot["gene_id"] == gene_id].groupby("datatype_id"):
            datatypes.append({
                "datatype": str(datatype),
                "score": float(rows["score"].max()),
                "datasources": sorted(rows["datasource_id"].astype(str).unique()),
            })
        datatypes.sort(key=lambda entry: (-entry["score"], entry["datatype"]))

    interactions = []
    if interactions_df is not None:
        touching = interactions_df[(interactions_df["gene_a"] == gene_id) | (interactions_df["gene_b"] == gene_id)]
        for row in touching.itertuples(index=False):
            if row.gene_a == gene_id:
                other_gene = row.gene_b
            else:
                other_gene = row.gene_a
            # Without graph_gene_index.json every touching interaction is shown.
            if graph_index is not None and other_gene not in graph_index:
                continue
            confidence = None
            if not pd.isna(row.confidence):
                confidence = float(row.confidence)
            interactions.append({
                "a": symbol_of(row.gene_a, symbols),
                "b": symbol_of(row.gene_b, symbols),
                "directed": bool(row.directed),
                "sign": SIGN_LABELS.get(int(row.sign), str(row.sign)),
                "evidence_type": str(row.evidence_type),
                "confidence": confidence,
                "to_seed": other_gene in seeds,
            })

    return {
        "shared_pathways": shared_pathways,
        "member_pathways": member_only_pathways,
        "datatypes": datatypes,
        "interactions": interactions,
        "pathways_available": gene_sets is not None,
        "datatypes_available": ot is not None,
        "interactions_available": interactions_df is not None,
        "interactions_unscoped": interactions_df is not None and graph_index is None,
    }


def escape(value):
    return html.escape(str(value))


def html_list(items):
    text = "<ul>"
    for item in items:
        text += f"<li>{item}</li>"
    return text + "</ul>"


def trace_html(gene_id, symbol, row, trace, rank_column):
    parts = [
        f'<details id="g-{escape(gene_id)}">',
        f"<summary><b>{escape(symbol)}</b> <span class=ens>{escape(gene_id)}</span> "
        f"&mdash; {escape(rank_column)} <b>{getattr(row, rank_column):.4f}</b></summary>",
        '<div class="trace">',
        '<p class="k">Pathways shared with seed genes</p>',
    ]
    if not trace["pathways_available"]:
        parts.append('<p class="none">gene_sets.parquet not in the run directory</p>')
    elif trace["shared_pathways"]:
        names = []
        for name, n_seeds in trace["shared_pathways"]:
            names.append(f"{escape(name)} <span class=\"src\">({n_seeds} seed gene{'s' if n_seeds != 1 else ''})</span>")
        parts.append(html_list(names))
    else:
        parts.append('<p class="none">none among the disease-relevant pathways</p>')

    if trace["member_pathways"]:
        parts.append('<p class="k">Also a member of (disease-relevant, no seed genes)</p>')
        names = []
        for name in trace["member_pathways"]:
            names.append(escape(name))
        parts.append(html_list(names))

    parts.append('<p class="k">Open Targets evidence datatypes</p>')
    if not trace["datatypes_available"]:
        parts.append('<p class="none">ot_disease_subset.parquet not in the run directory</p>')
    elif trace["datatypes"]:
        items = []
        for datatype in trace["datatypes"]:
            sources = escape(", ".join(datatype["datasources"]))
            items.append(f"{escape(datatype['datatype'])} &mdash; {datatype['score']:.3f} "
                         f"<span class=\"src\">({sources})</span>")
        parts.append(html_list(items))
    else:
        parts.append('<p class="none">no Open Targets association rows for this gene in this disease</p>')

    if trace["interactions_unscoped"]:
        parts.append('<p class="k">Interactions touching this gene</p>')
    else:
        parts.append('<p class="k">Interactions in the graph</p>')
    if not trace["interactions_available"]:
        parts.append('<p class="none">interactions.parquet not in the run directory</p>')
    elif trace["interactions"]:
        items = []
        for interaction in trace["interactions"]:
            arrow = "&mdash;"
            if interaction["directed"]:
                arrow = "&rarr;"
            confidence = ""
            if interaction["confidence"] is not None:
                confidence = f", confidence {interaction['confidence']:.2f}"
            to_seed = ""
            if interaction["to_seed"]:
                to_seed = '<span class="badge tgt">with seed</span>'
            items.append(f"{escape(interaction['a'])} {arrow} {escape(interaction['b'])} &mdash; "
                         f"{escape(interaction['sign'])}, {escape(interaction['evidence_type'])}{confidence}{to_seed}")
        parts.append(html_list(items))
    else:
        parts.append('<p class="none">no signed interaction edges touch this gene</p>')

    parts.append("</div></details>")
    return "".join(parts)


def render_html(run, n_traces=DEFAULT_N_TRACES):
    manifest = run["manifest"]
    candidates = run["annotated"].reset_index(drop=True)
    symbols = run["symbols"]
    disease = manifest["disease"]

    versions = []
    for source in manifest.get("sources", []):
        versions.append(f"{source['db']} {source['version']}")
    if run["graph_meta"].get("stage4_user_weights", {}).get("weighted"):
        graph_kind = "functional interactions on GSEA pathways, with user node/edge weights (Stage 4)"
    elif run["graph_meta"]:
        graph_kind = "Stage 3 (unweighted functional interactions on GSEA pathways)"
    else:
        graph_kind = "unknown (no graph_metadata.json)"

    rank_column = ranking_column(candidates.columns)
    has_rwr = "rwr_score" in candidates.columns
    has_bfwalk = "bfwalk_score" in candidates.columns

    # (column title, sort numerically, left-aligned text)
    headers = [("#", True, False), ("symbol", False, True), ("gene", False, True)]
    if has_rwr:
        headers.append(("rwr", True, False))
    if has_bfwalk:
        headers.append(("bfwalk", True, False))
    headers.append(("tractability", False, True))
    headers.append(("safety", False, True))
    header_html = ""
    for i, (title, numeric, is_text) in enumerate(headers):
        css_class = ""
        if is_text:
            css_class = ' class="txt"'
        header_html += f'<th{css_class} onclick="sortTable({i},{str(numeric).lower()})">{escape(title)}</th>'

    table_rows = []
    details = []
    rank = 0
    for row in candidates.itertuples(index=False):
        rank += 1
        gene_id = row.gene
        symbol = symbol_of(gene_id, symbols)
        has_trace = rank <= n_traces
        if has_trace:
            symbol_cell = (f'<td class="txt"><a href="#" onclick="openGene(\'g-{escape(gene_id)}\');return false;">'
                           f'{escape(symbol)}</a></td>')
        else:
            symbol_cell = f'<td class="txt">{escape(symbol)}</td>'
        cells = [
            f'<td data-v="{rank}">{rank}</td>',
            symbol_cell,
            f'<td class="txt ens">{escape(gene_id)}</td>',
        ]
        if has_rwr:
            cells.append(f'<td data-v="{row.rwr_score}">{row.rwr_score:.4f}</td>')
        if has_bfwalk:
            cells.append(f'<td data-v="{row.bfwalk_score}">{row.bfwalk_score:.4f}</td>')
        cells.append(f'<td class="txt">{escape(row.tractability)}</td>')
        cells.append(f'<td class="txt">{escape(row.safety)}</td>')
        table_rows.append("<tr>" + "".join(cells) + "</tr>")
        if has_trace:
            details.append(trace_html(gene_id, symbol, row, evidence_trace(gene_id, run), rank_column))

    missing_note = ""
    if run["missing"]:
        missing_note = (f'<p class="warn">Partial run: {escape(", ".join(run["missing"]))} not in the run directory — '
                        "the affected columns / trace sections are shown as unavailable below.</p>")
    trace_note = f"Evidence traces for all {len(details)} candidates."
    if len(details) < len(candidates):
        trace_note = (f"Evidence traces for the top {len(details)} of {len(candidates)} candidates "
                      f"(by <code>{escape(rank_column)}</code>; <code>--n-traces</code> to change).")
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    table_html = "\n".join(table_rows)
    details_html = "\n".join(details)

    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Targets from Pathways &mdash; {escape(disease['name'])}</title>
<style>{CSS}</style>
</head><body>
<h1>Targets from Pathways &mdash; candidate report</h1>
<p class="meta">
Disease <b>{escape(disease['name'])}</b> <span class=ens>{escape(disease['efo_id'])}</span> &middot;
{len(run["seeds"])} seed genes (genetic association) &middot; {len(candidates)} candidates &middot;
random seed {escape(manifest.get('seed'))} &middot; {escape(", ".join(versions))} &middot;
graph: {escape(graph_kind)} &middot; generated {generated}
</p>
{missing_note}
<p class="lead">
Candidates are ranked by <code>{escape(rank_column)}</code>, Stage 5's network propagation
from the seed genes; tractability and safety annotate them and do not
change the order. Click a header to sort; click a gene symbol to open its evidence trace.
</p>

<table id="ranktable">
<thead><tr>{header_html}</tr></thead>
<tbody id="ranktbody">
{table_html}
</tbody>
</table>

<h2>Evidence trace</h2>
<p class="lead">{trace_note}</p>
{details_html}

<script>{JAVASCRIPT}</script>
</body></html>
"""


def build_report(run_dir, out_path, n_traces=DEFAULT_N_TRACES):
    run = load_run(run_dir)
    out_path = Path(out_path)
    out_path.write_text(render_html(run, n_traces))
    logger.info("Report written (%d candidates): %s", len(run["annotated"]), out_path)
    return out_path


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="stage7_report.py",
        description="Stage 7 - output/report. Turns a finished run directory (manifest.json + "
        "candidates_annotated.tsv required; other files optional) into ONE self-contained interactive HTML "
        "file: sortable ranked table and click-to-expand evidence traces for the top candidates.",
    )
    parser.add_argument("--run-dir", required=True, type=Path, help="Directory with the Stage 1-6 outputs.")
    parser.add_argument("--out", type=Path, default=None, help="Output HTML path. Default: <run-dir>/report.html.")
    parser.add_argument("--n-traces", type=int, default=DEFAULT_N_TRACES,
                        help=f"Render evidence traces for this many top-ranked candidates. Default: {DEFAULT_N_TRACES}.")
    return parser.parse_args(argv)


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    args = parse_args(argv)
    if args.n_traces < 0:
        fail(f"--n-traces must be >= 0 (got {args.n_traces}).")
    build_report(args.run_dir, args.out or args.run_dir / "report.html", args.n_traces)


if __name__ == "__main__":
    main()
