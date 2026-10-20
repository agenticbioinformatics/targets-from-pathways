"""Stage 7 - HTML report.

Turns a run directory into one self-contained interactive HTML file (no server, no external assets):
a sortable ranked table, an evidence trace per candidate and a static benchmark panel.
Only manifest.json and candidates_annotated.tsv are required; a missing optional file only hides the
report sections that need it. The benchmark panel embeds benchmarking/benchmark_summary.txt
(or the .example.txt placeholder) and is never recomputed here.
"""

import argparse
import html
import json
import logging
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from stage0_schemas import AnnotatedCandidatesSchema, fail

PIPELINE_DIR = Path(__file__).resolve().parent
REPO_DIR = PIPELINE_DIR.parent

logger = logging.getLogger("stage7_report")

SIGN_LABELS = {1: "activating (+)", -1: "inhibiting (−)", 0: "unsigned (0)"}
# Stage 4's default --datatypes, used when graph_metadata.json does not record them.
DEFAULT_ADMITTED_DATATYPES = ["genetic_association", "known_drug"]

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
                 "interactions.parquet", "ot_disease_subset.parquet", "disease_pathways.tsv"]:
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
    stage4_info = graph_metadata.get("stage4_genetic_evidence", {})
    graph_index = None
    if optional["graph_gene_index.json"] is not None:
        graph_index = set(optional["graph_gene_index.json"])
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
        "admitted_datatypes": set(stage4_info.get("datatypes", DEFAULT_ADMITTED_DATATYPES)),
        "gene_sets": optional["gene_sets.parquet"],
        "interactions": optional["interactions.parquet"],
        "ot": optional["ot_disease_subset.parquet"],
        "disease_pathways": optional["disease_pathways.tsv"],
        "stage4_ran": "stage4_genetic_evidence" in graph_metadata,
        "missing": missing,
    }


def symbol_of(gene_id, symbols):
    return symbols.get(gene_id, gene_id)


def evidence_trace(gene_id, run):
    """Collect why a candidate is near the target: shared pathways, Open Targets datatypes and interactions.
    Stage 5 keeps only scores, so this is rebuilt here from the Stage 1-3 files."""
    target_id = run["manifest"]["resolved_target"]["gene_id"]
    symbols = run["symbols"]
    gene_sets = run["gene_sets"]
    ot = run["ot"]
    interactions_df = run["interactions"]
    graph_index = run["graph_index"]

    shared_pathways = []
    member_only_pathways = []
    if gene_sets is not None:
        set_names = dict(zip(gene_sets["set_id"], gene_sets["set_name"]))
        target_sets = set(gene_sets.loc[gene_sets["gene_id"] == target_id, "set_id"])
        disease_sets = set()
        if run["disease_pathways"] is not None:
            disease_sets = set(run["disease_pathways"]["set_id"])
        # Only pathways that Stage 3 used to build the graph.
        candidate_sets = set(gene_sets.loc[gene_sets["gene_id"] == gene_id, "set_id"]) & (disease_sets | target_sets)
        for set_id in candidate_sets:
            if set_id in target_sets:
                shared_pathways.append(set_names.get(set_id, set_id))
            else:
                member_only_pathways.append(set_names.get(set_id, set_id))
        shared_pathways.sort(key=str.lower)
        member_only_pathways.sort(key=str.lower)

    datatypes = []
    if ot is not None:
        for datatype, rows in ot[ot["gene_id"] == gene_id].groupby("datatype_id"):
            datatypes.append({
                "datatype": str(datatype),
                "score": float(rows["score"].max()),
                "datasources": sorted(rows["datasource_id"].astype(str).unique()),
                # "Feeds" the genetic-evidence score only if Stage 4 ran with this datatype.
                "feeds": run["stage4_ran"] and str(datatype) in run["admitted_datatypes"],
            })
        # Datatypes that feed the score first, then by score.
        datatypes.sort(key=lambda entry: (not entry["feeds"], -entry["score"]))

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
                "to_target": other_gene == target_id,
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


def benchmark_summary(benchmark_dir=None):
    """Return (summary text or None, label of its source)."""
    if benchmark_dir is None:
        benchmark_dir = REPO_DIR / "benchmarking"
    benchmark_dir = Path(benchmark_dir)
    real_summary = benchmark_dir / "benchmark_summary.txt"
    example_summary = benchmark_dir / "benchmark_summary.example.txt"
    if real_summary.exists():
        return real_summary.read_text().strip(), "benchmark_summary.txt"
    if example_summary.exists():
        return example_summary.read_text().strip(), "benchmark_summary.example.txt (EXAMPLE — not a real validation run)"
    return None, "none"


def escape(value):
    return html.escape(str(value))


def html_list(items):
    text = "<ul>"
    for item in items:
        text += f"<li>{item}</li>"
    return text + "</ul>"


def trace_html(gene_id, symbol, row, trace):
    parts = [
        f'<details id="g-{escape(gene_id)}">',
        f"<summary><b>{escape(symbol)}</b> <span class=ens>{escape(gene_id)}</span> "
        f"&mdash; composite <b>{row.composite_score:.3f}</b></summary>",
        '<div class="trace">',
        '<p class="k">Composite breakdown</p>',
        f"<p><code>{escape(row.composite_breakdown)}</code><br>weights <code>{escape(row.composite_weights)}</code></p>",
        '<p class="k">Shared pathways with the target</p>',
    ]
    if not trace["pathways_available"]:
        parts.append('<p class="none">gene_sets.parquet not in the run directory</p>')
    elif trace["shared_pathways"]:
        names = []
        for name in trace["shared_pathways"]:
            names.append(escape(name))
        parts.append(html_list(names))
    else:
        parts.append('<p class="none">none in the Stage 3 pathway union</p>')

    if trace["member_pathways"]:
        parts.append('<p class="k">Also a member of (disease-relevant; target is not)</p>')
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
            badge = ""
            if datatype["feeds"]:
                badge = '<span class="badge">feeds genetic-evidence score</span>'
            sources = escape(", ".join(datatype["datasources"]))
            items.append(f"{escape(datatype['datatype'])} &mdash; {datatype['score']:.3f}{badge} "
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
            to_target = ""
            if interaction["to_target"]:
                to_target = '<span class="badge tgt">to target</span>'
            items.append(f"{escape(interaction['a'])} {arrow} {escape(interaction['b'])} &mdash; "
                         f"{escape(interaction['sign'])}, {escape(interaction['evidence_type'])}{confidence}{to_target}")
        parts.append(html_list(items))
    else:
        parts.append('<p class="none">no signed interaction edges touch this gene</p>')

    parts.append("</div></details>")
    return "".join(parts)


def render_html(run, benchmark_text, benchmark_label):
    manifest = run["manifest"]
    candidates = run["annotated"].reset_index(drop=True)
    symbols = run["symbols"]
    target = manifest["resolved_target"]
    disease = manifest["disease"]
    target_symbol = escape(symbol_of(target["gene_id"], symbols))

    versions = []
    for source in manifest.get("sources", []):
        versions.append(f"{source['db']} {source['version']}")
    if run["stage4_ran"]:
        graph_kind = "Stage 4 (genetic-evidence-weighted)"
    elif run["graph_meta"]:
        graph_kind = "Stage 3 (structural)"
    else:
        graph_kind = "unknown (no graph_metadata.json)"

    has_rwr = "rwr_score" in candidates.columns
    has_bfwalk = "bfwalk_score" in candidates.columns
    has_genetic_evidence = "genetic_evidence_score" in candidates.columns

    # (column title, sort numerically, left-aligned text)
    headers = [("#", True, False), ("symbol", False, True), ("gene", False, True),
               ("composite", True, False), ("topology", True, False)]
    if has_rwr:
        headers.append(("rwr", True, False))
    if has_bfwalk:
        headers.append(("bfwalk", True, False))
    if has_genetic_evidence:
        headers.append(("genetic ev.", True, False))
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
        cells = [
            f'<td data-v="{rank}">{rank}</td>',
            f'<td class="txt"><a href="#" onclick="openGene(\'g-{escape(gene_id)}\');return false;">{escape(symbol)}</a></td>',
            f'<td class="txt ens">{escape(gene_id)}</td>',
            f'<td data-v="{row.composite_score}">{row.composite_score:.3f}</td>',
            f'<td data-v="{row.topology_score}">{row.topology_score:.3f}</td>',
        ]
        if has_rwr:
            cells.append(f'<td data-v="{row.rwr_score}">{row.rwr_score:.4f}</td>')
        if has_bfwalk:
            cells.append(f'<td data-v="{row.bfwalk_score}">{row.bfwalk_score:.4f}</td>')
        if has_genetic_evidence:
            cells.append(f'<td data-v="{row.genetic_evidence_score}">{row.genetic_evidence_score:.3f}</td>')
        cells.append(f'<td class="txt">{escape(row.tractability)}</td>')
        cells.append(f'<td class="txt">{escape(row.safety)}</td>')
        table_rows.append("<tr>" + "".join(cells) + "</tr>")
        details.append(trace_html(gene_id, symbol, row, evidence_trace(gene_id, run)))

    if benchmark_text is None:
        benchmark_html = ('<p class="none">No benchmark summary found. Run '
                          "<code>benchmarking/benchmark_validate.py</code> and re-generate this report.</p>")
    else:
        css_class = "meta"
        if "EXAMPLE" in benchmark_label:
            css_class = "warn"
        benchmark_html = f'<p class="{css_class}">Source: {escape(benchmark_label)}</p><pre>{escape(benchmark_text)}</pre>'

    missing_note = ""
    if run["missing"]:
        missing_note = (f'<p class="warn">Partial run: {escape(", ".join(run["missing"]))} not in the run directory — '
                        "the affected columns / trace sections are shown as unavailable below.</p>")
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    table_html = "\n".join(table_rows)
    details_html = "\n".join(details)

    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Target Rescue via Pathways &mdash; {target_symbol} / {escape(disease['name'])}</title>
<style>{CSS}</style>
</head><body>
<h1>Target Rescue via Pathways &mdash; candidate report</h1>
<p class="meta">
Target <b>{target_symbol}</b> <span class=ens>{escape(target['gene_id'])}</span>
(input <code>{escape(target['input'])}</code>) &middot;
Disease <b>{escape(disease['name'])}</b> <span class=ens>{escape(disease['efo_id'])}</span> &middot;
{len(candidates)} candidates &middot; seed {escape(manifest.get('seed'))} &middot; {escape(", ".join(versions))} &middot;
graph: {escape(graph_kind)} &middot; generated {generated}
</p>
{missing_note}
<p class="lead">
<code>composite_score</code> is a weighted average of the normalised components; each
candidate's <i>breakdown</i> below shows every component's weighted contribution
(the terms sum to the score). Click a header to sort; click a gene symbol to open
its evidence trace.
</p>

<table id="ranktable">
<thead><tr>{header_html}</tr></thead>
<tbody id="ranktbody">
{table_html}
</tbody>
</table>

<h2>Evidence trace</h2>
{details_html}

<h2>Benchmark validation</h2>
<p class="meta">Static &mdash; embedded from <code>benchmarking/</code>, not recomputed for this report.</p>
{benchmark_html}

<script>{JAVASCRIPT}</script>
</body></html>
"""


def build_report(run_dir, out_path, benchmark_dir=None):
    run = load_run(run_dir)
    benchmark_text, benchmark_label = benchmark_summary(benchmark_dir)
    out_path = Path(out_path)
    out_path.write_text(render_html(run, benchmark_text, benchmark_label))
    logger.info("Report written (%d candidates): %s", len(run["annotated"]), out_path)
    return out_path


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="stage7_report.py",
        description="Stage 7 - output/report. Turns a finished run directory (manifest.json + "
        "candidates_annotated.tsv required; other files optional) into ONE self-contained interactive HTML "
        "file: sortable ranked table, click-to-expand evidence trace per candidate, and a static "
        "benchmark-validation panel.",
    )
    parser.add_argument("--run-dir", required=True, type=Path, help="Directory with the Stage 1-6 outputs.")
    parser.add_argument("--out", type=Path, default=None, help="Output HTML path. Default: <run-dir>/report.html.")
    parser.add_argument("--benchmark-dir", type=Path, default=None,
                        help="Directory holding benchmark_summary.txt / .example.txt. Default: <repo>/benchmarking.")
    return parser.parse_args(argv)


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    args = parse_args(argv)
    build_report(args.run_dir, args.out or args.run_dir / "report.html", benchmark_dir=args.benchmark_dir)


if __name__ == "__main__":
    main()
