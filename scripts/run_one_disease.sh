#!/bin/bash
# One disease end to end: finish the data download, run the tests, run the pipeline, then the three leave-one-out
# modes. Every step runs alone, so nothing competes for memory, and stops if the machine runs short of memory.
#
#   scripts/run_one_disease.sh check      # read-only: inputs, disk, memory, other jobs
#   scripts/run_one_disease.sh all        # every step below, in order, stopping at the first failure
#   scripts/run_one_disease.sh <step>...  # download tests pipeline loo-standard loo-known-drug loo-refit
#
# Settings, as environment variables, with absolute paths (defaults in brackets):
#   DATA             source data [<repo>/data]
#   GSEA_CATALOGUE   DisPathCatalogue v2 directory, for loo-refit [required by that step]
#   CATALOGUE        disease graph catalogue [$DATA/disease_graph_catalogue]
#   ASSOCIATIONS     association_by_datatype_indirect [$DATA/association_by_datatype_indirect]
#   FI_FILE          [$DATA/reactome/FIsInGene_04142025_with_annotations.txt]
#   PIPELINE_DISEASE [EFO_0003767]    LOO_DISEASE [EFO_0000095]
#   OUT              results [$DATA/runs/one_disease]
#   PY               [<repo>/.venv/bin/python]
#   START_FREE_PCT   do not start a step below this % of free memory [45]
#   STOP_FREE_PCT    stop a running step below this % of free memory [15]
#   MIN_FREE_GB      do not start a step below this much free disk [3]

set -eu -o pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
DATA=${DATA:-$REPO/data}
GSEA_CATALOGUE=${GSEA_CATALOGUE:-}
CATALOGUE=${CATALOGUE:-$DATA/disease_graph_catalogue}
ASSOCIATIONS=${ASSOCIATIONS:-$DATA/association_by_datatype_indirect}
FI_FILE=${FI_FILE:-$DATA/reactome/FIsInGene_04142025_with_annotations.txt}
PIPELINE_DISEASE=${PIPELINE_DISEASE:-EFO_0003767}
LOO_DISEASE=${LOO_DISEASE:-EFO_0000095}
OUT=${OUT:-$DATA/runs/one_disease}
PY=${PY:-$REPO/.venv/bin/python}
START_FREE_PCT=${START_FREE_PCT:-45}
STOP_FREE_PCT=${STOP_FREE_PCT:-15}
MIN_FREE_GB=${MIN_FREE_GB:-3}

INDIRECT_DIR=$DATA/opentargets/26.06/association_by_datasource_indirect
INDIRECT_PATH=pub/databases/opentargets/platform/26.06/output/association_by_datasource_indirect/
TARGET_DIR=$DATA/opentargets/26.06/target
# One part file of the indirect table is at most about 100 MB.
PART_MB=100
# Written once every part is downloaded and read in full.
INDIRECT_COMPLETE=$INDIRECT_DIR/.download_complete
LOG_DIR=$OUT/logs

export PYTHONUNBUFFERED=1

die() { echo "ERROR: $*" >&2; exit 1; }

# Empty where memory_pressure does not exist (not macOS): the memory checks are then skipped.
free_memory_pct() { memory_pressure -Q 2>/dev/null | awk -F': ' '/free percentage/ {print $2 + 0}' || true; }
free_disk_gb() { df -Pk "$1" | awk 'NR == 2 {print int($4 / 1048576)}'; }
other_jobs() { pgrep -fl "leave_one_out\.py|run_pipeline\.py|stage[1-7]_[a-z_]*\.py|-m pytest" || true; }
n_indirect_parts() { find "$INDIRECT_DIR" -maxdepth 1 -name '*.snappy.parquet' 2>/dev/null | wc -l | tr -d ' '; }

check() {
    local problems=0 path free
    for path in "$PY" "$CATALOGUE/${LOO_DISEASE}_disease_graph" "$ASSOCIATIONS" "$FI_FILE" "$TARGET_DIR" \
                "$DATA/opentargets/26.06/disease" "$DATA/reactome/ReactomePathways.gmt" \
                "$GSEA_CATALOGUE/merged_gene_sets.gmt" "$GSEA_CATALOGUE/genetic_gene_scores_input.parquet"; do
        if [ -e "$path" ]; then echo "ok       $path"; else echo "MISSING  $path"; problems=1; fi
    done
    if [ -e "$INDIRECT_COMPLETE" ] || [ -e "$INDIRECT_DIR/.checksums.json" ]; then
        echo "ok       indirect table complete"
    else
        echo "todo     indirect table: $(n_indirect_parts) parts here, the download step finishes and validates it"
    fi
    echo "disk     $(free_disk_gb "$DATA") GB free (steps need at least $MIN_FREE_GB; the download needs about $PART_MB MB per missing part on top)"
    free=$(free_memory_pct)
    echo "memory   ${free:-unknown}% free (steps need at least $START_FREE_PCT%)"
    if [ -n "$free" ] && [ "$free" -lt "$START_FREE_PCT" ]; then echo "         close other apps first"; problems=1; fi
    if [ -n "$(other_jobs)" ]; then echo "RUNNING  other jobs:"; other_jobs; problems=1; fi
    if command -v pmset >/dev/null && ! pmset -g batt | grep -q "AC Power"; then
        echo "note     on battery: plug in, and keep the lid open (closing it sleeps the machine)"
    fi
    [ "$problems" -eq 0 ] || die "fix the items above before running."
    echo "Ready."
}

# Runs a command in its own process group, with its output in $LOG_DIR/<name>.log and on screen. Stops the whole
# group, worker processes included, if free memory falls below STOP_FREE_PCT.
PID=""
stop_group() { [ -z "$PID" ] || { kill -TERM -- "-$PID" 2>/dev/null || true; sleep 5; kill -KILL -- "-$PID" 2>/dev/null || true; }; }
trap 'stop_group; exit 130' INT TERM

run_step() {
    local name=$1 log free tail_pid rc=0 low_memory=0
    shift
    log=$LOG_DIR/$name.log
    mkdir -p "$LOG_DIR"
    [ -z "$(other_jobs)" ] || { other_jobs; die "$name: another job is running."; }
    [ "$(free_disk_gb "$DATA")" -ge "$MIN_FREE_GB" ] || die "$name: less than $MIN_FREE_GB GB of free disk."
    free=$(free_memory_pct)
    if [ -n "$free" ] && [ "$free" -lt "$START_FREE_PCT" ]; then
        die "$name: only $free% of memory is free (needs $START_FREE_PCT%). Close other apps and run this step again."
    fi
    echo "=== $name ($(date '+%H:%M:%S')): $*"
    echo "    log: $log"
    : > "$log"
    set -m
    (cd "$REPO" && exec "$@") >> "$log" 2>&1 &
    PID=$!
    set +m
    tail -f "$log" &
    tail_pid=$!
    while kill -0 "$PID" 2>/dev/null; do
        free=$(free_memory_pct)
        if [ -n "$free" ] && [ "$free" -lt "$STOP_FREE_PCT" ]; then
            low_memory=1
            echo "STOPPING $name: free memory is down to $free%." | tee -a "$log"
            stop_group
            break
        fi
        sleep 5
    done
    wait "$PID" || rc=$?
    PID=""
    kill "$tail_pid" 2>/dev/null || true
    wait "$tail_pid" 2>/dev/null || true
    [ "$low_memory" -eq 0 ] || die "$name stopped to protect the machine. See $log."
    [ "$rc" -eq 0 ] || die "$name failed with exit code $rc. See $log."
    echo "=== $name done ($(date '+%H:%M:%S'))"
}

download() {
    if [ -e "$INDIRECT_COMPLETE" ] || [ -e "$INDIRECT_DIR/.checksums.json" ]; then
        echo "=== download: already complete."
        return 0
    fi
    local url="https://ftp.ebi.ac.uk/$INDIRECT_PATH" listing names name n_missing=0
    if ! listing=$(curl -sS --fail --max-time 60 "$url"); then
        url="http://ftp.ebi.ac.uk/$INDIRECT_PATH"
        echo "HTTPS to ftp.ebi.ac.uk failed: using plain HTTP. Each part is read in full below to validate it."
        listing=$(curl -sS --fail --max-time 60 "$url") || die "cannot list $url"
    fi
    names=$(printf '%s\n' "$listing" | grep -oE 'href="[^"?]*\.snappy\.parquet"' | sed 's/href="//;s/"$//' | sort -u)
    [ -n "$names" ] || die "no parquet files listed at $url"
    mkdir -p "$INDIRECT_DIR"
    rm -f "$INDIRECT_DIR"/*.part
    for name in $names; do [ -f "$INDIRECT_DIR/$name" ] || n_missing=$((n_missing + 1)); done
    echo "=== download: $n_missing of $(printf '%s\n' "$names" | wc -l | tr -d ' ') parts missing."
    if [ "$(free_disk_gb "$DATA")" -lt $((MIN_FREE_GB + n_missing * PART_MB / 1000 + 1)) ]; then
        die "not enough disk for $n_missing parts: free some space first."
    fi
    for name in $names; do
        [ ! -f "$INDIRECT_DIR/$name" ] || continue
        echo "    $name"
        curl -sS --fail --retry 3 --max-time 1800 -o "$INDIRECT_DIR/$name.part" "$url$name" \
            || { rm -f "$INDIRECT_DIR/$name.part"; die "download failed: $name"; }
        mv "$INDIRECT_DIR/$name.part" "$INDIRECT_DIR/$name"
    done
    for name in $names; do [ -f "$INDIRECT_DIR/$name" ] || die "still missing: $name"; done
    # Batch by batch, so validating does not hold a whole table in memory.
    run_step validate-download "$PY" -c '
import glob, sys
import pyarrow.parquet as pq
rows = 0
for path in sorted(glob.glob(sys.argv[1] + "/*.snappy.parquet")):
    for batch in pq.ParquetFile(path).iter_batches(batch_size=100_000):
        rows += batch.num_rows
print("every part read in full,", rows, "rows")' "$INDIRECT_DIR"
    touch "$INDIRECT_COMPLETE"
}

run_tests() { run_step tests "$PY" -m pytest -q; }

# The complete table is on disk at this point, so --no-download pins all of it. --allow-low-coverage: Stage 1 stopped
# on Reactome's Ensembl2Reactome mapping at 89.85% in the last run, just under its 90% floor.
pipeline() {
    [ -e "$INDIRECT_COMPLETE" ] || [ -e "$INDIRECT_DIR/.checksums.json" ] \
        || die "pipeline: run the download step first; a partial table would silently give partial evidence."
    run_step "pipeline-$PIPELINE_DISEASE" "$PY" pipeline/run_pipeline.py --disease "$PIPELINE_DISEASE" \
        --data-dir "$DATA" --out-dir "$OUT/pipeline_$PIPELINE_DISEASE" --no-download --allow-low-coverage \
        --method bfwalk,rwr --report
}

loo_standard() {
    run_step loo-standard "$PY" validation/leave_one_out.py --catalogue "$CATALOGUE" --out-dir "$OUT/loo_standard" \
        --diseases "$LOO_DISEASE" --jobs 1
}

loo_known_drug() {
    run_step loo-known-drug "$PY" validation/leave_one_out.py --catalogue "$CATALOGUE" \
        --out-dir "$OUT/loo_known_drug" --diseases "$LOO_DISEASE" --holdout known_drug --associations "$ASSOCIATIONS"
}

loo_refit() {
    [ -n "$GSEA_CATALOGUE" ] || die "loo-refit: set GSEA_CATALOGUE to the DisPathCatalogue v2 directory."
    run_step loo-refit "$PY" validation/leave_one_out.py --catalogue "$CATALOGUE" --out-dir "$OUT/loo_refit" \
        --diseases "$LOO_DISEASE" --jobs 1 --refit-gsea --gsea-catalogue-dir "$GSEA_CATALOGUE" \
        --target-dir "$TARGET_DIR" --fi-file "$FI_FILE"
}

[ "$#" -gt 0 ] || { sed -n 2,7p "$0" | sed 's/^# \{0,1\}//'; exit 2; }
if [ "$1" = all ]; then set -- check download tests pipeline loo-standard loo-known-drug loo-refit; fi

# Keeps the machine awake for as long as this script runs (not with the lid closed).
if command -v caffeinate >/dev/null; then caffeinate -i -w $$ & fi

for step in "$@"; do
    case $step in
        check) check ;;
        download) download ;;
        tests) run_tests ;;
        pipeline) pipeline ;;
        loo-standard) loo_standard ;;
        loo-known-drug) loo_known_drug ;;
        loo-refit) loo_refit ;;
        *) die "unknown step: $step" ;;
    esac
done
echo "Finished: $*. Results under $OUT."
