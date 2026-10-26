#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────
# [deposit] created for the code deposit — runs the raw-data part of the
# GSE165409 pipeline (scripts/01 … 07, 13) on the tiny test subset in this
# folder, to show that the scripts run end-to-end without the ~72 GB FASTQ and
# ~36 GB reference. Numbers produced here are NOT meaningful biologically.
#
# Usage:  bash run_test.sh WORK_DIR            (or: PROJ=WORK_DIR bash run_test.sh)
#   WORK_DIR must be outside the deposit; the subset is copied there and used as PROJ.
#   env: THREADS (default 4), CONDA_ENVS (default ~/miniforge3/envs),
#        RUN_DOWNSTREAM=1  additionally run 08_deseq2.R and 09_rmats_figures.R on the test counts
#
# Steps:  build_test_reference.sh → 01 FastQC → 02 fastp → 03 STAR → 04 RSeQC → 05 featureCounts
#         (+ author-count comparison) → 06 salmon → 07 rMATS → 13 dupRadar
#         (10_sashimi.sh is figure-only and not part of the deposit)
# Logs:   WORK_DIR/logs/test_<step>.log, summary in WORK_DIR/logs/TEST_STATUS
# ─────────────────────────────────────────────────────────────
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATASET="$(cd "$HERE/.." && pwd)"
SCRIPTS="$DATASET/scripts"
DEPOSIT="$(cd "$DATASET/../.." && pwd)"

WORK="${1:-${PROJ:-}}"
[[ -n "$WORK" ]] || { echo "usage: bash run_test.sh WORK_DIR   (WORK_DIR outside the deposit)" >&2; exit 2; }
mkdir -p "$WORK"; WORK="$(cd "$WORK" && pwd -P)"
case "$WORK/" in
    "$(cd "$DEPOSIT" && pwd -P)/"*) echo "ERROR: WORK_DIR ($WORK) is inside the deposit; choose a directory outside it" >&2; exit 2 ;;
esac

export PROJ="$WORK"
export THREADS="${THREADS:-4}"
export SORT_MEM="${SORT_MEM:-1G}"
export CONDA_ENVS="${CONDA_ENVS:-$HOME/miniforge3/envs}"
export SKIP_QUALIMAP="${SKIP_QUALIMAP:-1}"
export TINY_TEST=1          # lets 08_deseq2.R use varianceStabilizingTransformation (<1000 genes)

# ── stage inputs into WORK (the deposit itself is only read) ──
mkdir -p "$WORK"/{raw_data/fastq,raw_data/geo,meta,ref,logs}
cp "$HERE"/raw_data/fastq/*.fastq.gz "$WORK/raw_data/fastq/"
cp "$HERE/meta/samples.tsv" "$WORK/meta/samples.tsv"
cp "$DATASET"/raw_data/geo/* "$WORK/raw_data/geo/"      # author counts for the 05 comparison
cd "$WORK"                                               # stray files (e.g. Rplots.pdf) land here

STATUS="$WORK/logs/TEST_STATUS"
step() {
    local name="$1"; shift
    local t0; t0=$(date +%s)
    if "$@" > "$WORK/logs/test_${name}.log" 2>&1; then
        echo "PASS  $name  $(( $(date +%s) - t0 ))s" | tee -a "$STATUS"
    else
        echo "FAIL  $name  $(( $(date +%s) - t0 ))s  (see logs/test_${name}.log)" | tee -a "$STATUS"
        tail -20 "$WORK/logs/test_${name}.log" >&2
        exit 1
    fi
}

: > "$STATUS"
step 00_reference  bash "$HERE/build_test_reference.sh"
step 01_fastqc     bash "$SCRIPTS/01_fastqc_raw.sh"
step 02_trim       bash "$SCRIPTS/02_trim.sh"
step 03_star       bash "$SCRIPTS/03_star_align.sh"
step 04_qc_align   bash "$SCRIPTS/04_post_align_qc.sh"
step 05_counts     bash "$SCRIPTS/05_featurecounts.sh"
step 06_salmon     bash "$SCRIPTS/06_salmon.sh"
step 07_rmats      bash "$SCRIPTS/07_rmats.sh"
step 13_dupradar   bash "$SCRIPTS/13_dupradar.sh"
if [[ "${RUN_DOWNSTREAM:-0}" == "1" ]]; then
    step 08_deseq2        "$CONDA_ENVS/lmna_r/bin/Rscript" "$SCRIPTS/08_deseq2.R"
    step 09_rmats_figures "$CONDA_ENVS/lmna_r/bin/Rscript" "$SCRIPTS/09_rmats_figures.R"
fi
echo "test finished -> $WORK (status: $STATUS)"
