#!/usr/bin/env bash
# [deposit] created for the code deposit: end-to-end smoke test of the raw-data steps on the tiny test subset
# (1 Mb chr1 window around LMNA, ~6,000 read pairs per sample; see make_test_subset.sh).
#
# Runs, with PROJ pointing at a WORK copy OUTSIDE the deposit (nothing is written next to this script):
#   01 (download + md5 logic, against one tiny public ENA run; needs network; DOWNLOAD_TEST=0 skips it)
#   02 (pre-placed mini genome/GTF/transcripts .gz -> downloads skipped; unzip + STAR + salmon indices)
#   02b (GTF -> BED12) -> 03 -> 04 -> 05 -> 06 (+06b dupRadar, Qualimap) -> 07 -> 08 -> 13   ([deposit] 15_sashimi.sh removed: figure-only)
#   and finally run_pipeline.sh with ONLY=07 to check the runner wiring.
#
# Usage:  bash test/run_test.sh                     (WORK defaults to a new dir under ${TMPDIR:-/tmp})
#         WORK=/path/outside/deposit THREADS=4 bash test/run_test.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"          # .../GSE198983_human_AoSMC/test
DS="$(dirname "$HERE")"                                        # dataset dir (holds scripts/)
SCRIPTS="$DS/scripts"
WORK="${WORK:-$(mktemp -d "${TMPDIR:-/tmp}/aosmc_test.XXXXXX")}"
mkdir -p "$WORK"; WORK="$(cd "$WORK" && pwd)"
case "$WORK/" in "$DS/"*) echo "ERROR: WORK must be outside the dataset dir ($DS)" >&2; exit 1 ;; esac

# test-mode settings (all are ordinary env overrides understood by the scripts)
export PROJ="$WORK/proj"
export THREADS="${THREADS:-4}"
export STAR_GENOME_EXTRA="${STAR_GENOME_EXTRA:---genomeSAindexNbases 8}"   # 1 Mb genome: min(14, log2(1e6)/2-1)≈9
export PICARD_MEM="${PICARD_MEM:--Xmx4g}"
export SKIP_QUALIMAP="${SKIP_QUALIMAP:-0}"                                 # exercise the Qualimap branch too
DOWNLOAD_TEST="${DOWNLOAD_TEST:-1}"

echo "test PROJ = $PROJ"
mkdir -p "$PROJ"
cp -r "$HERE/ref" "$HERE/raw_data" "$HERE/meta" "$PROJ/"

step() { local t0; t0=$(date +%s); echo; echo "════ $1 ════"; shift; "$@"; echo "──── ok ($(( $(date +%s) - t0 )) s)"; }

# ── 01: download logic against a tiny public run (GET + md5, then SKIP on rerun) ──
if [[ "$DOWNLOAD_TEST" == "1" ]]; then
    DL="$WORK/dl_test"; mkdir -p "$DL/meta"
    cp "$HERE/meta_download/ena_filereport.tsv" "$DL/meta/"
    step "01 download (GET)"  env PROJ="$DL" PAR=2 bash "$SCRIPTS/01_download_fastq.sh"
    step "01 download (SKIP)" env PROJ="$DL" PAR=2 bash "$SCRIPTS/01_download_fastq.sh"
fi

# ── 02..13 on the mini genome ──
step "02 reference + indices" bash "$SCRIPTS/02_reference.sh"
step "02b BED12" bash -c 'set -a; source "$0/00_config.sh"; set +a
                          "$RSCRIPT" "$0/02b_make_bed12.R" "$GTF" "$BED12" "$BED12_HK"' "$SCRIPTS"
step "03 FastQC raw"    bash "$SCRIPTS/03_fastqc_raw.sh"
step "04 fastp trim"    bash "$SCRIPTS/04_trim.sh"
step "05 STAR align"    bash "$SCRIPTS/05_star_align.sh"
step "06 post-align QC" bash "$SCRIPTS/06_post_align_qc.sh"
step "07 featureCounts" bash "$SCRIPTS/07_featurecounts.sh"
step "08 salmon"        bash "$SCRIPTS/08_salmon.sh"
step "13 rMATS"         bash "$SCRIPTS/13_rmats.sh"

# ── runner wiring: re-run step 07 through run_pipeline.sh (scripts dir = runner's own dir) ──
rm -f "$PROJ/results/07_counts/counts_matrix.tsv"
step "run_pipeline.sh ONLY=07" env ONLY=07 bash "$SCRIPTS/run_pipeline.sh"
[[ -s "$PROJ/results/07_counts/counts_matrix.tsv" ]]

echo
echo "════ summary ════"
R="$PROJ/results"
echo "BAMs:        $(ls "$R"/05_star/*.Aligned.sortedByCoord.out.bam | wc -l)"
echo "counts:      $(( $(wc -l < "$R/07_counts/counts_matrix.tsv") - 1 )) genes"
echo "salmon:      $(ls "$R"/08_salmon/*/quant.sf | wc -l) quant.sf"
echo "rMATS:";  column -t "$R/13_rmats/event_summary.tsv"
echo "ALL TEST STEPS PASSED  (outputs in $WORK)"
