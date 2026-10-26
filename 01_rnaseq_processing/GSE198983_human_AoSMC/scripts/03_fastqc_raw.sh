#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────
# STEP 03 — 원본 FASTQ QC (FastQC + MultiQC)
#   품질/어댑터/중복/GC 분포를 트리밍 전에 확인합니다.
# ─────────────────────────────────────────────────────────────
source "$(dirname "$0")/00_config.sh"

OUT="$RESULT_DIR/03_fastqc_raw"
mkdir -p "$OUT"

FQ=("$FASTQ_DIR"/*.fastq.gz)
log "FastQC 대상 ${#FQ[@]}개 파일"

TODO=()
for f in "${FQ[@]}"; do
    b=$(basename "$f" .fastq.gz)
    [[ -s "$OUT/${b}_fastqc.html" ]] || TODO+=("$f")
done

if [[ ${#TODO[@]} -gt 0 ]]; then
    log "실행 대상 ${#TODO[@]}개"
    "$FASTQC" -t "$THREADS" -o "$OUT" "${TODO[@]}" 2>&1 | tee -a "$LOG_DIR/03_fastqc_raw.log"
else
    log "이미 전부 완료 — 건너뜀"
fi

log "MultiQC 통합"
"$MULTIQC" "$OUT" -o "$OUT" -n multiqc_raw --force 2>&1 | tail -3 | tee -a "$LOG_DIR/03_fastqc_raw.log"

log "완료 → $OUT/multiqc_raw.html"
