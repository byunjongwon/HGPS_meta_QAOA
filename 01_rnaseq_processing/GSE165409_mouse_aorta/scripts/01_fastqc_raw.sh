#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────
# STEP 1. Raw read QC — FastQC (72 파일) + MultiQC 통합
#
# 산출 그림:
#   - per-base quality boxplot
#   - per-sequence quality histogram
#   - per-base sequence content
#   - GC content 분포
#   - sequence duplication levels
#   - adapter content curve
#   - MultiQC 종합 히트맵 (36 run x 지표)
#
# 예상 시간: 20코어에서 30~40분
# ─────────────────────────────────────────────────────────────
source "$(dirname "$0")/00_config.sh"

OUT="$RESULT_DIR/01_fastqc_raw"
mkdir -p "$OUT"

log "FastQC 시작 — 대상 $(ls "$FASTQ_DIR"/*.fastq.gz | wc -l) 파일"

# 이미 처리한 파일은 건너뛴다 (중단 후 재실행 대비)
TODO=()
for f in "$FASTQ_DIR"/*.fastq.gz; do
    base=$(basename "$f" .fastq.gz)
    [[ -s "$OUT/${base}_fastqc.zip" ]] || TODO+=("$f")
done

if [[ ${#TODO[@]} -gt 0 ]]; then
    log "처리할 파일 ${#TODO[@]}개"
    "$FASTQC" -t "$THREADS" -o "$OUT" --noextract "${TODO[@]}" 2>&1 | tee -a "$LOG_DIR/01_fastqc_raw.log"
else
    log "모든 파일이 이미 처리됨 — FastQC 건너뜀"
fi

log "MultiQC 통합 리포트 생성"
"$MULTIQC" "$OUT" -o "$OUT" -n multiqc_raw --force 2>&1 | tee -a "$LOG_DIR/01_fastqc_raw.log"

log "완료 → $OUT/multiqc_raw.html"
