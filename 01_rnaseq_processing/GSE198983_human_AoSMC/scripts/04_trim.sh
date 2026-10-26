#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────
# STEP 04 — fastp 트리밍 (어댑터 자동 검출 + 품질 필터) + 트림 후 FastQC
#
# 옵션 근거:
#   --detect_adapter_for_pe  PE 어댑터 자동 검출 (TruSeq 등 하드코딩 불필요)
#   -q 15 -u 40              Q15 미만 염기가 40%를 넘는 read 제거 (fastp 기본)
#   -l 36                    36bp 미만은 정렬 신뢰도가 낮아 제거
#   --trim_poly_g            NovaSeq(2색 화학)의 polyG artifact 제거 ★필수
#   --correct                PE overlap 영역 염기 보정
# ─────────────────────────────────────────────────────────────
source "$(dirname "$0")/00_config.sh"

OUT="$RESULT_DIR/04_trimmed"
REP="$RESULT_DIR/04_fastp_reports"
QC="$RESULT_DIR/04_fastqc_trimmed"
mkdir -p "$OUT" "$REP" "$QC"

sample_rows | while IFS=$'\t' read -r sample gsm srr group replicate; do
    O1="$OUT/${srr}_1.trim.fastq.gz"
    O2="$OUT/${srr}_2.trim.fastq.gz"
    if [[ -s "$O1" && -s "$O2" ]]; then
        log "  $sample ($srr) — 이미 완료, 건너뜀"
        continue
    fi

    I1="$FASTQ_DIR/${srr}_1.fastq.gz"
    I2="$FASTQ_DIR/${srr}_2.fastq.gz"
    require_file "$I1"; require_file "$I2"

    log "  $sample ($group, $srr) 트리밍 중"
    "$FASTP" \
        -i "$I1" -I "$I2" \
        -o "$O1" -O "$O2" \
        --detect_adapter_for_pe \
        -q 15 -u 40 -l 36 \
        --trim_poly_g \
        --correction \
        -w 16 \
        -j "$REP/${sample}.fastp.json" \
        -h "$REP/${sample}.fastp.html" \
        2>&1 | tail -20 | tee -a "$LOG_DIR/04_trim.log"
done

# ── 트림 후 FastQC ────────────────────────────────────────────
TRIMMED=("$OUT"/*.trim.fastq.gz)
TODO=()
for f in "${TRIMMED[@]}"; do
    b=$(basename "$f" .fastq.gz)
    [[ -s "$QC/${b}_fastqc.html" ]] || TODO+=("$f")
done
if [[ ${#TODO[@]} -gt 0 ]]; then
    log "트림 후 FastQC ${#TODO[@]}개"
    "$FASTQC" -t "$THREADS" -o "$QC" "${TODO[@]}" 2>&1 | tail -3 | tee -a "$LOG_DIR/04_trim.log"
fi

log "MultiQC 통합 (fastp + FastQC)"
"$MULTIQC" "$REP" "$QC" -o "$REP" -n multiqc_trimmed --force 2>&1 | tail -3 | tee -a "$LOG_DIR/04_trim.log"

log "완료 → $REP/multiqc_trimmed.html"
