#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────
# STEP 2. 어댑터 트리밍 — fastp (run 36개, paired-end)
#
# 산출 그림:
#   - fastp HTML 리포트 (before/after 품질 곡선, 어댑터 제거량)
#   - MultiQC가 fastp JSON을 읽어 통합 요약
#
# ⚠️ rMATS 주의:
#   트리밍하면 read 길이가 가변이 됩니다. rMATS 실행 시 반드시
#   --variable-read-length 옵션을 붙여야 합니다 (07 스크립트에 반영됨).
#   대안으로 --length_required 150 을 주면 길이가 고정되지만
#   어댑터가 낀 read를 통째로 버리게 되어 권장하지 않습니다.
#
# 예상 시간: 36 run x 약 2분 = 1시간 내외
# 디스크: 트리밍 결과 약 65~70 G 추가
# ─────────────────────────────────────────────────────────────
source "$(dirname "$0")/00_config.sh"

OUT="$RESULT_DIR/02_trimmed"
REPORT="$RESULT_DIR/02_fastp_reports"
mkdir -p "$OUT" "$REPORT"

# samples.tsv에서 run 목록을 뽑아 중복 없이 처리
RUNS=$(sample_rows | awk -F'\t' '{print $7"\n"$8}' | sort -u)
log "트리밍 대상 run: $(echo "$RUNS" | wc -l)개"

for run in $RUNS; do
    R1_IN="$FASTQ_DIR/${run}_1.fastq.gz"
    R2_IN="$FASTQ_DIR/${run}_2.fastq.gz"
    R1_OUT="$OUT/${run}_1.trim.fastq.gz"
    R2_OUT="$OUT/${run}_2.trim.fastq.gz"

    if [[ -s "$R1_OUT" && -s "$R2_OUT" ]]; then
        log "  $run — 이미 완료, 건너뜀"
        continue
    fi

    require_file "$R1_IN"; require_file "$R2_IN"
    log "  $run 트리밍 중..."

    # --cut_tail: STEP 1에서 R2의 3' 말단이 심하게 망가진 것이 확인됨
    #   cycle 150 → 평균 Phred 2.6, N 94.9% / cycle 145-149 → Phred 15, N 52.9%
    #   이걸 잘라내지 않으면 read 전체가 필터로 버려져 페어가 통째로 손실됨.
    #   슬라이딩 윈도우로 3' 말단만 제거해 read를 살린다.
    "$FASTP" \
        -i "$R1_IN" -I "$R2_IN" \
        -o "$R1_OUT" -O "$R2_OUT" \
        --detect_adapter_for_pe \
        --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 \
        --n_base_limit 5 \
        --qualified_quality_phred 15 \
        --unqualified_percent_limit 40 \
        --length_required 36 \
        --trim_poly_g \
        --thread 16 \
        --html "$REPORT/${run}_fastp.html" \
        --json "$REPORT/${run}_fastp.json" \
        2>> "$LOG_DIR/02_trim.log"
done

log "MultiQC (fastp 통합)"
"$MULTIQC" "$REPORT" -o "$REPORT" -n multiqc_fastp --force 2>&1 | tee -a "$LOG_DIR/02_trim.log"

log "완료 → $OUT (트리밍 fastq), $REPORT/multiqc_fastp.html"
