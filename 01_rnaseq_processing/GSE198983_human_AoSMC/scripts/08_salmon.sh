#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────
# STEP 08 — Salmon transcript 수준 정량 (selective alignment)
#
# featureCounts(유전자 수준)와 별개로 transcript 수준 추정치를 만듭니다.
#   · DTU(IsoformSwitchAnalyzeR) 입력
#   · tximport → DESeq2 유전자 수준 재집계 (길이 보정, featureCounts 교차검증)
#
# --validateMappings  selective alignment (정확도 ↑, 기본값이지만 명시)
# --seqBias --gcBias  서열/GC 편향 보정
# --numGibbsSamples   DTU 불확실성 추정용 posterior sample
# ─────────────────────────────────────────────────────────────
source "$(dirname "$0")/00_config.sh"

IN_DIR="$RESULT_DIR/04_trimmed"
OUT="$RESULT_DIR/08_salmon"
mkdir -p "$OUT"

require_file "$SALMON_INDEX/info.json"

# ── 라이브러리 타입 자동 판정 ────────────────────────────────
STRAND=$(cat "$RESULT_DIR/07_counts/strandness_used.txt" 2>/dev/null || echo 0)
case "$STRAND" in
    2) LIBTYPE="ISR" ;;
    1) LIBTYPE="ISF" ;;
    *) LIBTYPE="A"   ;;   # A = 자동 판정
esac
log "Salmon libType: $LIBTYPE (featureCounts -s $STRAND 기준)"

sample_rows | while IFS=$'\t' read -r sample gsm srr group replicate; do
    SDIR="$OUT/$sample"
    if [[ -s "$SDIR/quant.sf" ]]; then
        log "  $sample — 이미 완료, 건너뜀"
        continue
    fi
    R1="$IN_DIR/${srr}_1.trim.fastq.gz"
    R2="$IN_DIR/${srr}_2.trim.fastq.gz"
    require_file "$R1"; require_file "$R2"

    log "  $sample ($group) 정량 중"
    "$SALMON" quant \
        -i "$SALMON_INDEX" \
        -l "$LIBTYPE" \
        -1 "$R1" -2 "$R2" \
        -o "$SDIR" \
        -p "$THREADS" \
        --validateMappings \
        --seqBias --gcBias \
        --numGibbsSamples 30 \
        2>&1 | tail -5 | tee -a "$LOG_DIR/08_salmon.log"
done

# ── tx2gene 매핑 (GTF에서 추출) ──────────────────────────────
T2G="$REF_DIR/tx2gene.v${GENCODE_REL}.tsv"
if [[ ! -s "$T2G" ]]; then
    log "tx2gene 생성"
    awk -F'\t' '$3=="transcript"{
        match($9, /transcript_id "[^"]+"/); t=substr($9, RSTART+15, RLENGTH-16);
        match($9, /gene_id "[^"]+"/);       g=substr($9, RSTART+9,  RLENGTH-10);
        match($9, /gene_name "[^"]+"/);     n=(RSTART>0)?substr($9, RSTART+11, RLENGTH-12):g;
        print t"\t"g"\t"n
    }' "$GTF" > "$T2G"
fi
log "tx2gene: $(wc -l < "$T2G") transcript"

log "MultiQC (Salmon 매핑률)"
"$MULTIQC" "$OUT" -o "$OUT" -n multiqc_salmon --force 2>&1 | tail -3 | tee -a "$LOG_DIR/08_salmon.log"

log "완료 → $OUT"
