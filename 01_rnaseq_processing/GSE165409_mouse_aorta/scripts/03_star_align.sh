#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────
# STEP 3. STAR 정렬 — 샘플당 run 2개를 한 번에 넣어 병합
#
# 옵션 근거:
#   ENCODE long-RNA 표준 + rMATS 권장 조합
#   --outSAMstrandField intronMotif  : 라이브러리가 unstranded일 때
#                                      XS 태그를 넣어 StringTie/rMATS 호환
#   --quantMode GeneCounts           : STAR가 직접 유전자 카운트 생성
#                                      (featureCounts와 교차검증용)
#                TranscriptomeSAM    : salmon(alignment mode)/isoform 분석용
#   --outSAMattributes NH HI AS nM NM MD jM jI MC ch
#                                      : rMATS/sashimi에 필요한 junction 태그 포함
#
# 산출:
#   {sample}.Aligned.sortedByCoord.out.bam (+ .bai)
#   {sample}.Aligned.toTranscriptome.out.bam
#   {sample}.ReadsPerGene.out.tab
#   {sample}.SJ.out.tab              ← junction 그림 재료
#   {sample}.Log.final.out           ← mapping statistics 그림 재료
#
# 예상 시간: 샘플당 8~15분 x 18 = 3~4시간
# 디스크: BAM 약 90 G
# ⚠️ 메모리: STAR 인덱스 로딩에 약 32 G. 동시 실행 금지 (순차 처리).
# ─────────────────────────────────────────────────────────────
source "$(dirname "$0")/00_config.sh"

# 트리밍 결과를 쓸지, 원본을 쓸지 선택 (기본: 트리밍)
INPUT_DIR="${INPUT_DIR:-$RESULT_DIR/02_trimmed}"
SUFFIX="${SUFFIX:-.trim.fastq.gz}"

OUT="$RESULT_DIR/03_star"
mkdir -p "$OUT"

require_file "$STAR_INDEX/SA"
require_file "$STAR_INDEX/SAindex"

log "STAR 정렬 시작 — 인덱스: $STAR_INDEX"
log "입력: $INPUT_DIR"

sample_rows | while IFS=$'\t' read -r sample gsm group genotype age_group counts_col run1 run2; do
    BAM="$OUT/${sample}.Aligned.sortedByCoord.out.bam"
    if [[ -s "$BAM" && -s "${BAM}.bai" ]]; then
        log "  $sample — 이미 완료, 건너뜀"
        continue
    fi

    R1="$INPUT_DIR/${run1}_1${SUFFIX},$INPUT_DIR/${run2}_1${SUFFIX}"
    R2="$INPUT_DIR/${run1}_2${SUFFIX},$INPUT_DIR/${run2}_2${SUFFIX}"

    for f in ${R1//,/ } ${R2//,/ }; do require_file "$f"; done

    log "  $sample ($group) 정렬 중 — run: $run1 + $run2"

    "$STAR" \
        --runThreadN "$THREADS" \
        --genomeDir "$STAR_INDEX" \
        --readFilesIn "$R1" "$R2" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$OUT/${sample}." \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS nM NM MD jM jI MC ch \
        --outSAMstrandField intronMotif \
        --outSAMattrIHstart 0 \
        --quantMode TranscriptomeSAM GeneCounts \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outSAMattrRGline ID:"$sample" SM:"$sample" LB:"$sample" PL:ILLUMINA \
        2>&1 | tee -a "$LOG_DIR/03_star_align.log"

    "$SAMTOOLS" index -@ "$THREADS" "$BAM"
    log "  $sample 완료 — $(du -h "$BAM" | cut -f1)"
done

log "정렬 통계 요약 (MultiQC가 Log.final.out을 읽음)"
"$MULTIQC" "$OUT" -o "$OUT" -n multiqc_star --force 2>&1 | tee -a "$LOG_DIR/03_star_align.log"

log "완료 → $OUT"
