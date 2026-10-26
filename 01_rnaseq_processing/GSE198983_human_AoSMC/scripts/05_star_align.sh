#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────
# STEP 05 — STAR 정렬 (샘플당 run 1개)
#
# 옵션 근거: ENCODE long-RNA 표준 + rMATS 권장 조합
#   --outSAMstrandField intronMotif  unstranded일 때 XS 태그 부여 (StringTie/rMATS 호환)
#   --quantMode GeneCounts           STAR 자체 카운트 (featureCounts 교차검증용)
#              TranscriptomeSAM      isoform 분석용 transcriptome BAM
#   --outSAMattributes ... jM jI      rMATS/sashimi에 필요한 junction 태그
#
# 산출: {sample}.Aligned.sortedByCoord.out.bam(+.bai), .toTranscriptome.out.bam,
#       .ReadsPerGene.out.tab, .SJ.out.tab, .Log.final.out
# ⚠️ 메모리: 인덱스 로딩에 ~32G. 순차 처리(동시 실행 금지).
# ─────────────────────────────────────────────────────────────
source "$(dirname "$0")/00_config.sh"

INPUT_DIR="${INPUT_DIR:-$RESULT_DIR/04_trimmed}"
SUFFIX="${SUFFIX:-.trim.fastq.gz}"

OUT="$RESULT_DIR/05_star"
mkdir -p "$OUT"

require_file "$STAR_INDEX/SA"
require_file "$STAR_INDEX/SAindex"

log "STAR 정렬 시작 — 인덱스: $STAR_INDEX / 입력: $INPUT_DIR"

sample_rows | while IFS=$'\t' read -r sample gsm srr group replicate; do
    BAM="$OUT/${sample}.Aligned.sortedByCoord.out.bam"
    if [[ -s "$BAM" && -s "${BAM}.bai" ]]; then
        log "  $sample — 이미 완료, 건너뜀"
        continue
    fi

    R1="$INPUT_DIR/${srr}_1${SUFFIX}"
    R2="$INPUT_DIR/${srr}_2${SUFFIX}"
    require_file "$R1"; require_file "$R2"

    log "  $sample ($group) 정렬 중 — $srr"
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
        2>&1 | tee -a "$LOG_DIR/05_star_align.log"

    "$SAMTOOLS" index -@ "$THREADS" "$BAM"
    log "  $sample 완료 — $(du -h "$BAM" | cut -f1)"
done

log "MultiQC (Log.final.out 기반 정렬 통계)"
"$MULTIQC" "$OUT" -o "$OUT" -n multiqc_star --force 2>&1 | tail -3 | tee -a "$LOG_DIR/05_star_align.log"

log "완료 → $OUT"
