#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────
# STEP 06 — 정렬 후 QC (counts 행렬로는 만들 수 없는 지표들)
#
#   [RSeQC]  bam_stat / infer_experiment(★방향성 판정) / read_distribution
#            geneBody_coverage(5'-3' bias) / junction_annotation
#            junction_saturation(★splicing 분석 깊이 근거) / inner_distance
#   [Picard] CollectRnaSeqMetrics / MarkDuplicates
#   [dupRadar] 발현량 대비 중복률 — PCR 중복 vs 자연 중복 구분
#   [Qualimap] 종합 BAM QC (SKIP_QUALIMAP=0 일 때만)
#   [MultiQC] 전부 통합
# ─────────────────────────────────────────────────────────────
source "$(dirname "$0")/00_config.sh"

BAM_DIR="$RESULT_DIR/05_star"
OUT="$RESULT_DIR/06_qc_align"
mkdir -p "$OUT"/{bam_stat,infer,read_dist,genebody,junction,inner_dist,qualimap,picard,dupradar}

# ── BED12 준비 ────────────────────────────────────────────────
if [[ ! -s "$BED12" ]]; then
    log "BED12 생성 중 (GTF → BED12)"
    "$RSCRIPT" "$(dirname "$0")/02b_make_bed12.R" "$GTF" "$BED12" "$BED12_HK" \
        2>&1 | tee -a "$LOG_DIR/06_qc_align.log"
fi
require_file "$BED12"

BAMS=("$BAM_DIR"/*.Aligned.sortedByCoord.out.bam)
log "대상 BAM ${#BAMS[@]}개"

# ── 샘플별 RSeQC ──────────────────────────────────────────────
for bam in "${BAMS[@]}"; do
    s=$(basename "$bam" .Aligned.sortedByCoord.out.bam)
    log "  [$s] RSeQC"
    [[ -s "$OUT/bam_stat/${s}.txt" ]]  || "$RSEQC_BIN/bam_stat.py"  -i "$bam" > "$OUT/bam_stat/${s}.txt" 2>/dev/null
    [[ -s "$OUT/infer/${s}.txt" ]]     || "$RSEQC_BIN/infer_experiment.py" -i "$bam" -r "$BED12" > "$OUT/infer/${s}.txt" 2>/dev/null
    [[ -s "$OUT/read_dist/${s}.txt" ]] || "$RSEQC_BIN/read_distribution.py" -i "$bam" -r "$BED12" > "$OUT/read_dist/${s}.txt" 2>/dev/null
    [[ -s "$OUT/junction/${s}.junction.xls" ]] || \
        "$RSEQC_BIN/junction_annotation.py" -i "$bam" -r "$BED12" -o "$OUT/junction/${s}" 2>/dev/null
    [[ -s "$OUT/junction/${s}.junctionSaturation_plot.pdf" ]] || \
        "$RSEQC_BIN/junction_saturation.py" -i "$bam" -r "$BED12" -o "$OUT/junction/${s}" 2>/dev/null
    [[ -s "$OUT/inner_dist/${s}.inner_distance_plot.pdf" ]] || \
        "$RSEQC_BIN/inner_distance.py" -i "$bam" -r "$BED12" -o "$OUT/inner_dist/${s}" 2>/dev/null
done

# ── geneBody coverage (축소 BED로 전체 샘플 일괄) ─────────────
if [[ ! -s "$OUT/genebody/all.geneBodyCoverage.curves.pdf" ]]; then
    log "geneBody coverage (8샘플 일괄)"
    printf "%s\n" "${BAMS[@]}" > "$OUT/genebody/bam_list.txt"
    "$RSEQC_BIN/geneBody_coverage.py" -i "$OUT/genebody/bam_list.txt" \
        -r "$BED12_HK" -o "$OUT/genebody/all" 2>&1 | tail -5 | tee -a "$LOG_DIR/06_qc_align.log"
fi

# ── 방향성 판정 (★ Picard 및 이후 단계 파라미터) ─────────────
log "────────── 라이브러리 방향성 판정 ──────────"
{
    echo -e "sample\tfailed\tfwd(1++,1--,2+-,2-+)\trev(1+-,1-+,2++,2--)"
    for f in "$OUT"/infer/*.txt; do
        s=$(basename "$f" .txt)
        fail=$(grep -oP 'Fraction of reads failed.*: \K[0-9.]+' "$f" || echo NA)
        fwd=$(grep -oP '"1\+\+,1--,2\+-,2-\+":\s*\K[0-9.]+' "$f" || echo NA)
        rev=$(grep -oP '"1\+-,1-\+,2\+\+,2--":\s*\K[0-9.]+' "$f" || echo NA)
        echo -e "$s\t$fail\t$fwd\t$rev"
    done
} | tee "$OUT/strandedness_summary.tsv" | column -t

read -r AVG_FWD AVG_REV < <(awk -F'\t' 'NR>1 && $3!="NA" {f+=$3; r+=$4; n++} END{printf "%.3f %.3f\n", f/n, r/n}' "$OUT/strandedness_summary.tsv")
if   (( $(echo "$AVG_REV > 0.75" | bc -l) )); then PICARD_STRAND="SECOND_READ_TRANSCRIPTION_STRAND"
elif (( $(echo "$AVG_FWD > 0.75" | bc -l) )); then PICARD_STRAND="FIRST_READ_TRANSCRIPTION_STRAND"
else PICARD_STRAND="NONE"; fi
log "평균 fwd=$AVG_FWD rev=$AVG_REV → Picard STRAND_SPECIFICITY=$PICARD_STRAND"

cat <<'NOTE'

── 해석 ────────────────────────────────────────────────────
  fwd ≈ rev ≈ 0.5       → unstranded
      featureCounts -s 0 / salmon -l IU / rMATS --libType fr-unstranded
  rev >> fwd (0.75 이상) → reverse stranded (dUTP, 가장 흔함)
      featureCounts -s 2 / salmon -l ISR / rMATS --libType fr-firststrand
  fwd >> rev (0.75 이상) → forward stranded
      featureCounts -s 1 / salmon -l ISF / rMATS --libType fr-secondstrand
  → 이후 STEP 07/08/13 이 이 파일을 읽어 자동 반영합니다.
────────────────────────────────────────────────────────────
NOTE

# ── Picard ────────────────────────────────────────────────────
# refFlat / rRNA interval 준비
REFFLAT="$REF_DIR/gencode.v${GENCODE_REL}.refFlat.txt"
if [[ ! -s "$REFFLAT" ]]; then
    log "refFlat 생성 (BED12 → refFlat)"
    # BED12의 blockSizes(11)/blockStarts(12)는 txStart 기준 상대값입니다.
    # Picard refFlat은 절대 좌표의 exonStarts/exonEnds를 요구하므로 변환합니다.
    awk 'BEGIN{OFS="\t"}
    {
        n = split($11, sz, ","); split($12, st, ",")
        starts = ""; ends = ""
        for (i = 1; i <= $10; i++) {
            s = $2 + st[i]; e = s + sz[i]
            starts = starts s ","
            ends   = ends   e ","
        }
        print $4, $4, $1, $6, $2, $3, $7, $8, $10, starts, ends
    }' "$BED12" > "$REFFLAT"
fi

RRNA_INTERVAL="$REF_DIR/rRNA.interval_list"
if [[ ! -s "$RRNA_INTERVAL" ]]; then
    log "rRNA interval_list 생성"
    "$SAMTOOLS" view -H "${BAMS[0]}" | grep '^@SQ' > "$RRNA_INTERVAL"
    awk -F'\t' '$3=="gene" && /gene_type "rRNA(_pseudogene)?"/ {
        match($9, /gene_name "[^"]+"/); g=substr($9, RSTART+11, RLENGTH-12);
        print $1"\t"$4"\t"$5"\t"$7"\t"g
    }' "$GTF" >> "$RRNA_INTERVAL"
fi

# ⚠️ conda의 picard 래퍼는 힙을 -Xmx2g로 고정합니다.
#    GENCODE v47 전체(387,944 transcript) 모델을 올리는 CollectRnaSeqMetrics는
#    2 GB로는 OOM으로 죽습니다(샘플당 10분 낭비 후 실패).
#    래퍼가 -Xm* 인자를 JVM 옵션으로 넘겨주므로 여기서 명시적으로 키웁니다.
PICARD_MEM="${PICARD_MEM:--Xmx24g}"

for bam in "${BAMS[@]}"; do
    s=$(basename "$bam" .Aligned.sortedByCoord.out.bam)
    if [[ ! -s "$OUT/picard/${s}.RnaSeqMetrics.txt" ]]; then
        log "  [$s] Picard CollectRnaSeqMetrics"
        "$PICARD" "$PICARD_MEM" CollectRnaSeqMetrics \
            -I "$bam" -O "$OUT/picard/${s}.RnaSeqMetrics.txt" \
            --REF_FLAT "$REFFLAT" --STRAND_SPECIFICITY "$PICARD_STRAND" \
            --RIBOSOMAL_INTERVALS "$RRNA_INTERVAL" \
            --CHART_OUTPUT "$OUT/picard/${s}.coverage.pdf" \
            2>&1 | tail -2 | tee -a "$LOG_DIR/06_qc_align.log" || log "  [$s] Picard 실패 (건너뜀)"
    fi
    if [[ ! -s "$OUT/picard/${s}.markdup.txt" ]]; then
        log "  [$s] Picard MarkDuplicates"
        "$PICARD" "$PICARD_MEM" MarkDuplicates \
            -I "$bam" -O "$OUT/picard/${s}.markdup.bam" \
            -M "$OUT/picard/${s}.markdup.txt" \
            --ASSUME_SORT_ORDER coordinate --REMOVE_DUPLICATES false \
            --TMP_DIR "$TMP_DIR" \
            2>&1 | tail -2 | tee -a "$LOG_DIR/06_qc_align.log" || log "  [$s] MarkDuplicates 실패 (건너뜀)"
    fi
done

# ── dupRadar (markdup BAM 필요) ──────────────────────────────
if compgen -G "$OUT/picard/*.markdup.bam" > /dev/null; then
    log "dupRadar 실행"
    "$RSCRIPT" "$(dirname "$0")/06b_dupradar.R" "$OUT/picard" "$GTF" "$OUT/dupradar" "$THREADS" \
        2>&1 | tail -10 | tee -a "$LOG_DIR/06_qc_align.log" || log "dupRadar 실패 (건너뜀)"
fi

# ── Qualimap (기본 건너뜀 — RSeQC와 정보 중복, 샘플당 10분+) ──
SKIP_QUALIMAP="${SKIP_QUALIMAP:-1}"
if [[ "$SKIP_QUALIMAP" == "1" ]]; then
    log "Qualimap 건너뜀 (SKIP_QUALIMAP=1)"
else
    for bam in "${BAMS[@]}"; do
        s=$(basename "$bam" .Aligned.sortedByCoord.out.bam)
        [[ -s "$OUT/qualimap/$s/qualimapReport.html" ]] && continue
        log "  [$s] Qualimap"
        unset DISPLAY
        "$QUALIMAP" rnaseq -bam "$bam" -gtf "$GTF" -outdir "$OUT/qualimap/$s" \
            -pe --java-mem-size=16G 2>&1 | tail -2 | tee -a "$LOG_DIR/06_qc_align.log"
    done
fi

log "MultiQC 통합"
"$MULTIQC" "$OUT" "$BAM_DIR" "$RESULT_DIR/04_fastp_reports" -o "$OUT" \
    -n multiqc_qc_align --force 2>&1 | tail -3 | tee -a "$LOG_DIR/06_qc_align.log"

log "완료 → $OUT/multiqc_qc_align.html"

# ---------------------------------------------------------------------------
# PDF -> PNG 변환 (RSeQC / sashimi 등 외부 도구가 만든 PDF를 300dpi PNG 로 바꾸고 삭제)
# 사용자는 PowerPoint 삽입을 위해 PNG 만 사용한다.
# ---------------------------------------------------------------------------
# [deposit] convert PDFs under $RESULT_DIR (was scripts/../results, i.e. wrong tree when PROJ is overridden)
find "${RESULTS:-$RESULT_DIR}" -name "*.pdf" 2>/dev/null | while read -r f; do
  b="${f%.pdf}"
  if [ "$(pdfinfo "$f" 2>/dev/null | awk '/^Pages:/{print $2}')" = "1" ]; then
    pdftoppm -png -r 300 -f 1 -l 1 -singlefile "$f" "$b" && rm -f "$f"
  else
    pdftoppm -png -r 300 "$f" "$b" && rm -f "$f"
  fi
done
