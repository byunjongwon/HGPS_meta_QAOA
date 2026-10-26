#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────
# STEP 4. 정렬 후 QC — counts 행렬로는 절대 못 만드는 그림들
#
# 산출 그림:
#   [RSeQC]
#     bam_stat.py            → mapping statistics
#     infer_experiment.py    → 라이브러리 방향성 판정 (★ 이후 단계 파라미터 결정)
#     read_distribution.py   → exon/intron/UTR/intergenic stacked bar
#     geneBody_coverage.py   → 5'→3' 커버리지 곡선 (RNA 분해 검출)
#     junction_annotation.py → known vs novel junction pie
#     junction_saturation.py → 시퀀싱 깊이가 splicing 분석에 충분한지 (★ 필수 근거)
#     inner_distance.py      → paired-end fragment 길이 분포
#   [Qualimap]
#     rnaseq 모듈           → 종합 BAM QC 리포트
#   [MultiQC]
#     위 전부를 한 장으로 통합
#
# 예상 시간: 1~2시간 (geneBody_coverage가 대부분)
# ─────────────────────────────────────────────────────────────
source "$(dirname "$0")/00_config.sh"

BAM_DIR="$RESULT_DIR/03_star"
OUT="$RESULT_DIR/04_qc_align"
mkdir -p "$OUT"/{bam_stat,infer,read_dist,genebody,junction,inner_dist,qualimap}

BED12_SUB="$REF_DIR/gencode.vM25.genebody2000.bed12"

# ── BED12 준비 ────────────────────────────────────────────────
if [[ ! -s "$BED12" ]]; then
    log "BED12 생성 중 (GTF → BED12, 수 분 소요)"
    "$RSCRIPT" "$(dirname "$0")/03b_make_bed12.R" "$GTF" "$BED12" "$BED12_SUB" \
        2>&1 | tee -a "$LOG_DIR/04_qc_align.log"
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

# ── geneBody coverage (축소 BED로 전체 샘플 한 번에) ──────────
if [[ ! -s "$OUT/genebody/all.geneBodyCoverage.curves.pdf" ]]; then
    log "geneBody coverage (18샘플 일괄, 오래 걸림)"
    printf "%s\n" "${BAMS[@]}" > "$OUT/genebody/bam_list.txt"
    "$RSEQC_BIN/geneBody_coverage.py" -i "$OUT/genebody/bam_list.txt" \
        -r "$BED12_SUB" -o "$OUT/genebody/all" 2>&1 | tee -a "$LOG_DIR/04_qc_align.log"
fi

# ── Qualimap ──────────────────────────────────────────────────
# 샘플당 약 11분(단일 스레드, BAM 전면 스캔)이라 18샘플에 3시간이 걸립니다.
# Qualimap이 주는 정보(exon/intron 분포, 5'-3' bias, 커버리지 균일도)는
# 위의 read_distribution.py / geneBody_coverage.py 결과와 사실상 중복이라
# 기본적으로 건너뜁니다. 필요하면 SKIP_QUALIMAP=0 으로 실행하세요.
SKIP_QUALIMAP="${SKIP_QUALIMAP:-1}"
if [[ "$SKIP_QUALIMAP" == "1" ]]; then
    log "Qualimap 건너뜀 (SKIP_QUALIMAP=1). RSeQC 결과로 대체합니다."
else
    for bam in "${BAMS[@]}"; do
        s=$(basename "$bam" .Aligned.sortedByCoord.out.bam)
        [[ -s "$OUT/qualimap/$s/qualimapReport.html" ]] && continue
        log "  [$s] Qualimap"
        unset DISPLAY
        "$QUALIMAP" rnaseq -bam "$bam" -gtf "$GTF" -outdir "$OUT/qualimap/$s" \
            -pe --java-mem-size=16G 2>&1 | tail -2 | tee -a "$LOG_DIR/04_qc_align.log"
    done
fi

# ── 방향성 판정 요약 (★ 이후 단계에 그대로 쓰임) ──────────────
log "────────── 라이브러리 방향성 판정 ──────────"
{
    echo "sample	failed	fwd(1++,1--,2+-,2-+)	rev(1+-,1-+,2++,2--)"
    for f in "$OUT"/infer/*.txt; do
        s=$(basename "$f" .txt)
        fail=$(grep -oP 'Fraction of reads failed.*: \K[0-9.]+' "$f" || echo NA)
        fwd=$(grep -oP '"1\+\+,1--,2\+-,2-\+":\s*\K[0-9.]+' "$f" || echo NA)
        rev=$(grep -oP '"1\+-,1-\+,2\+\+,2--":\s*\K[0-9.]+' "$f" || echo NA)
        echo "$s	$fail	$fwd	$rev"
    done
} | tee "$OUT/strandedness_summary.tsv" | column -t

cat <<'NOTE'

── 해석 방법 ───────────────────────────────────────────────
  fwd ≈ rev ≈ 0.5       → unstranded
      featureCounts -s 0 / salmon -l IU / rMATS --libType fr-unstranded
  rev >> fwd (0.9 이상) → reverse stranded (dUTP, 가장 흔함)
      featureCounts -s 2 / salmon -l ISR / rMATS --libType fr-firststrand
  fwd >> rev (0.9 이상) → forward stranded
      featureCounts -s 1 / salmon -l ISF / rMATS --libType fr-secondstrand

  판정 결과를 00_config.sh 의 STRANDNESS 변수에 반영하세요.
────────────────────────────────────────────────────────────
NOTE

log "MultiQC 통합"
"$MULTIQC" "$OUT" "$BAM_DIR" -o "$OUT" -n multiqc_qc_align --force 2>&1 | tee -a "$LOG_DIR/04_qc_align.log"

log "완료 → $OUT/multiqc_qc_align.html"
