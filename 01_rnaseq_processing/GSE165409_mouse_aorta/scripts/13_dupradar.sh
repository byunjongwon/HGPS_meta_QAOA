#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────
# STEP 13. 중복 read 판별 — samtools markdup + dupRadar
#
# STEP 1의 FastQC에서 duplication FAIL이 다수 나왔는데, RNA-seq에서는
# 고발현 유전자 때문에 자연스럽게 생기는 중복이 대부분입니다.
# dupRadar는 "발현량 대비 중복률"을 그려 진짜 PCR duplicate인지
# 생물학적 중복인지 구분합니다.
#
# 산출 그림:
#   dupradar_{sample}.pdf   발현량 vs 중복률 산점도
#   dupradar_summary.pdf    18샘플 요약 (낮은 발현 구간 중복률)
#
# 예상 시간: 샘플당 5~8분 (markdup 포함)
# ─────────────────────────────────────────────────────────────
source "$(dirname "$0")/00_config.sh"

BAM_DIR="$RESULT_DIR/03_star"
OUT="$RESULT_DIR/13_dupradar"
MARKED="$OUT/marked_bam"
mkdir -p "$MARKED" "$OUT/figures"

STRAND=$(cat "$RESULT_DIR/05_counts/strandness_used.txt" 2>/dev/null || echo 2)
log "dupRadar 시작 (strandedness=$STRAND)"

sample_rows | while IFS=$'\t' read -r sample gsm group genotype age_group counts_col run1 run2; do
    IN="$BAM_DIR/${sample}.Aligned.sortedByCoord.out.bam"
    MK="$MARKED/${sample}.markdup.bam"
    [[ -s "$IN" ]] || { log "  $sample BAM 없음, 건너뜀"; continue; }
    if [[ -s "$MK" ]]; then log "  $sample — markdup 완료됨"; continue; fi

    log "  $sample — 중복 표시 중"
    # markdup은 name-sorted → fixmate → coord-sorted 순서를 요구
    "$SAMTOOLS" sort -n -@ "$THREADS" -m 2G -o "$OUT/tmp_n.bam" "$IN"
    "$SAMTOOLS" fixmate -m -@ "$THREADS" "$OUT/tmp_n.bam" "$OUT/tmp_f.bam"
    "$SAMTOOLS" sort -@ "$THREADS" -m 2G -o "$OUT/tmp_c.bam" "$OUT/tmp_f.bam"
    "$SAMTOOLS" markdup -@ "$THREADS" "$OUT/tmp_c.bam" "$MK"
    "$SAMTOOLS" index -@ "$THREADS" "$MK"
    rm -f "$OUT/tmp_n.bam" "$OUT/tmp_f.bam" "$OUT/tmp_c.bam"
done

log "dupRadar 분석 (R)"
"$RSCRIPT" - "$OUT" "$GTF" "$STRAND" "$THREADS" "$SAMPLES" <<'RS' 2>&1 | tee -a "$LOG_DIR/13_dupradar.log"
suppressPackageStartupMessages({ library(dupRadar); library(data.table); library(ggplot2) })
a <- commandArgs(trailingOnly = TRUE)
OUT <- a[1]; GTF <- a[2]; STRAND <- as.integer(a[3]); TH <- as.integer(a[4])
meta <- fread(a[5])
FIG <- file.path(OUT, "figures"); dir.create(FIG, showWarnings = FALSE)

summ <- list()
for (s in meta$sample) {
    bam <- file.path(OUT, "marked_bam", paste0(s, ".markdup.bam"))
    if (!file.exists(bam)) next
    dm <- tryCatch(analyzeDuprates(bam, GTF, STRAND, paired = TRUE, threads = TH),
                   error = function(e) { message("  실패 ", s, ": ", conditionMessage(e)); NULL })
    if (is.null(dm)) next
    pdf(file.path(FIG, sprintf("dupradar_%s.pdf", s)), 6, 5)
    duprateExpDensPlot(DupMat = dm, main = s)
    dev.off()
    fit <- duprateExpFit(DupMat = dm)
    summ[[s]] <- data.table(sample = s,
                            group = meta[sample == s, group],
                            intercept = fit$intercept,
                            slope = fit$slope)
    cat(sprintf("  %-16s intercept %.2f%%  slope %.4f\n", s, fit$intercept, fit$slope))
}
if (length(summ)) {
    d <- rbindlist(summ)
    fwrite(d, file.path(OUT, "dupradar_summary.tsv"), sep = "\t")
    pdf(file.path(FIG, "dupradar_summary.pdf"), 7, 5)
    print(ggplot(d, aes(group, intercept, color = group)) +
          geom_boxplot(outlier.shape = NA) + geom_jitter(width = .15, size = 3) +
          theme_bw(base_size = 12) +
          labs(y = "저발현 구간 중복률 (%)", x = NULL,
               title = "dupRadar — PCR duplicate 지표",
               subtitle = "intercept가 낮으면 중복은 대부분 고발현 유전자에서 유래(정상)"))
    dev.off()
}
RS

log "완료 → $OUT"
