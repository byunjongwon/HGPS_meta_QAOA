#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────
# dupRadar — 발현량 대비 중복률
#
# RNA-seq에서 FastQC duplication FAIL은 대부분 고발현 유전자에서 오는
# 자연스러운 중복입니다. dupRadar는 발현량 축을 놓고 중복률을 그려
# 진짜 PCR duplicate(저발현 구간에서도 중복률이 높음)인지 구분합니다.
#
# 사용: 06b_dupradar.R <markdup_bam_dir> <gtf> <outdir> <threads>
# ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
    library(dupRadar); library(data.table); library(ggplot2)
})
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))),
                 "_plot_setup.R"))

a <- commandArgs(trailingOnly = TRUE)
BAMDIR <- a[1]; GTF <- a[2]; OUT <- a[3]; TH <- as.integer(a[4])
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

meta <- fread(Sys.getenv("SAMPLES"))

# strandedness는 STEP 06이 직접 만든 판정 결과에서 읽습니다.
# (STEP 07의 strandness_used.txt는 이 시점에 아직 없어서, 그걸 참조하면
#  stranded 라이브러리에서도 항상 0으로 떨어집니다.)
sfile <- file.path(Sys.getenv("RESULT_DIR"), "06_qc_align", "strandedness_summary.tsv")
STRAND <- 0L
if (file.exists(sfile)) {
    s <- fread(sfile)
    fwd <- mean(suppressWarnings(as.numeric(s[[3]])), na.rm = TRUE)
    rev <- mean(suppressWarnings(as.numeric(s[[4]])), na.rm = TRUE)
    STRAND <- if (!is.na(rev) && rev > 0.75) 2L else if (!is.na(fwd) && fwd > 0.75) 1L else 0L
    cat(sprintf("strandedness 판정: fwd=%.3f rev=%.3f -> %d\n", fwd, rev, STRAND))
} else {
    cat("strandedness 파일 없음 — 0(unstranded)으로 진행\n")
}

summ <- list()
for (s in meta$sample) {
    bam <- file.path(BAMDIR, paste0(s, ".markdup.bam"))
    if (!file.exists(bam)) { cat("  BAM 없음:", s, "\n"); next }
    dm <- tryCatch(analyzeDuprates(bam, GTF, STRAND, paired = TRUE, threads = TH),
                   error = function(e) { message("  실패 ", s, ": ", conditionMessage(e)); NULL })
    if (is.null(dm)) next

    png(file.path(OUT, sprintf("dupradar_%s.png", s)), 6, 5, units = "in", res = 300, type = "cairo")
    duprateExpDensPlot(DupMat = dm, main = s)
    dev.off()

    fit <- duprateExpFit(DupMat = dm)
    summ[[s]] <- data.table(sample = s, group = meta[sample == s, group],
                            intercept = fit$intercept, slope = fit$slope)
    cat(sprintf("  %-14s intercept %.2f%%  slope %.4f\n", s, fit$intercept, fit$slope))
}

if (length(summ)) {
    d <- rbindlist(summ)
    fwrite(d, file.path(OUT, "dupradar_summary.tsv"), sep = "\t")
    png(file.path(OUT, "dupradar_summary.png"), 7, 5, units = "in", res = 300, type = "cairo")
    print(ggplot(d, aes(group, intercept, color = group)) +
          geom_boxplot(outlier.shape = NA) + geom_jitter(width = .15, size = 3) +
          theme_bw(base_size = 12) +
          labs(y = "저발현 구간 중복률 (%)", x = NULL,
               title = "dupRadar — PCR duplicate 지표",
               subtitle = "intercept가 낮으면 중복은 대부분 고발현 유전자 유래(정상)"))
    dev.off()
    cat("\n완료 →", OUT, "\n")
}
