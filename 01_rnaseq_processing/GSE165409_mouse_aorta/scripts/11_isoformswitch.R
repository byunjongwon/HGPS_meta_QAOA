#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────
# STEP 11. IsoformSwitchAnalyzeR — DTU (differential transcript usage)
#
# salmon 결과를 transcript 수준에서 재분석해, 유전자 발현은 그대로인데
# isoform 비율만 바뀌는 경우를 찾습니다. rMATS가 "event" 단위라면
# 이쪽은 "transcript" 단위 관점입니다.
#
# 산출 그림:
#   01_switch_summary.pdf       비교별 isoform switch 유전자 수
#   02_volcano_dIF.pdf          dIF volcano (x=dIF, y=-log10 q)
#   03_switch_vs_expression.pdf isoform switch vs 유전자 발현 변화
#   04_consequence_enrichment.pdf switch가 초래하는 결과(도메인 손실 등)
#   05_top_switch_plots/        유전자별 switch plot (구조+발현+ORF)
#   06_lmna_switch.pdf          ★ Lmna isoform switch
#
# 산출 표: switching_genes_{비교}.tsv
# ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
    library(IsoformSwitchAnalyzeR); library(data.table); library(ggplot2)
})

# UTF-8(Δ, 한글) 렌더링을 위한 cairo PDF 장치 설정
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)[1])), "_plot_setup.R"))

# [deposit] default PROJ = this dataset dir (scripts/..), still overridable by PROJ (was ~/lmna_rnaseq)
.this <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
PROJ <- Sys.getenv("PROJ", normalizePath(file.path(dirname(.this), "..")))
OUT  <- file.path(PROJ, "results/11_isoformswitch")
FIG  <- file.path(OUT, "figures")
SW   <- file.path(FIG, "05_top_switch_plots")
dir.create(SW, recursive = TRUE, showWarnings = FALSE)

meta   <- fread(file.path(PROJ, "meta/samples.tsv"))
QUANT  <- file.path(PROJ, "results/06_salmon/quant")
GTF    <- file.path(PROJ, "ref/gencode.vM25.primary_assembly.annotation.gtf")
TX_FA  <- file.path(PROJ, "ref/gencode.vM25.transcripts.fa")
# [deposit] the deposit ships the GENCODE files gzipped; use .gz when the uncompressed file is absent
if (!file.exists(GTF))   GTF   <- paste0(GTF, ".gz")
if (!file.exists(TX_FA)) TX_FA <- paste0(TX_FA, ".gz")

stopifnot(dir.exists(QUANT), file.exists(GTF), file.exists(TX_FA))

# ── salmon 결과 불러오기 ──────────────────────────────────────
# [deposit] the deposit omits salmon's Gibbs samples (aux_info/bootstrap/, ~250 MB). They are not
#   used by IsoformSwitchAnalyzeR, but tximport tries to read them because meta_info.json lists them
#   and stops if the file is missing. Only when they are absent, make tximport skip them
#   (dropInfReps = TRUE); counts/abundances are identical.
if (!all(file.exists(file.path(list.dirs(QUANT, recursive = FALSE),
                               "aux_info", "bootstrap", "bootstraps.gz")))) {
    .txi <- tximport::tximport
    formals(.txi)$dropInfReps <- TRUE
    utils::assignInNamespace("tximport", .txi, ns = "tximport")
    message("[deposit] salmon bootstraps not found -> tximport(dropInfReps = TRUE)")
}
cat("salmon quant 불러오는 중...\n")
quant <- importIsoformExpression(parentDir = QUANT, addIsofomIdAsColumn = TRUE)

design <- data.frame(sampleID = colnames(quant$abundance)[-1])
design$condition <- meta$group[match(design$sampleID, meta$sample)]
stopifnot(!any(is.na(design$condition)))
print(table(design$condition))

comparisons <- data.frame(
    condition_1 = c("LMNA",  "WTold", "LMNA"),
    condition_2 = c("WTyng", "WTyng", "WTold")
)

cat("SwitchAnalyzeRlist 생성 중...\n")
sw <- importRdata(
    isoformCountMatrix   = quant$counts,
    isoformRepExpression = quant$abundance,
    designMatrix         = design,
    comparisonsToMake    = comparisons,
    isoformExonAnnoation = GTF,
    isoformNtFasta       = TX_FA,
    showProgress         = FALSE,
    ignoreAfterPeriod    = FALSE,
    fixStringTieAnnotationProblem = FALSE
)

cat("필터링 + DTU 검정 (satuRn)...\n")
sw <- preFilter(sw, geneExpressionCutoff = 1, isoformExpressionCutoff = 0,
                removeSingleIsoformGenes = TRUE)
sw <- isoformSwitchTestSatuRn(sw, reduceToSwitchingGenes = FALSE)

cat("ORF 예측 + switch 결과 분석...\n")
sw <- analyzeORF(sw, orfMethod = "longest", showProgress = FALSE)
sw <- analyzeAlternativeSplicing(sw, quiet = TRUE)

# analyzeSwitchConsequences의 "ORF_seq_similarity"는 아미노산 서열을 요구합니다.
# analyzeORF()는 ORF 좌표만 넣어주므로 extractSequence()로 실제 서열을 뽑아
# switchAnalyzeRlist에 채워 넣어야 합니다. (없으면 그 자리에서 에러)
cat("ORF 서열 추출 중...\n")
sw <- extractSequence(sw, onlySwitchingGenes = TRUE,
                      extractNTseq = TRUE, extractAAseq = TRUE,
                      removeShortAAseq = TRUE, removeLongAAseq = FALSE,
                      addToSwitchAnalyzeRlist = TRUE, writeToFile = FALSE)

sw <- analyzeSwitchConsequences(
    sw,
    consequencesToAnalyze = c("intron_retention", "ORF_seq_similarity",
                              "NMD_status", "tss", "tts"),
    showProgress = FALSE
)

saveRDS(sw, file.path(OUT, "switchlist.rds"))

# ── 표 저장 ───────────────────────────────────────────────────
res <- as.data.table(sw$isoformFeatures)
fwrite(res, file.path(OUT, "isoform_features_all.tsv"), sep = "\t")
for (cp in unique(res$condition_1)) {
    sub <- res[condition_1 == cp]
    for (c2 in unique(sub$condition_2)) {
        d <- sub[condition_2 == c2 & isoform_switch_q_value < 0.05 & abs(dIF) > 0.1]
        fwrite(d[order(isoform_switch_q_value)],
               file.path(OUT, sprintf("switching_genes_%s_vs_%s.tsv", cp, c2)), sep = "\t")
    }
}

# ── 01 switch summary ─────────────────────────────────────────
pdf(file.path(FIG, "01_switch_summary.pdf"), 8, 5)
print(extractSwitchSummary(sw, filterForConsequences = FALSE))
tryCatch(print(extractSwitchSummary(sw, filterForConsequences = TRUE)),
         error = function(e) message("  consequence 요약 생략: ", conditionMessage(e)))
dev.off()

# ── 02 dIF volcano ────────────────────────────────────────────
res[, comp := paste0(condition_1, "_vs_", condition_2)]
res[, sig := !is.na(isoform_switch_q_value) &
              isoform_switch_q_value < 0.05 & abs(dIF) > 0.1]
pdf(file.path(FIG, "02_volcano_dIF.pdf"), 11, 4.5)
print(ggplot(res[!is.na(isoform_switch_q_value)],
             aes(dIF, -log10(pmax(isoform_switch_q_value, 1e-300)), color = sig)) +
      geom_point(size = 0.7, alpha = 0.5) +
      geom_vline(xintercept = c(-0.1, 0.1), lty = 2, col = "grey40") +
      geom_hline(yintercept = -log10(0.05), lty = 2, col = "grey40") +
      facet_wrap(~comp, nrow = 1) +
      scale_color_manual(values = c(`FALSE` = "grey75", `TRUE` = "#D55E00")) +
      theme_bw(base_size = 11) +
      labs(x = "dIF (differential isoform fraction)", y = "-log10(q)",
           title = "Differential transcript usage"))
dev.off()

# ── 03 switch vs 유전자 발현 ──────────────────────────────────
pdf(file.path(FIG, "03_switch_vs_expression.pdf"), 11, 4.5)
print(ggplot(res[!is.na(gene_log2_fold_change) & !is.na(dIF)],
             aes(gene_log2_fold_change, dIF, color = sig)) +
      geom_point(size = 0.7, alpha = 0.5) +
      geom_hline(yintercept = 0, lty = 2, col = "grey50") +
      geom_vline(xintercept = 0, lty = 2, col = "grey50") +
      facet_wrap(~comp, nrow = 1) +
      scale_color_manual(values = c(`FALSE` = "grey75", `TRUE` = "#0072B2")) +
      theme_bw(base_size = 11) +
      labs(x = "gene log2 fold change", y = "dIF",
           title = "발현 변화와 isoform 사용 변화는 독립인가"))
dev.off()

# ── 04 consequence enrichment ─────────────────────────────────
tryCatch({
    pdf(file.path(FIG, "04_consequence_enrichment.pdf"), 8, 6)
    print(extractConsequenceEnrichment(sw, returnResult = FALSE))
    dev.off()
}, error = function(e) message("  consequence enrichment 생략: ", conditionMessage(e)))

# ── 05 상위 switch plot ───────────────────────────────────────
top <- res[sig == TRUE][order(isoform_switch_q_value)]
top_genes <- unique(top$gene_name)
top_genes <- head(top_genes[!is.na(top_genes)], 15)
for (g in top_genes) {
    tryCatch({
        pdf(file.path(SW, sprintf("switch_%s.pdf", gsub("[^A-Za-z0-9_.-]", "_", g))), 10, 7)
        switchPlot(sw, gene = g)
        dev.off()
    }, error = function(e) { try(dev.off(), silent = TRUE)
        message("  switchPlot 실패 (", g, "): ", conditionMessage(e)) })
}
cat("switch plot", length(top_genes), "개 생성\n")

# ── 06 ★ Lmna ─────────────────────────────────────────────────
if ("Lmna" %in% res$gene_name) {
    tryCatch({
        pdf(file.path(FIG, "06_lmna_switch.pdf"), 10, 7)
        switchPlot(sw, gene = "Lmna")
        dev.off()
        cat("★ Lmna switch plot 생성\n")
    }, error = function(e) { try(dev.off(), silent = TRUE)
        message("  Lmna switchPlot 실패: ", conditionMessage(e)) })
    fwrite(res[gene_name == "Lmna"], file.path(OUT, "lmna_isoforms.tsv"), sep = "\t")
}

cat("완료 → ", FIG, "\n")
