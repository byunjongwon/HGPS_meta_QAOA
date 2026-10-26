#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────
# STEP 16 — IsoformSwitchAnalyzeR: 차등 isoform 사용(DTU)
#
# rMATS는 "이 exon이 더/덜 포함되는가"를 봅니다.
# DTU는 "유전자 총 발현은 그대로인데 isoform 비율이 바뀌었는가"를 봅니다.
# 둘은 상보적입니다.
#
# 산출:
#   switch 요약표, isoform fraction(IF) 변화, switch plot,
#   기능적 귀결(ORF 손실, NMD 민감성, domain 손실) 분석
#
# ⚠️ 외부 도구(CPAT/Pfam/SignalP)는 웹서버 제출이 필요해 생략합니다.
#    → analyzeSwitchConsequences 중 서열 기반 항목은 제외하고,
#      ORF/NMD/intron retention 등 내부 계산 항목만 씁니다.
# ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
    library(IsoformSwitchAnalyzeR); library(data.table); library(ggplot2)
})
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])),
                 "_plot_setup.R"))

# [deposit] PROJ default = dataset dir (parent of this scripts/ dir; was ~/progerin_rnaseq); overridable by PROJ
PROJ    <- Sys.getenv("PROJ", normalizePath(file.path(dirname(sub("--file=", "",
               grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "..")))
RESULTS <- file.path(PROJ, "results")
REF     <- Sys.getenv("REF_DIR", file.path(PROJ, "ref"))
GTF     <- Sys.getenv("GTF")
TXFA    <- Sys.getenv("TXFA")
# [deposit] the deposit ships only gencode.v47.*.gz: prefer the uncompressed file (set by 00_config.sh)
#           if present, else use the .gz next to it (rtracklayer/Biostrings read .gz directly)
if (!nzchar(GTF))  GTF  <- file.path(REF, "gencode.v47.primary_assembly.annotation.gtf")
if (!nzchar(TXFA)) TXFA <- file.path(REF, "gencode.v47.transcripts.fa")
if (!file.exists(GTF)  && file.exists(paste0(GTF, ".gz")))  GTF  <- paste0(GTF, ".gz")
if (!file.exists(TXFA) && file.exists(paste0(TXFA, ".gz"))) TXFA <- paste0(TXFA, ".gz")
stopifnot(file.exists(GTF), file.exists(TXFA))
cat("GTF :", GTF, "\nTXFA:", TXFA, "\n")
OUT     <- file.path(RESULTS, "16_isoformswitch")
FIG     <- file.path(OUT, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
theme_set(theme_bw(base_size = 12))

# STEP 09가 확정한 샘플 구성(배치 포함)을 따릅니다.
.cd <- file.path(RESULTS, "09_deseq2/coldata_used.tsv")
meta <- if (file.exists(.cd)) fread(.cd) else fread(file.path(PROJ, "meta/samples.tsv"))
SALMON <- file.path(RESULTS, "08_salmon")

# [deposit] the deposit omits salmon's Gibbs samples (aux_info/bootstrap/bootstraps.gz, ~380 MB). tximport would
#           still try to read them (meta_info.json: num_bootstraps=30) although ISA never uses inferential
#           replicates (only counts/abundance/length). If they are missing, read from a temp mirror whose
#           meta_info.json says num_bootstraps=0; quant.sf is the same file, so results are unchanged.
.sdirs <- list.dirs(SALMON, recursive = FALSE)
.sdirs <- .sdirs[file.exists(file.path(.sdirs, "quant.sf"))]
.noboot <- vapply(.sdirs, function(d) {
    mi <- file.path(d, "aux_info", "meta_info.json")
    file.exists(mi) && any(grepl('"num_bootstraps": *[1-9]', readLines(mi))) &&
        !file.exists(file.path(d, "aux_info", "bootstrap", "bootstraps.gz"))
}, logical(1))
if (any(.noboot)) {
    .mirror <- file.path(tempdir(), "salmon_noboot")
    for (d in .sdirs) {
        m <- file.path(.mirror, basename(d)); dir.create(file.path(m, "aux_info"), recursive = TRUE)
        file.symlink(normalizePath(file.path(d, "quant.sf")), file.path(m, "quant.sf"))
        for (f in c("cmd_info.json", "lib_format_counts.json"))
            if (file.exists(file.path(d, f))) file.copy(file.path(d, f), m)
        writeLines(sub('"num_bootstraps": *[0-9]+', '"num_bootstraps": 0',
                       readLines(file.path(d, "aux_info", "meta_info.json"))),
                   file.path(m, "aux_info", "meta_info.json"))
    }
    cat("[deposit] salmon bootstraps absent -> reading quant.sf without inferential replicates\n")
    SALMON <- .mirror
}

# ── Salmon 결과 읽기 ──────────────────────────────────────────
quant <- importIsoformExpression(
    parentDir          = SALMON,
    addIsofomIdAsColumn = TRUE,   # 첫 컬럼에 isoform_id (ISA 4.x 인자명)
    showProgress       = FALSE
)
cat(sprintf("isoform %d개 x 샘플 %d개\n", nrow(quant$abundance), ncol(quant$abundance) - 1))

design <- data.frame(sampleID = colnames(quant$abundance)[-1])
design$condition <- meta$group[match(design$sampleID, meta$sample)]
# 배치를 공변량으로 넣습니다. STEP 09에서 확인한 대로 이 데이터는 배치가
# 분산의 75.6%를 차지하므로, 빼놓으면 DTU도 배치 신호에 묻힙니다.
if ("batch" %in% names(meta)) {
    design$batch <- meta$batch[match(design$sampleID, meta$sample)]
    cat("배치 공변량 포함\n")
}
stopifnot(!any(is.na(design$condition)))
print(design)

sw <- importRdata(
    isoformCountMatrix   = quant$counts,
    isoformRepExpression = quant$abundance,
    designMatrix         = design,
    isoformExonAnnoation = GTF,
    isoformNtFasta       = TXFA,
    comparisonsToMake    = data.frame(condition_1 = "Control", condition_2 = "Progerin"),
    showProgress         = FALSE,
    ignoreAfterPeriod    = FALSE,
    fixStringTieAnnotationProblem = FALSE
)

# ── 필터링 + DTU 검정 ─────────────────────────────────────────
# preFilter는 showProgress 인자가 없습니다 (ISA 4.x)
sw <- preFilter(sw, geneExpressionCutoff = 1, isoformExpressionCutoff = 0,
                removeSingleIsoformGenes = TRUE, quiet = TRUE)
cat(sprintf("필터 후 isoform: %d\n", nrow(sw$isoformFeatures)))

sw <- isoformSwitchTestDEXSeq(sw, reduceToSwitchingGenes = FALSE, showProgress = FALSE)

feat <- as.data.table(sw$isoformFeatures)
fwrite(feat, file.path(OUT, "isoform_features.tsv"), sep = "\t")

sig <- feat[isoform_switch_q_value < 0.05 & abs(dIF) > 0.1]
cat(sprintf("유의 DTU isoform: %d개 (유전자 %d개)\n", nrow(sig), uniqueN(sig$gene_name)))
fwrite(sig[order(isoform_switch_q_value)], file.path(OUT, "significant_DTU.tsv"), sep = "\t")

# ── 01 dIF volcano ────────────────────────────────────────────
v <- feat[!is.na(isoform_switch_q_value) & !is.na(dIF)]
v[, signif := isoform_switch_q_value < 0.05 & abs(dIF) > 0.1]
png(file.path(FIG, "01_dIF_volcano.png"), 7, 6, units = "in", res = 300, type = "cairo")
print(ggplot(v, aes(dIF, -log10(pmax(isoform_switch_q_value, 1e-300)), color = signif)) +
      geom_point(size = .8, alpha = .6) +
      scale_color_manual(values = c("FALSE" = "grey75", "TRUE" = "#D55E00")) +
      geom_vline(xintercept = c(-.1, .1), linetype = 2) +
      geom_hline(yintercept = -log10(0.05), linetype = 2) +
      labs(x = "dIF (isoform fraction 변화, Progerin − Control)",
           y = "−log10 q", color = "유의", title = "차등 isoform 사용 (DTU)"))
dev.off()

# ── 02 유전자 발현 변화 vs isoform 사용 변화 ─────────────────
png(file.path(FIG, "02_dIF_vs_logFC.png"), 7, 6, units = "in", res = 300, type = "cairo")
print(ggplot(v, aes(gene_log2_fold_change, dIF, color = signif)) +
      geom_point(size = .8, alpha = .6) +
      scale_color_manual(values = c("FALSE" = "grey75", "TRUE" = "#D55E00")) +
      geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 0, linetype = 2) +
      labs(x = "유전자 log2FC", y = "dIF", color = "유의",
           title = "발현 변화와 isoform 사용 변화는 독립적",
           subtitle = "세로축만 큰 점 = 총 발현은 그대로인데 isoform 비율만 바뀐 유전자"))
dev.off()

# ── 03 기능적 귀결 (서열 예측 도구 없이 가능한 항목만) ───────
sw2 <- tryCatch({
    s <- analyzeORF(sw, showProgress = FALSE)
    s <- analyzeAlternativeSplicing(s, showProgress = FALSE)
    analyzeSwitchConsequences(
        s,
        # ORF_seq_similarity는 extractSequence()로 아미노산 서열을 먼저
        # 넣어야 해서 제외합니다. 나머지는 좌표 기반이라 그대로 계산됩니다.
        consequencesToAnalyze = c("intron_retention", "NMD_status",
                                  "tss", "tts", "last_exon"),
        dIFcutoff = 0.1, showProgress = FALSE)
}, error = function(e) { message("consequence 분석 실패: ", conditionMessage(e)); NULL })

if (!is.null(sw2)) {
    saveRDS(sw2, file.path(OUT, "switchList.rds"))

    cons <- tryCatch(extractConsequenceSummary(sw2, asFractionTotal = FALSE, plot = FALSE),
                     error = function(e) NULL)
    if (!is.null(cons)) {
        fwrite(as.data.table(cons), file.path(OUT, "consequence_summary.tsv"), sep = "\t")
        png(file.path(FIG, "03_consequences.png"), 8, 6, units = "in", res = 300, type = "cairo")
        print(extractConsequenceSummary(sw2, asFractionTotal = FALSE, plot = TRUE,
                                        returnResult = FALSE))
        dev.off()
    }

    # splicing 유형 요약
    sp <- tryCatch(extractSplicingSummary(sw2, asFractionTotal = FALSE, plot = FALSE),
                   error = function(e) NULL)
    if (!is.null(sp)) {
        fwrite(as.data.table(sp), file.path(OUT, "splicing_summary.tsv"), sep = "\t")
        png(file.path(FIG, "04_splicing_summary.png"), 8, 6, units = "in", res = 300, type = "cairo")
        print(extractSplicingSummary(sw2, asFractionTotal = FALSE, plot = TRUE,
                                     returnResult = FALSE))
        dev.off()
    }

    # ── 05 상위 switch 유전자 개별 plot ──────────────────────
    top <- head(sig[order(isoform_switch_q_value)], 12)
    genes <- unique(top$gene_name)
    n_ok <- 0
    png(file.path(FIG, "05_switchplots.png"), 11, 7, units = "in", res = 300, type = "cairo")
    for (g in genes) {
        ok <- tryCatch({
            print(switchPlot(sw2, gene = g,
                             condition1 = "Control", condition2 = "Progerin"))
            TRUE
        }, error = function(e) { message("  switchPlot 실패 (", g, "): ", conditionMessage(e)); FALSE })
        if (ok) n_ok <- n_ok + 1
    }
    dev.off()
    cat(sprintf("switch plot %d개 생성\n", n_ok))

    # LMNA 개별
    if ("LMNA" %in% feat$gene_name) {
        png(file.path(FIG, "06_LMNA_switch.png"), 11, 7, units = "in", res = 300, type = "cairo")
        tryCatch(print(switchPlot(sw2, gene = "LMNA",
                                  condition1 = "Control", condition2 = "Progerin")),
                 error = function(e) message("  LMNA switchPlot 실패: ", conditionMessage(e)))
        dev.off()
    }
} else {
    saveRDS(sw, file.path(OUT, "switchList.rds"))
}

cat("\n완료 →", FIG, "\n")
