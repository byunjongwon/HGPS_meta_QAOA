#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────
# STEP 09 — DESeq2 차등발현 분석 + 그림 일괄 생성
#   비교: Progerin vs Control (n=4 vs 4)
#
# 산출 그림 (figures/):
#   01_library_size / 02_counts_boxplot / 03_pca / 03b_pca_scree
#   04_sample_distance / 05_sample_correlation / 06_dispersion
#   07_pvalue_hist / 08_ma / 09_volcano / 10_heatmap_top50
#   11_topgenes / 12_marker_panels (VSMC/SASP/senescence/cellcycle)
#   13_lmna_locus / 14_deg_barplot
#
# 산출 표:
#   DE_Progerin_vs_Control.tsv, normalized_counts.tsv, vst_counts.tsv,
#   tpm_matrix.tsv, dds.rds, vsd.rds
# ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
    library(DESeq2); library(data.table); library(ggplot2); library(ggrepel)
    library(pheatmap); library(RColorBrewer); library(EnhancedVolcano)
    library(org.Hs.eg.db); library(AnnotationDbi); library(matrixStats)
})
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])),
                 "_plot_setup.R"))

# [deposit] PROJ default = dataset dir (parent of this scripts/ dir; was ~/progerin_rnaseq); overridable by PROJ
PROJ    <- Sys.getenv("PROJ", normalizePath(file.path(dirname(sub("--file=", "",
               grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "..")))
RESULTS <- file.path(PROJ, "results")
OUT     <- file.path(RESULTS, "09_deseq2")
FIG     <- file.path(OUT, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_bw(base_size = 12))
GRP_COL <- c(Control = "#0072B2", Progerin = "#D55E00")
FDR <- 0.05; LFC <- 1

# ── 입력 ──────────────────────────────────────────────────────
meta <- fread(file.path(PROJ, "meta/samples.tsv"))
cts  <- fread(file.path(RESULTS, "07_counts/counts_matrix.tsv"))

cts[, gene_id := sub("\\..*$", "", gene_id)]
cts <- cts[, lapply(.SD, sum), by = gene_id]        # PAR_Y 중복 합치기
mat <- as.matrix(cts[, -1]); rownames(mat) <- cts$gene_id
mat <- mat[, meta$sample]
storage.mode(mat) <- "integer"

coldata <- data.frame(group = factor(meta$group, levels = c("Control", "Progerin")),
                      replicate = factor(meta$replicate),
                      row.names = meta$sample)
stopifnot(identical(rownames(coldata), colnames(mat)))
cat(sprintf("입력: %d 유전자 x %d 샘플\n", nrow(mat), ncol(mat)))

# ─────────────────────────────────────────────────────────────
# 반복 구조 진단과 설계 결정 (2026-08-08, 데이터 기반)
#
# GEO 메타데이터에는 짝 정보가 없습니다(8개 모두 별도 BioSample,
# "four independent replicates per condition"으로 신고). 실제 구조는
# 상관·클러스터링·PCA로 확인했습니다.
# ─────────────────────────────────────────────────────────────
# 1차 진단에서 replicate (1,4)와 (2,3)이 한 덩어리로 묶였습니다.
#   Control : 1↔4 r=0.9956, 2↔3 r=0.9947  (짝 밖 0.956~0.958)
#   Progerin: 1↔4 r=0.9824, 2↔3 r=0.9936  (짝 밖 0.907~0.930)
#
# PCA를 보니 이 구조의 정체가 드러났습니다:
#   PC1(75.6%)이 그룹이 아니라 짝으로 갈립니다 — (1,4) 전부 음수,
#   (2,3) 전부 양수. 그룹을 나누는 건 PC2(19.4%)입니다.
#   → 같은 세포의 반복 측정이 아니라 "두 번의 독립 실험 배치"이고,
#     각 배치에 Control 2 + Progerin 2가 들어간 완전 균형 2x2 설계입니다.
#
# 배치가 그룹과 교락되지 않고 균형을 이루므로, 올바른 처리는
# 합치는 것(collapse)이 아니라 모델에 배치를 넣는 것(blocking)입니다.
#   ~ group          n=8 : DEG 540  (배치 분산이 오차로 들어가 검정력 손실)
#   ~ batch + group  n=8 : DEG 893  ← 채택
#   collapse         n=2v2: DEG  65  (반복을 버려 과소검출)
# 배치 효과 자체는 2,244개 유전자로 그룹 효과보다 큽니다.
PAIRS <- list(A = c(1, 4), B = c(2, 3))   # replicate 번호 기준 배치

# DESIGN: "blocked"(기본, ~batch+group) | "collapse"(n=2v2) | "simple"(~group)
DESIGN   <- Sys.getenv("DESIGN", "blocked")
COLLAPSE <- DESIGN == "collapse"

# ── 진단: 짝 내부 상관 vs 짝 외부 상관 ───────────────────────
diag_dds <- DESeqDataSetFromMatrix(mat, coldata, design = ~ group)
diag_dds <- diag_dds[rowSums(counts(diag_dds) >= 10) >= 4, ]
diag_vsd <- vst(diag_dds, blind = TRUE)
cm <- cor(assay(diag_vsd), method = "spearman")

pair_of <- setNames(rep(NA_character_, nrow(meta)), meta$sample)
for (p in names(PAIRS)) pair_of[meta$sample[meta$replicate %in% PAIRS[[p]]]] <- p

cmp <- data.table()
for (i in 1:(ncol(cm) - 1)) for (j in (i + 1):ncol(cm)) {
    si <- colnames(cm)[i]; sj <- colnames(cm)[j]
    if (coldata[si, "group"] != coldata[sj, "group"]) next   # 그룹 내 쌍만
    cmp <- rbind(cmp, data.table(
        s1 = si, s2 = sj, group = as.character(coldata[si, "group"]),
        relation = if (pair_of[si] == pair_of[sj]) sprintf("짝 내부 (%s)", paste(sapply(names(PAIRS), function(p) paste(PAIRS[[p]], collapse="-")), collapse=", ")) else "짝 외부",
        spearman = cm[i, j]))
}
fwrite(cmp, file.path(OUT, "replicate_structure_check.tsv"), sep = "\t")

cat("\n── 반복 구조 진단 (그룹 내 샘플 쌍의 Spearman 상관) ──\n")
print(cmp[order(-spearman)])
agg <- cmp[, .(n = .N, mean_r = mean(spearman), min_r = min(spearman),
               max_r = max(spearman)), by = relation]
print(agg)
if (nrow(agg) == 2) {
    d <- agg[relation %like% "내부", mean_r] - agg[relation %like% "외부", mean_r]
    cat(sprintf("\n짝 내부 − 짝 외부 평균 상관 차이: %+.4f\n", d))
    cat("  → 기술적 반복이라면 짝 내부가 뚜렷하게 높아야 합니다(보통 +0.01 이상, r>0.99 수준).\n")
    cat("  → 차이가 미미하거나 음수면 지정한 짝 가정은 데이터로 지지되지 않습니다.\n")
}

png(file.path(FIG, "00_replicate_structure.png"), 7, 5, units = "in", res = 300, type = "cairo")
print(ggplot(cmp, aes(relation, spearman, color = group)) +
      geom_boxplot(outlier.shape = NA, color = "grey50") +
      geom_jitter(width = .12, size = 3) +
      labs(x = NULL, y = "Spearman 상관 (VST)",
           title = "반복 구조 진단 — 그룹 내 샘플 쌍 상관",
           subtitle = "기술적 반복이면 '짝 내부'가 뚜렷하게 높아야 합니다"))
dev.off()

# ── 설계 적용 ────────────────────────────────────────────────
mat_full <- mat; coldata_full <- coldata          # 원본 8샘플 보관
coldata$batch <- factor(pair_of[rownames(coldata)])
coldata_full$batch <- coldata$batch

if (COLLAPSE) {
    grp_lab <- factor(paste0(as.character(coldata$group), "_", pair_of[rownames(coldata)]))
    tmp <- DESeqDataSetFromMatrix(mat, coldata, design = ~ group)
    tmp <- collapseReplicates(tmp, groupby = grp_lab, renameCols = TRUE)

    mat <- counts(tmp)
    storage.mode(mat) <- "integer"
    coldata <- data.frame(
        group = factor(sub("_[AB]$", "", colnames(mat)), levels = c("Control", "Progerin")),
        row.names = colnames(mat))
    cat(sprintf("\ncollapseReplicates 적용: %d 샘플 → %d 샘플 (조건당 %d개)\n",
                ncol(mat_full), ncol(mat), ncol(mat) / 2))
    print(coldata)
}

# 설계식 결정
DESIGN_FORMULA <- switch(DESIGN,
    blocked  = ~ batch + group,   # 균형 2x2 — 배치를 blocking factor로
    collapse = ~ group,           # 배치별로 합쳤으므로 batch는 이미 흡수됨
    simple   = ~ group,           # 배치 무시 (비교용)
    stop("DESIGN은 blocked | collapse | simple 중 하나여야 합니다: ", DESIGN))
cat(sprintf("\n설계: DESIGN=%s / formula = %s / 샘플 %d개\n",
            DESIGN, paste(deparse(DESIGN_FORMULA), collapse = ""), ncol(mat)))
print(table(coldata$group, if ("batch" %in% names(coldata)) coldata$batch else ""))

# 이후 단계(11/12/17)가 동일한 샘플 구성을 쓰도록 기록
fwrite(data.table(sample = rownames(coldata), group = as.character(coldata$group),
                  batch = if ("batch" %in% names(coldata)) as.character(coldata$batch) else NA),
       file.path(OUT, "coldata_used.tsv"), sep = "\t")

# ── TPM 행렬 (길이 보정, 시각화/스코어링용) ───────────────────
glen <- fread(file.path(RESULTS, "07_counts/gene_lengths.tsv"),
              header = FALSE, col.names = c("gene_id", "length"))
glen[, gene_id := sub("\\..*$", "", gene_id)]
glen <- glen[, .(length = max(length)), by = gene_id]
lv <- setNames(glen$length, glen$gene_id)[rownames(mat)]
rpk <- mat / (lv / 1000)
tpm <- t(t(rpk) / colSums(rpk, na.rm = TRUE) * 1e6)
fwrite(data.table(gene_id = rownames(tpm), as.data.frame(tpm)),
       file.path(OUT, "tpm_matrix.tsv"), sep = "\t")

# ── DESeq2 ────────────────────────────────────────────────────
MIN_GRP <- min(table(coldata$group))                # collapse 여부에 따라 2 또는 4
dds  <- DESeqDataSetFromMatrix(mat, coldata, design = DESIGN_FORMULA)
keep <- rowSums(counts(dds) >= 10) >= MIN_GRP       # 최소 그룹 크기
dds  <- dds[keep, ]
cat(sprintf("필터 후: %d 유전자\n", nrow(dds)))

dds <- DESeq(dds)
vsd <- vst(dds, blind = TRUE)

fwrite(data.table(gene_id = rownames(counts(dds, normalized = TRUE)),
                  as.data.frame(counts(dds, normalized = TRUE))),
       file.path(OUT, "normalized_counts.tsv"), sep = "\t")
fwrite(data.table(gene_id = rownames(assay(vsd)), as.data.frame(assay(vsd))),
       file.path(OUT, "vst_counts.tsv"), sep = "\t")

sym <- mapIds(org.Hs.eg.db, keys = rownames(dds), keytype = "ENSEMBL",
              column = "SYMBOL", multiVals = "first")
ann <- data.frame(group = coldata$group, row.names = rownames(coldata))

# ── 01 라이브러리 크기 ────────────────────────────────────────
d <- data.frame(sample = colnames(mat), total = colSums(mat) / 1e6, group = coldata$group)
png(file.path(FIG, "01_library_size.png"), 7, 5, units = "in", res = 300, type = "cairo")
print(ggplot(d, aes(reorder(sample, total), total, fill = group)) +
      geom_col() + coord_flip() + scale_fill_manual(values = GRP_COL) +
      labs(x = NULL, y = "Assigned read pairs (M)", title = "Library size"))
dev.off()

# ── 02 정규화 전/후 분포 ──────────────────────────────────────
mk <- function(m, stage) {
    dt <- as.data.table(as.table(m)); setnames(dt, c("gene", "sample", "value"))
    dt[, stage := stage][]
}
box <- rbind(mk(log2(mat + 1), "Raw (log2)"), mk(assay(vsd), "VST"))
box[, group := coldata[as.character(sample), "group"]]
png(file.path(FIG, "02_counts_boxplot.png"), 10, 5, units = "in", res = 300, type = "cairo")
print(ggplot(box, aes(sample, value, fill = group)) +
      geom_boxplot(outlier.size = 0.2) + facet_wrap(~stage, scales = "free_y") +
      scale_fill_manual(values = GRP_COL) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
      labs(x = NULL, y = "Expression", title = "Normalization effect"))
dev.off()

# ── 03 PCA + scree ────────────────────────────────────────────
pv  <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 500)
pca <- prcomp(t(assay(vsd)[pv, ]))
pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pd  <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                  group = coldata$group, sample = rownames(coldata),
                  batch = if ("batch" %in% names(coldata)) coldata$batch else factor("-"))
png(file.path(FIG, "03_pca.png"), 7.5, 6, units = "in", res = 300, type = "cairo")
print(ggplot(pd, aes(PC1, PC2, color = group, shape = batch, label = sample)) +
      geom_point(size = 4) + geom_text_repel(size = 2.8, show.legend = FALSE) +
      scale_color_manual(values = GRP_COL) + scale_shape_manual(values = c(16, 17, 15)) +
      labs(x = sprintf("PC1 (%.1f%%)", pct[1]), y = sprintf("PC2 (%.1f%%)", pct[2]),
           title = "PCA (top 500 variable genes, VST)",
           subtitle = "색=그룹, 모양=배치"))
dev.off()

# ── 03c 배치 구조 진단 PCA (원본 8샘플, 설계 근거 그림) ──────
if (ncol(mat_full) > ncol(mat) || DESIGN != "collapse") {
    d8 <- DESeqDataSetFromMatrix(mat_full, coldata_full, design = ~ group)
    d8 <- d8[rowSums(counts(d8) >= 10) >= 4, ]
    v8 <- assay(vst(d8, blind = TRUE))
    pv8 <- head(order(rowVars(v8), decreasing = TRUE), 500)
    p8 <- prcomp(t(v8[pv8, ]))
    pct8 <- round(100 * p8$sdev^2 / sum(p8$sdev^2), 1)
    pd8 <- data.frame(PC1 = p8$x[, 1], PC2 = p8$x[, 2], group = coldata_full$group,
                      batch = coldata_full$batch, sample = rownames(coldata_full))
    fwrite(as.data.table(pd8), file.path(OUT, "pca_8samples.tsv"), sep = "\t")
    png(file.path(FIG, "03c_pca_batch_structure.png"), 7.5, 6, units = "in", res = 300, type = "cairo")
    print(ggplot(pd8, aes(PC1, PC2, color = group, shape = batch, label = sample)) +
          geom_point(size = 4.5) + geom_text_repel(size = 2.8, show.legend = FALSE) +
          scale_color_manual(values = GRP_COL) + scale_shape_manual(values = c(16, 17)) +
          labs(x = sprintf("PC1 (%.1f%%)", pct8[1]), y = sprintf("PC2 (%.1f%%)", pct8[2]),
               title = "설계 근거 — 원본 8샘플 PCA",
               subtitle = sprintf("PC1(%.1f%%)은 그룹이 아니라 배치로 갈립니다. 그룹은 PC2(%.1f%%).",
                                  pct8[1], pct8[2])))
    dev.off()
    cat(sprintf("\n배치 구조 PCA: PC1 %.1f%% / PC2 %.1f%%\n", pct8[1], pct8[2]))
}

png(file.path(FIG, "03b_pca_scree.png"), 6, 4, units = "in", res = 300, type = "cairo")
print(ggplot(data.frame(PC = factor(paste0("PC", seq_along(pct)), levels = paste0("PC", seq_along(pct))),
                        var = pct), aes(PC, var)) +
      geom_col(fill = "grey40") + labs(y = "Variance explained (%)", x = NULL, title = "Scree plot"))
dev.off()

# ── 04 거리 / 05 상관 ─────────────────────────────────────────
png(file.path(FIG, "04_sample_distance.png"), 7, 6, units = "in", res = 300, type = "cairo")
pheatmap(as.matrix(dist(t(assay(vsd)))), annotation_col = ann,
         annotation_colors = list(group = GRP_COL),
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample-to-sample distance (VST)")
dev.off()

png(file.path(FIG, "05_sample_correlation.png"), 7, 6, units = "in", res = 300, type = "cairo")
pheatmap(cor(assay(vsd), method = "spearman"), annotation_col = ann,
         annotation_colors = list(group = GRP_COL),
         col = colorRampPalette(brewer.pal(9, "YlOrRd"))(255),
         display_numbers = TRUE, number_format = "%.3f",
         main = "Spearman correlation (VST)")
dev.off()

# ── 06 dispersion ─────────────────────────────────────────────
png(file.path(FIG, "06_dispersion.png"), 6, 5, units = "in", res = 300, type = "cairo")
plotDispEsts(dds, main = "Dispersion estimates")
dev.off()

# ── 차등발현 ──────────────────────────────────────────────────
nm  <- "Progerin_vs_Control"
res <- results(dds, contrast = c("group", "Progerin", "Control"), alpha = FDR)
# 2군 설계이므로 coef 기반 apeglm 사용 가능 (ashr보다 권장)
res_s <- lfcShrink(dds, coef = "group_Progerin_vs_Control", type = "apeglm")

dt <- data.table(gene_id = rownames(res_s), symbol = sym[rownames(res_s)],
                 baseMean = res_s$baseMean, log2FoldChange = res_s$log2FoldChange,
                 lfcSE = res_s$lfcSE, pvalue = res$pvalue, padj = res$padj)
setorder(dt, padj, na.last = TRUE)
fwrite(dt, file.path(OUT, sprintf("DE_%s.tsv", nm)), sep = "\t")

sig <- dt[!is.na(padj) & padj < FDR & abs(log2FoldChange) > LFC]
fwrite(sig, file.path(OUT, sprintf("DEG_sig_%s.tsv", nm)), sep = "\t")
cat(sprintf("\n%s: DEG %d개 (up %d / down %d)  [FDR<%g, |log2FC|>%g]\n",
            nm, nrow(sig), sum(sig$log2FoldChange > 0), sum(sig$log2FoldChange < 0), FDR, LFC))

# ── 07 p-value 분포 ───────────────────────────────────────────
png(file.path(FIG, "07_pvalue_hist.png"), 6, 4, units = "in", res = 300, type = "cairo")
print(ggplot(dt[!is.na(pvalue)], aes(pvalue)) +
      geom_histogram(bins = 50, fill = "grey40") +
      labs(title = "p-value distribution", y = "genes",
           subtitle = "평평한 배경 + 0 근처 봉우리 = 건전한 검정"))
dev.off()

# ── 08 MA ─────────────────────────────────────────────────────
png(file.path(FIG, "08_ma.png"), 6, 5, units = "in", res = 300, type = "cairo")
plotMA(res_s, ylim = c(-8, 8), main = sprintf("MA plot — %s (apeglm)", nm))
dev.off()

# ── 09 Volcano ────────────────────────────────────────────────
lab <- ifelse(is.na(dt$symbol), dt$gene_id, dt$symbol)
png(file.path(FIG, "09_volcano.png"), 8, 7, units = "in", res = 300, type = "cairo")
print(EnhancedVolcano(dt, lab = lab, x = "log2FoldChange", y = "padj",
                      pCutoff = FDR, FCcutoff = LFC, title = nm,
                      subtitle = sprintf("FDR<%g, |log2FC|>%g  (n=%d DEG)", FDR, LFC, nrow(sig)),
                      labSize = 3, pointSize = 1.5, drawConnectors = TRUE, max.overlaps = 20))
dev.off()

# ── 10 상위 DEG 히트맵 ────────────────────────────────────────
if (nrow(sig) >= 5) {
    top <- head(sig$gene_id, 50)
    m <- t(scale(t(assay(vsd)[top, , drop = FALSE])))
    rownames(m) <- ifelse(is.na(sym[top]), top, sym[top])
    png(file.path(FIG, "10_heatmap_top50.png"), 7, 10, units = "in", res = 300, type = "cairo")
    pheatmap(m, annotation_col = ann, annotation_colors = list(group = GRP_COL),
             fontsize_row = 7, main = "Top 50 DEG (z-score, VST)")
    dev.off()
}

# ── 11 상위 유전자 개별 counts ────────────────────────────────
if (nrow(sig) >= 9) {
    top9 <- head(sig$gene_id, 9)
    pl <- rbindlist(lapply(top9, function(g)
        data.table(gene = ifelse(is.na(sym[g]), g, sym[g]),
                   count = counts(dds, normalized = TRUE)[g, ], group = coldata$group)))
    png(file.path(FIG, "11_topgenes.png"), 9, 7, units = "in", res = 300, type = "cairo")
    print(ggplot(pl, aes(group, count, color = group)) +
          geom_boxplot(outlier.shape = NA) + geom_jitter(width = .15, size = 2) +
          facet_wrap(~gene, scales = "free_y") + scale_y_log10() +
          scale_color_manual(values = GRP_COL) +
          labs(y = "Normalized count (log10)", x = NULL, title = "Top DEG"))
    dev.off()
}

# ── 12 마커 패널 (이 논문의 핵심 주장 검증) ───────────────────
PANELS <- list(
    `VSMC contractile (수축형)`   = c("MYH11", "ACTA2", "CNN1", "TAGLN", "LMOD1", "MYOCD", "SRF"),
    `VSMC synthetic (합성/증식형)`= c("KLF4", "SPP1", "VCAM1", "COL1A1", "COL3A1", "FN1", "MMP2", "MMP9"),
    `Senescence / cell cycle 정지`= c("CDKN1A", "CDKN2A", "CDKN2B", "TP53", "RB1", "LMNB1", "GLB1"),
    `SASP (염증성 분비)`          = c("IL6", "IL1A", "IL1B", "CXCL8", "CXCL1", "CCL2", "TNF", "IGFBP3", "SERPINE1"),
    `Proliferation`               = c("MKI67", "PCNA", "CCNA2", "CCNB1", "CCNE1", "CDK1", "TOP2A", "E2F1"),
    `Replication stress / DDR`    = c("ATR", "ATM", "CHEK1", "CHEK2", "H2AX", "RPA1", "RRM2B", "CDC25A"),
    `Fibrosis / ECM remodeling`   = c("TGFB1", "CTGF", "CCN2", "POSTN", "ELN", "LOX", "TIMP1"),
    `Lamin / nuclear envelope`    = c("LMNA", "LMNB1", "LMNB2", "LBR", "EMD", "SUN1", "SUN2", "BANF1")
)
sym2id <- setNames(names(sym), sym)
norm_c <- counts(dds, normalized = TRUE)

panel_dt <- rbindlist(lapply(names(PANELS), function(p) {
    gs <- PANELS[[p]]; ids <- sym2id[gs]; ids <- ids[!is.na(ids)]
    if (!length(ids)) return(NULL)
    rbindlist(lapply(names(ids), function(g)
        data.table(panel = p, gene = g, sample = colnames(norm_c),
                   count = norm_c[ids[[g]], ], group = coldata$group)))
}))
fwrite(panel_dt, file.path(OUT, "marker_panel_counts.tsv"), sep = "\t")

png(file.path(FIG, "12_marker_panels.png"), 11, 7, units = "in", res = 300, type = "cairo")
for (p in unique(panel_dt$panel)) {
    sub <- panel_dt[panel == p]
    print(ggplot(sub, aes(group, count, color = group)) +
          geom_boxplot(outlier.shape = NA) + geom_jitter(width = .15, size = 2) +
          facet_wrap(~gene, scales = "free_y", ncol = 4) + scale_y_log10() +
          scale_color_manual(values = GRP_COL) +
          labs(y = "Normalized count (log10)", x = NULL, title = p))
}
dev.off()

# 마커 패널 통합 히트맵
pids <- unique(unlist(lapply(PANELS, function(gs) sym2id[gs])))
pids <- pids[!is.na(pids) & pids %in% rownames(vsd)]
if (length(pids) > 3) {
    m <- t(scale(t(assay(vsd)[pids, , drop = FALSE])))
    rownames(m) <- sym[pids]
    rowann <- data.frame(panel = rep(names(PANELS), sapply(PANELS, function(gs) {
        ids <- sym2id[gs]; sum(!is.na(ids) & ids %in% pids)
    })), row.names = NULL)
    png(file.path(FIG, "12b_marker_heatmap.png"), 7, 12, units = "in", res = 300, type = "cairo")
    pheatmap(m, annotation_col = ann, annotation_colors = list(group = GRP_COL),
             fontsize_row = 7, cluster_rows = TRUE,
             main = "표현형 마커 패널 (z-score, VST)")
    dev.off()
}

# ── 13 LMNA 유전자 수준 ───────────────────────────────────────
lmna_id <- sym2id[["LMNA"]]
if (!is.na(lmna_id)) {
    pl <- data.table(count = norm_c[lmna_id, ], group = coldata$group)
    png(file.path(FIG, "13_lmna_gene.png"), 5, 5, units = "in", res = 300, type = "cairo")
    print(ggplot(pl, aes(group, count, color = group)) +
          geom_boxplot(outlier.shape = NA) + geom_jitter(width = .15, size = 3) +
          scale_color_manual(values = GRP_COL) +
          labs(y = "Normalized count", x = NULL, title = "LMNA (gene level)",
               subtitle = "렌티바이러스 progerin은 cDNA라 내인성 LMNA와 구분됩니다.\nSTEP 13 rMATS의 exon 11 cryptic junction을 함께 보세요."))
    dev.off()
}

# ── 14 DEG 개수 요약 ──────────────────────────────────────────
bd <- data.table(direction = c("Up in Progerin", "Down in Progerin"),
                 n = c(sum(sig$log2FoldChange > 0), sum(sig$log2FoldChange < 0)))
png(file.path(FIG, "14_deg_barplot.png"), 5, 4, units = "in", res = 300, type = "cairo")
print(ggplot(bd, aes(direction, n, fill = direction)) +
      geom_col(width = .6) + geom_text(aes(label = n), vjust = -0.4) +
      scale_fill_manual(values = c("Up in Progerin" = "#D55E00", "Down in Progerin" = "#0072B2")) +
      labs(x = NULL, y = "DEG 수", title = sprintf("DEG (FDR<%g, |log2FC|>%g)", FDR, LFC)) +
      theme(legend.position = "none"))
dev.off()

saveRDS(dds, file.path(OUT, "dds.rds"))
saveRDS(vsd, file.path(OUT, "vsd.rds"))

# ─────────────────────────────────────────────────────────────
# 민감도 분석 — collapse 하지 않은 원본 설계로도 돌려 비교
#   목적: 4개 데이터셋 교집합을 낼 때 이 데이터셋의 DEG 수가
#         반복 구조 가정에 얼마나 좌우되는지 정량화합니다.
# ─────────────────────────────────────────────────────────────
# 세 설계를 모두 돌려 DEG 수와 효과크기 일치도를 비교합니다.
#   blocked  ~ batch + group  (n=8)  — 균형 2x2의 표준 처리
#   simple   ~ group          (n=8)  — 배치 무시
#   collapse ~ group          (n=2v2)— 배치별로 합쳐 반복을 버림
# 4개 데이터셋 교집합에서 이 데이터의 기여도가 설계에 얼마나
# 좌우되는지 정량화하는 것이 목적입니다.
# ─────────────────────────────────────────────────────────────
fit_design <- function(tag) {
    if (tag == "collapse") {
        gl <- factor(paste0(as.character(coldata_full$group), "_",
                            pair_of[rownames(coldata_full)]))
        t2 <- collapseReplicates(
            DESeqDataSetFromMatrix(mat_full, coldata_full, design = ~ group),
            groupby = gl, renameCols = TRUE)
        am <- counts(t2); storage.mode(am) <- "integer"
        ac <- data.frame(group = factor(sub("_[AB]$", "", colnames(am)),
                                        levels = c("Control", "Progerin")),
                         row.names = colnames(am))
        fm <- ~ group
    } else {
        am <- mat_full; ac <- coldata_full
        fm <- if (tag == "blocked") ~ batch + group else ~ group
    }
    d <- DESeqDataSetFromMatrix(am, ac, design = fm)
    d <- d[rowSums(counts(d) >= 10) >= min(table(ac$group)), ]
    d <- DESeq(d, quiet = TRUE)
    r <- results(d, contrast = c("group", "Progerin", "Control"), alpha = FDR)
    rs <- lfcShrink(d, coef = "group_Progerin_vs_Control", type = "apeglm", quiet = TRUE)
    dt <- data.table(gene_id = rownames(rs), symbol = sym[rownames(rs)],
                     baseMean = rs$baseMean, log2FoldChange = rs$log2FoldChange,
                     pvalue = r$pvalue, padj = r$padj)
    setorder(dt, padj, na.last = TRUE)
    fwrite(dt, file.path(OUT, sprintf("DE_%s_%s.tsv", nm, tag)), sep = "\t")
    dt
}

DESIGNS <- c("blocked", "simple", "collapse")
res_all <- list()
for (tg in DESIGNS) {
    res_all[[tg]] <- tryCatch(fit_design(tg),
        error = function(e) { message("설계 ", tg, " 실패: ", conditionMessage(e)); NULL })
}
res_all <- res_all[!sapply(res_all, is.null)]

lbl <- c(blocked = "~ batch + group (n=8)", simple = "~ group (n=8)",
         collapse = "collapse (n=2 vs 2)")
comp <- rbindlist(lapply(names(res_all), function(tg) {
    s <- res_all[[tg]][!is.na(padj) & padj < FDR & abs(log2FoldChange) > LFC]
    data.table(design = lbl[[tg]], tag = tg, n_DEG = nrow(s),
               up = sum(s$log2FoldChange > 0), down = sum(s$log2FoldChange < 0))
}))
comp[, primary := tag == DESIGN]
fwrite(comp, file.path(OUT, "design_sensitivity.tsv"), sep = "\t")
cat("\n── 설계별 DEG 비교 ──\n"); print(comp)

sets <- lapply(res_all, function(d)
    d[!is.na(padj) & padj < FDR & abs(log2FoldChange) > LFC, gene_id])
if (length(sets) >= 2) {
    core <- Reduce(intersect, sets)
    fwrite(data.table(gene_id = core, symbol = sym[core]),
           file.path(OUT, "DEG_core_all_designs.tsv"), sep = "\t")
    cat(sprintf("세 설계 공통 DEG: %d개 (교집합용 핵심 목록)\n", length(core)))
}

png(file.path(FIG, "15_design_sensitivity.png"), 7, 4.5, units = "in", res = 300, type = "cairo")
print(ggplot(comp, aes(reorder(design, -n_DEG), n_DEG, fill = primary)) +
      geom_col(width = .6) + geom_text(aes(label = n_DEG), vjust = -0.3) +
      scale_fill_manual(values = c("TRUE" = "#D55E00", "FALSE" = "grey70"),
                        labels = c("TRUE" = "채택", "FALSE" = "비교용"), name = NULL) +
      labs(x = NULL, y = "DEG 수",
           title = sprintf("설계별 DEG 수 (FDR<%g, |log2FC|>%g)", FDR, LFC),
           subtitle = "배치가 그룹과 균형을 이루므로 blocking이 표준 처리입니다") +
      theme(axis.text.x = element_text(angle = 12, hjust = 1)))
dev.off()

if (length(res_all) >= 2) {
    ref <- res_all[[DESIGN]]
    for (tg in setdiff(names(res_all), DESIGN)) {
        mg <- merge(ref[, .(gene_id, a = log2FoldChange)],
                    res_all[[tg]][, .(gene_id, b = log2FoldChange)], by = "gene_id")
        rr <- cor(mg$a, mg$b, use = "complete.obs")
        cat(sprintf("log2FC 상관 %s vs %s: r = %.4f\n", DESIGN, tg, rr))
        png(file.path(FIG, sprintf("16_lfc_%s_vs_%s.png", DESIGN, tg)), 6, 6, units = "in", res = 300, type = "cairo")
        print(ggplot(mg, aes(a, b)) + geom_point(size = .5, alpha = .3) +
              geom_abline(slope = 1, intercept = 0, color = "#D55E00", linetype = 2) +
              labs(x = sprintf("log2FC (%s)", lbl[[DESIGN]]), y = sprintf("log2FC (%s)", lbl[[tg]]),
                   title = "설계 간 효과크기 일치도",
                   subtitle = sprintf("Pearson r = %.4f", rr)))
        dev.off()
    }
}

cat("\n완료 →", FIG, "\n")
