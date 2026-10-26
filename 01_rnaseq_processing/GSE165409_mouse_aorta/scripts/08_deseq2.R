#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────
# STEP 8. DESeq2 — 차등발현 분석 + 그림 일괄 생성
#
# 산출 그림:
#   01_library_size.pdf        라이브러리 크기 bar
#   02_counts_boxplot.pdf      정규화 전/후 분포
#   03_pca.pdf                 PCA (그룹별 색)
#   04_sample_distance.pdf     샘플 거리 히트맵
#   05_sample_correlation.pdf  상관 히트맵
#   06_dispersion.pdf          dispersion estimate
#   07_pvalue_hist.pdf         p-value 분포 (검정 건전성 확인)
#   08_ma_{비교}.pdf           MA plot (apeglm shrinkage)
#   09_volcano_{비교}.pdf      Volcano
#   10_heatmap_{비교}.pdf      상위 DEG 히트맵
#   11_topgenes_{비교}.pdf     상위 유전자 개별 counts
#   12_upset_DEG.pdf           세 비교 간 DEG 중첩
#   13_lmna_gene.pdf           Lmna 유전자 발현 (progerin 맥락)
#
# 산출 표:
#   DE_{비교}.tsv, normalized_counts.tsv, vst_counts.tsv
# ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
    library(DESeq2); library(data.table); library(ggplot2); library(ggrepel)
    library(pheatmap); library(RColorBrewer); library(ComplexHeatmap)
    library(EnhancedVolcano); library(UpSetR); library(org.Mm.eg.db); library(AnnotationDbi)
})

# UTF-8(Δ, 한글) 렌더링을 위한 cairo PDF 장치 설정
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)[1])), "_plot_setup.R"))

# [deposit] default PROJ = this dataset dir (scripts/..), still overridable by PROJ (was ~/lmna_rnaseq)
.this <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
PROJ    <- Sys.getenv("PROJ", normalizePath(file.path(dirname(.this), "..")))
RESULTS <- file.path(PROJ, "results")
OUT     <- file.path(RESULTS, "08_deseq2")
FIG     <- file.path(OUT, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_bw(base_size = 12))
GRP_COL <- c(LMNA = "#D55E00", WTyng = "#0072B2", WTold = "#009E73")

# ── 입력 ──────────────────────────────────────────────────────
meta <- fread(file.path(PROJ, "meta/samples.tsv"))
cts  <- fread(file.path(RESULTS, "05_counts/counts_matrix.tsv"))

cts[, gene_id := sub("\\..*$", "", gene_id)]
cts <- cts[, lapply(.SD, sum), by = gene_id]
mat <- as.matrix(cts[, -1]); rownames(mat) <- cts$gene_id
mat <- mat[, meta$sample]                      # 순서 정렬
storage.mode(mat) <- "integer"

coldata <- data.frame(
    group = factor(meta$group, levels = c("WTyng", "WTold", "LMNA")),
    row.names = meta$sample
)
stopifnot(identical(rownames(coldata), colnames(mat)))

cat(sprintf("입력: %d 유전자 x %d 샘플\n", nrow(mat), ncol(mat)))

# ── DESeq2 ────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(mat, coldata, design = ~ group)
keep <- rowSums(counts(dds) >= 10) >= 6        # 최소 그룹 크기
dds <- dds[keep, ]
cat(sprintf("필터 후: %d 유전자\n", nrow(dds)))

dds <- DESeq(dds)
# [deposit] optional TINY_TEST=1 (test subset only): vst() needs >= 1000 genes (nsub), the test subset
#   has ~50, so use varianceStabilizingTransformation() there. Default (unset) = original vst().
vsd <- if (Sys.getenv("TINY_TEST") == "1" && nrow(dds) < 1000) {
    varianceStabilizingTransformation(dds, blind = TRUE)
} else vst(dds, blind = TRUE)

fwrite(data.table(gene_id = rownames(counts(dds, normalized = TRUE)),
                  as.data.frame(counts(dds, normalized = TRUE))),
       file.path(OUT, "normalized_counts.tsv"), sep = "\t")
fwrite(data.table(gene_id = rownames(assay(vsd)), as.data.frame(assay(vsd))),
       file.path(OUT, "vst_counts.tsv"), sep = "\t")

sym <- mapIds(org.Mm.eg.db, keys = rownames(dds), keytype = "ENSEMBL",
              column = "SYMBOL", multiVals = "first")

# ── 01 라이브러리 크기 ────────────────────────────────────────
d <- data.frame(sample = colnames(mat), total = colSums(mat) / 1e6,
                group = coldata$group)
pdf(file.path(FIG, "01_library_size.pdf"), 8, 5)
print(ggplot(d, aes(reorder(sample, total), total, fill = group)) +
      geom_col() + coord_flip() + scale_fill_manual(values = GRP_COL) +
      labs(x = NULL, y = "Assigned reads (M)", title = "Library size"))
dev.off()

# ── 02 정규화 전/후 분포 ──────────────────────────────────────
raw_l  <- reshape2::melt(log2(mat + 1));           raw_l$stage  <- "Raw (log2)"
norm_l <- reshape2::melt(assay(vsd));              norm_l$stage <- "VST"
box <- rbind(raw_l, norm_l)
colnames(box)[1:3] <- c("gene", "sample", "value")
box$group <- coldata[as.character(box$sample), "group"]
pdf(file.path(FIG, "02_counts_boxplot.pdf"), 11, 5)
print(ggplot(box, aes(sample, value, fill = group)) +
      geom_boxplot(outlier.size = 0.2) + facet_wrap(~stage, scales = "free_y") +
      scale_fill_manual(values = GRP_COL) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
      labs(x = NULL, y = "Expression", title = "Normalization effect"))
dev.off()

# ── 03 PCA ────────────────────────────────────────────────────
pv   <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 500)
pca  <- prcomp(t(assay(vsd)[pv, ]))
pct  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pd   <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                   group = coldata$group, sample = rownames(coldata))
pdf(file.path(FIG, "03_pca.pdf"), 7, 6)
print(ggplot(pd, aes(PC1, PC2, color = group, label = sample)) +
      geom_point(size = 4) + geom_text_repel(size = 2.6, show.legend = FALSE) +
      scale_color_manual(values = GRP_COL) +
      labs(x = sprintf("PC1 (%.1f%%)", pct[1]), y = sprintf("PC2 (%.1f%%)", pct[2]),
           title = "PCA (top 500 variable genes, VST)"))
dev.off()

# ── 04 샘플 거리 / 05 상관 ────────────────────────────────────
sd_mat <- as.matrix(dist(t(assay(vsd))))
ann <- data.frame(group = coldata$group, row.names = rownames(coldata))
pdf(file.path(FIG, "04_sample_distance.pdf"), 8, 7)
pheatmap(sd_mat, annotation_col = ann, annotation_colors = list(group = GRP_COL),
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample-to-sample distance (VST)")
dev.off()

cor_mat <- cor(assay(vsd), method = "spearman")
pdf(file.path(FIG, "05_sample_correlation.pdf"), 8, 7)
pheatmap(cor_mat, annotation_col = ann, annotation_colors = list(group = GRP_COL),
         col = colorRampPalette(brewer.pal(9, "YlOrRd"))(255),
         main = "Spearman correlation (VST)")
dev.off()

# ── 06 dispersion ─────────────────────────────────────────────
pdf(file.path(FIG, "06_dispersion.pdf"), 6, 5)
plotDispEsts(dds, main = "Dispersion estimates")
dev.off()

# ── 비교별 분석 ───────────────────────────────────────────────
CONTRASTS <- list(
    LMNA_vs_WTyng  = c("group", "LMNA",  "WTyng"),
    WTold_vs_WTyng = c("group", "WTold", "WTyng"),
    LMNA_vs_WTold  = c("group", "LMNA",  "WTold")
)
FDR <- 0.05; LFC <- 1
deg_sets <- list(); pval_all <- list()

for (nm in names(CONTRASTS)) {
    ct <- CONTRASTS[[nm]]
    res <- results(dds, contrast = ct, alpha = FDR)
    # apeglm은 coef 기반이라 contrast에는 ashr 사용
    res_s <- lfcShrink(dds, contrast = ct, res = res, type = "ashr")

    dt <- data.table(gene_id = rownames(res_s), symbol = sym[rownames(res_s)],
                     as.data.frame(res_s))
    setorder(dt, padj, na.last = TRUE)
    fwrite(dt, file.path(OUT, sprintf("DE_%s.tsv", nm)), sep = "\t")

    sig <- dt[!is.na(padj) & padj < FDR & abs(log2FoldChange) > LFC]
    deg_sets[[nm]] <- sig$gene_id
    pval_all[[nm]] <- data.frame(pvalue = res$pvalue, contrast = nm)
    cat(sprintf("%-16s DEG %d개 (up %d / down %d)\n", nm, nrow(sig),
                sum(sig$log2FoldChange > 0), sum(sig$log2FoldChange < 0)))

    # 08 MA
    pdf(file.path(FIG, sprintf("08_ma_%s.pdf", nm)), 6, 5)
    plotMA(res_s, ylim = c(-6, 6), main = sprintf("MA plot — %s", nm))
    dev.off()

    # 09 Volcano
    lab <- ifelse(is.na(dt$symbol), dt$gene_id, dt$symbol)
    pdf(file.path(FIG, sprintf("09_volcano_%s.pdf", nm)), 8, 7)
    print(EnhancedVolcano(dt, lab = lab, x = "log2FoldChange", y = "padj",
                          pCutoff = FDR, FCcutoff = LFC,
                          title = nm, subtitle = sprintf("FDR<%g, |log2FC|>%g", FDR, LFC),
                          labSize = 3, pointSize = 1.5, drawConnectors = TRUE,
                          max.overlaps = 20))
    dev.off()

    # 10 DEG 히트맵
    if (nrow(sig) >= 5) {
        top <- head(sig$gene_id, 50)
        m <- assay(vsd)[top, , drop = FALSE]
        m <- t(scale(t(m)))
        rownames(m) <- ifelse(is.na(sym[top]), top, sym[top])
        pdf(file.path(FIG, sprintf("10_heatmap_%s.pdf", nm)), 8, 10)
        pheatmap(m, annotation_col = ann, annotation_colors = list(group = GRP_COL),
                 show_colnames = TRUE, fontsize_row = 7,
                 main = sprintf("Top 50 DEG — %s (z-score)", nm))
        dev.off()
    }

    # 11 상위 유전자 개별 plot
    if (nrow(sig) >= 6) {
        top6 <- head(sig$gene_id, 6)
        pl <- rbindlist(lapply(top6, function(g) {
            data.table(gene = ifelse(is.na(sym[g]), g, sym[g]),
                       count = counts(dds, normalized = TRUE)[g, ],
                       group = coldata$group)
        }))
        pdf(file.path(FIG, sprintf("11_topgenes_%s.pdf", nm)), 9, 6)
        print(ggplot(pl, aes(group, count, color = group)) +
              geom_boxplot(outlier.shape = NA) + geom_jitter(width = .15, size = 2) +
              facet_wrap(~gene, scales = "free_y") + scale_y_log10() +
              scale_color_manual(values = GRP_COL) +
              labs(y = "Normalized count (log10)", x = NULL, title = nm))
        dev.off()
    }
}

# ── 07 p-value 분포 ───────────────────────────────────────────
pv_dt <- rbindlist(lapply(names(pval_all), function(n) as.data.table(pval_all[[n]])))
pdf(file.path(FIG, "07_pvalue_hist.pdf"), 10, 4)
print(ggplot(pv_dt[!is.na(pvalue)], aes(pvalue)) +
      geom_histogram(bins = 50, fill = "grey40") + facet_wrap(~contrast) +
      labs(title = "p-value distribution", y = "genes"))
dev.off()

# ── 12 UpSet ──────────────────────────────────────────────────
if (all(lengths(deg_sets) > 0)) {
    pdf(file.path(FIG, "12_upset_DEG.pdf"), 8, 6)
    print(upset(fromList(deg_sets), order.by = "freq", text.scale = 1.3))
    dev.off()
}

# ── 13 Lmna 유전자 (progerin 맥락) ────────────────────────────
lmna_id <- names(sym)[which(sym == "Lmna")]
if (length(lmna_id) > 0) {
    g <- lmna_id[1]
    pl <- data.table(count = counts(dds, normalized = TRUE)[g, ], group = coldata$group)
    pdf(file.path(FIG, "13_lmna_gene.pdf"), 5, 5)
    print(ggplot(pl, aes(group, count, color = group)) +
          geom_boxplot(outlier.shape = NA) + geom_jitter(width = .15, size = 3) +
          scale_color_manual(values = GRP_COL) +
          labs(y = "Normalized count", x = NULL,
               title = "Lmna (gene level)",
               subtitle = "유전자 수준에서는 progerin이 보이지 않습니다.\nSTEP 7 rMATS의 exon 11 cryptic junction을 보세요."))
    dev.off()
}

saveRDS(dds, file.path(OUT, "dds.rds"))
saveRDS(vsd, file.path(OUT, "vsd.rds"))
cat("\n완료 → ", FIG, "\n")
