#!/usr/bin/env Rscript
# =============================================================================
# GSE306264 : HGPS 환자 섬유아세포 20라인  vs  정상 섬유아세포 10라인
#   + 부수: TetON GFP-progerin isogenic (Dox 6일 vs -Dox), 같은 세포주 3반복
#
# 핵심 설계 판단
#   - 각 세포주의 R1/R2/R3 (matrix 컬럼 ...pt1/pt2/pt3) 는 "같은 세포주의 반복"이다.
#     생물학적 단위는 '사람(세포주)' 이므로 collapseReplicates 로 합쳐
#     HGPS 20 vs WT 10 (n=30) 으로 검정한다.  <- GSE137083 / GSE198983 에서 겪은
#     pseudo-replication 을 여기서는 처음부터 차단한다.
#   - 비교용으로 합치지 않은 naive 결과(90 샘플)도 같이 계산해 차이를 보고한다.
# =============================================================================
suppressPackageStartupMessages({
  library(DESeq2); library(ggplot2); library(ggrepel); library(pheatmap)
  library(RColorBrewer); library(matrixStats)
})

# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
.this <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
DEPOSIT <- Sys.getenv("PROGERIA_DEPOSIT", normalizePath(file.path(dirname(.this), "..", "..", "..")))
BASE <- file.path(DEPOSIT, "02_per_dataset_DE", "GSE306264")  # [deposit] was ~/car_project/progeria_meta/GSE306264
RES  <- file.path(BASE, "results"); FIG <- file.path(RES, "figures")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
set.seed(1)

# ------------------------------------------------------------------ load ----
f <- file.path(BASE, "raw_data", "GSE306264_featureCounts_all_gene_name.txt.gz")
raw <- read.delim(gzfile(f), check.names = FALSE)
cat(sprintf("featureCounts 원본: %d 행 x %d 열\n", nrow(raw), ncol(raw)))

glen <- raw$Length
cts  <- as.matrix(raw[, 7:ncol(raw)])
rownames(cts) <- raw$Geneid
storage.mode(cts) <- "integer"
cat(sprintf("count matrix: %d 유전자 x %d 샘플\n", nrow(cts), ncol(cts)))

# --------------------------------------------------------------- coldata ----
# meta/coldata.tsv 는 GEO SOFT 원본에서 생성했고, 컬럼 순서 <-> GSM genotype 일치를
# 이미 검증했다 (스크립트 00 참고). 여기서는 순서만 다시 확인한다.
col <- read.delim(file.path(BASE, "meta", "coldata.tsv"), stringsAsFactors = FALSE)
rownames(col) <- col$sample
stopifnot(identical(rownames(col), colnames(cts)))
col$genotype <- col$condition

cat("\n=== 샘플 배정 ===\n")
print(table(col$arm, col$genotype))
cat("\nWT 세포주   :", paste(unique(col$line[col$genotype == "WT"]), collapse = ", "), "\n")
cat("HGPS 세포주 :", paste(unique(col$line[col$genotype == "HGPS"]), collapse = ", "), "\n")

# =============================================================================
# PART 1. 환자 비교  HGPS(20라인) vs WT(10라인)   -- 세포주 단위로 collapse
# =============================================================================
cat("\n\n############ PART 1. 환자 비교 (세포주 단위) ############\n")
pc <- col[col$arm == "patient", ]
pm <- cts[, rownames(pc)]

dds_naive <- DESeqDataSetFromMatrix(pm, pc, design = ~ genotype)
dds_naive$genotype <- relevel(factor(dds_naive$genotype), ref = "WT")

# --- collapse: 같은 세포주의 3반복을 합산 -> 30 컬럼 (HGPS 20 / WT 10)
dds <- collapseReplicates(dds_naive, groupby = factor(pc$line))
colData(dds)$genotype <- factor(colData(dds)$genotype, levels = c("WT", "HGPS"))
cat(sprintf("collapse: %d 샘플 -> %d 세포주\n", ncol(dds_naive), ncol(dds)))
print(table(colData(dds)$genotype))

keep <- rowSums(counts(dds) >= 10) >= 10   # 최소 10개 세포주에서 발현
cat(sprintf("발현 필터: %d / %d 유전자\n", sum(keep), nrow(dds)))
dds <- dds[keep, ]
dds <- DESeq(dds)

vsd  <- vst(dds, blind = TRUE); vm <- assay(vsd)

# ---- QC : PCA
pca <- prcomp(t(vm[order(rowVars(vm), decreasing = TRUE)[1:2000], ]))
pv  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
cat(sprintf("\n=== PCA 분산설명 PC1-3 : %.1f%% / %.1f%% / %.1f%% ===\n", pv[1], pv[2], pv[3]))
cat("   (GSE137083 은 PC1=98.2%% 였다. 정상 데이터면 40~70%% 대여야 한다)\n")

png(file.path(FIG, "QC_PCA_patient.png"), width = 7, height = 5.5, units = "in", res = 300, type = "cairo")
print(ggplot(data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2],
                        genotype = colData(dds)$genotype, line = colnames(dds)),
             aes(PC1, PC2, color = genotype, label = line)) +
  geom_point(size = 3.5) + geom_text_repel(size = 2.6, max.overlaps = 30, show.legend = FALSE) +
  scale_color_manual(values = c(WT = "#4C72B0", HGPS = "#C44E52")) +
  labs(x = sprintf("PC1 (%.1f%%)", pv[1]), y = sprintf("PC2 (%.1f%%)", pv[2]),
       title = "GSE306264: 20 HGPS vs 10 control fibroblast lines") +
  theme_bw(base_size = 11))
dev.off()

# ---- DE
res  <- results(dds, contrast = c("genotype", "HGPS", "WT"), alpha = 0.05)
resS <- lfcShrink(dds, coef = "genotype_HGPS_vs_WT", type = "apeglm", res = res)
df <- as.data.frame(resS); df$gene <- rownames(df)
df <- df[order(df$padj, -abs(df$log2FoldChange)), ]
FDR <- 0.05; LFC <- 1
df$sig <- ifelse(!is.na(df$padj) & df$padj < FDR & df$log2FoldChange >=  LFC, "UP",
          ifelse(!is.na(df$padj) & df$padj < FDR & df$log2FoldChange <= -LFC, "DOWN", "NS"))

cat("\n================ DE : HGPS vs WT (세포주 단위, 20 vs 10) ================\n")
cat(sprintf("검정 유전자          : %d\n", sum(!is.na(df$padj))))
cat(sprintf("padj<0.05            : %d\n", sum(df$padj < FDR, na.rm = TRUE)))
cat(sprintf("padj<0.05 & |LFC|>=1 : %d  (UP %d / DOWN %d)\n",
            sum(df$sig != "NS"), sum(df$sig == "UP"), sum(df$sig == "DOWN")))
cat(sprintf("|log2FC| > 5 인 유전자 : %d   (GSE137083 은 540개였다)\n",
            sum(abs(df$log2FoldChange) > 5, na.rm = TRUE)))
cat("\n--- 상위 30 (padj) ---\n")
print(head(df[, c("gene","baseMean","log2FoldChange","padj","sig")], 30), digits = 3)

# ---- naive(합치지 않은 90샘플) 결과와 비교 : pseudo-replication 영향 정량화
dn <- dds_naive[keep, ]; dn <- DESeq(dn, quiet = TRUE)
rn <- results(dn, contrast = c("genotype","HGPS","WT"), alpha = 0.05)
cat(sprintf("\n[비교] collapse 안 했을 때(90샘플) padj<0.05 : %d\n", sum(rn$padj < 0.05, na.rm = TRUE)))
cat(sprintf("[비교] collapse 했을 때(30세포주) padj<0.05  : %d\n", sum(df$padj < 0.05, na.rm = TRUE)))
cat("   -> 차이가 클수록 '반복을 잘못 세면 얼마나 부풀려지는지'를 보여준다\n")

# ---- 마커 점검
markers <- list(
  `Lamin/HGPS` = c("LMNA","LMNB1","LMNB2","ZMPSTE24"),
  `Senescence` = c("CDKN2A","CDKN1A","TP53","RB1","LMNB1","MKI67"),
  `SASP`       = c("IL6","IL1A","IL1B","CXCL8","CCL2","MMP1","MMP3","SERPINE1","IGFBP3"),
  `Interferon` = c("ISG15","IFI6","IFIT1","IFIT3","MX1","OAS1","STAT1","IRF7"),
  `ECM`        = c("COL1A1","COL1A2","COL3A1","FN1","ACTA2","TGFB1","CCN2"),
  `Cell cycle` = c("CCNB1","CDK1","PCNA","MCM2","TOP2A","E2F1")
)
cat("\n================ 마커 점검 ================\n")
for (k in names(markers)) {
  g <- intersect(markers[[k]], df$gene); if (!length(g)) next
  s <- df[g, c("gene","log2FoldChange","padj")]; s <- s[order(-abs(s$log2FoldChange)), ]
  cat(sprintf("\n[%s]\n", k)); print(s, digits = 3, row.names = FALSE)
}

write.table(df, file.path(RES, "DE_HGPS_vs_WT_patient.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- volcano
df$lab <- ifelse(df$sig != "NS" & rank(df$padj, ties.method="first") <= 30, df$gene, NA)
png(file.path(FIG, "volcano_patient.png"), width = 7, height = 6, units = "in", res = 300, type = "cairo")
print(ggplot(df[!is.na(df$padj), ], aes(log2FoldChange, -log10(padj), color = sig)) +
  geom_point(alpha = .5, size = 1.1) +
  scale_color_manual(values = c(UP="#C44E52", DOWN="#4C72B0", NS="grey80")) +
  geom_vline(xintercept = c(-LFC, LFC), linetype = 2, color = "grey40") +
  geom_hline(yintercept = -log10(FDR), linetype = 2, color = "grey40") +
  geom_text_repel(aes(label = lab), size = 2.6, max.overlaps = 40, show.legend = FALSE) +
  labs(title = "GSE306264: HGPS vs control fibroblasts (20 vs 10 donor lines)",
       subtitle = sprintf("FDR<0.05 & |log2FC|>=1  |  UP %d / DOWN %d",
                          sum(df$sig=="UP"), sum(df$sig=="DOWN")),
       x = "log2 fold change (HGPS / WT)", y = "-log10 adjusted p") +
  theme_bw(base_size = 11))
dev.off()

# =============================================================================
# PART 2. isogenic  TetON GFP-progerin  Dox6d vs noDox  (3 vs 3, 같은 세포)
# =============================================================================
cat("\n\n############ PART 2. isogenic (TetON GFP-progerin) ############\n")
ic <- col[col$arm == "isogenic", ]
ic$genotype <- factor(ic$genotype, levels = c("noDox", "Dox6d"))
im <- cts[, rownames(ic)]
keep2 <- rowSums(im >= 10) >= 3
ddi <- DESeqDataSetFromMatrix(im[keep2, ], ic, design = ~ genotype)
ddi <- DESeq(ddi, quiet = TRUE)
ri  <- lfcShrink(ddi, coef = "genotype_Dox6d_vs_noDox", type = "apeglm")
di  <- as.data.frame(ri); di$gene <- rownames(di)
di  <- di[order(di$padj, -abs(di$log2FoldChange)), ]
cat(sprintf("발현 필터 통과 %d 유전자\n", sum(keep2)))
cat(sprintf("padj<0.05            : %d\n", sum(di$padj < 0.05, na.rm = TRUE)))
cat(sprintf("padj<0.05 & |LFC|>=1 : %d\n",
            sum(di$padj < 0.05 & abs(di$log2FoldChange) >= 1, na.rm = TRUE)))
cat("\n--- 상위 20 ---\n")
print(head(di[, c("gene","baseMean","log2FoldChange","padj")], 20), digits = 3)
cat("\n[sanity] LMNA (progerin transgene 유도 확인):\n")
print(di[intersect(c("LMNA","GFP","CDKN1A","CDKN2A"), di$gene),
         c("gene","log2FoldChange","padj")], digits = 3, row.names = FALSE)

write.table(di, file.path(RES, "DE_isogenic_Dox6d_vs_noDox.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- 환자 결과와 isogenic 결과의 일치도
common <- intersect(df$gene[df$sig != "NS"], di$gene[!is.na(di$padj) & di$padj < 0.05])
if (length(common) > 5) {
  agree <- sign(df[common, "log2FoldChange"]) == sign(di[common, "log2FoldChange"])
  cat(sprintf("\n=== 환자 DEG 중 isogenic 에서도 유의한 유전자 : %d개, 방향 일치 %d개 (%.1f%%) ===\n",
              length(common), sum(agree), 100*mean(agree)))
  cr <- cor(df[common, "log2FoldChange"], di[common, "log2FoldChange"], method = "spearman")
  cat(sprintf("    log2FC Spearman 상관 : %.3f\n", cr))
}

saveRDS(list(dds = dds, vsd = vsd, res = df, iso = di, coldata = col),
        file.path(RES, "GSE306264_objects.rds"))
cat(sprintf("\n저장 완료 -> %s\n", RES))
