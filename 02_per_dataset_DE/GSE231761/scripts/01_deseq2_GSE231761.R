#!/usr/bin/env Rscript
# =============================================================================
# GSE231761 : HGPS 환자 유래 혈관평활근세포(VSMC) vs 정상 VSMC
#   4 세포주 (대조 0901C/BJ1C, HGPS 0031C/AG1B) x Passage 7/14 x 2 반복 = 16
#
# 설계 판단
#   - 같은 (세포주, passage) 의 반복 2개는 기술적 반복 -> collapse (16 -> 8)
#   - design = ~ passage + condition
#     passage 를 블록으로 잡아 계대 효과를 제거하고 질환 효과만 본다.
#   - 생물학적 단위(세포주)는 군당 2개뿐이므로 검정력은 낮다. 이는 메타분석에서
#     lfcSE 를 통해 자동으로 낮은 가중치로 반영된다.
# =============================================================================
suppressPackageStartupMessages({
  library(DESeq2); library(ggplot2); library(ggrepel); library(matrixStats)
  library(org.Hs.eg.db); library(AnnotationDbi)
})

# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
.this <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
DEPOSIT <- Sys.getenv("PROGERIA_DEPOSIT", normalizePath(file.path(dirname(.this), "..", "..", "..")))
BASE <- file.path(DEPOSIT, "02_per_dataset_DE", "GSE231761")  # [deposit] was ~/car_project/progeria_meta/GSE231761
RAW  <- file.path(BASE, "raw_data"); RES <- file.path(BASE, "results")
FIG  <- file.path(RES, "figures"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
set.seed(1)

# ------------------------------------------------------------ merge files ---
files <- sort(list.files(RAW, pattern = "^GSM.*\\.txt\\.gz$", full.names = TRUE))
cat(sprintf("featureCounts 파일 %d개\n", length(files)))
stopifnot(length(files) == 16)

read_fc <- function(p) {
  d <- read.delim(gzfile(p), skip = 1, check.names = FALSE)
  setNames(d[[ncol(d)]], d$Geneid)
}
lst <- lapply(files, read_fc)
gid <- names(lst[[1]])
stopifnot(all(vapply(lst, function(x) identical(names(x), gid), logical(1))))
cts <- do.call(cbind, lst)
rownames(cts) <- gid

# 파일명에서 샘플 정보 파싱: GSM7300149_0901C_P7_1.txt.gz
bn  <- sub("\\.txt\\.gz$", "", basename(files))
prt <- do.call(rbind, strsplit(bn, "_"))
col <- data.frame(
  sample    = paste(prt[,2], prt[,3], prt[,4], sep = "_"),
  gsm       = prt[,1],
  line      = prt[,2],
  passage   = prt[,3],
  tech_rep  = prt[,4],
  condition = ifelse(prt[,2] %in% c("0901C", "BJ1C"), "Control", "HGPS"),
  stringsAsFactors = FALSE
)
colnames(cts) <- col$sample
rownames(col) <- col$sample
cat("\n=== 샘플 구성 ===\n"); print(table(col$condition, col$line))
print(table(col$condition, col$passage))
write.table(col, file.path(BASE, "meta", "coldata.tsv"), sep = "\t", quote = FALSE)

# --------------------------------------------------- Ensembl -> gene symbol --
ens <- sub("\\..*$", "", rownames(cts))
sym <- mapIds(org.Hs.eg.db, keys = ens, column = "SYMBOL", keytype = "ENSEMBL",
              multiVals = "first")
cat(sprintf("\nEnsembl %d개 중 symbol 매핑 %d개\n", length(ens), sum(!is.na(sym))))

# --------------------------------------------------------------- collapse ---
storage.mode(cts) <- "integer"
dds0 <- DESeqDataSetFromMatrix(cts, col, design = ~ passage + condition)
grp  <- factor(paste(col$line, col$passage, sep = "_"))
dds  <- collapseReplicates(dds0, groupby = grp)
cat(sprintf("\ncollapse: %d -> %d (세포주 x passage)\n", ncol(dds0), ncol(dds)))
colData(dds)$condition <- factor(colData(dds)$condition, levels = c("Control", "HGPS"))
colData(dds)$passage   <- factor(colData(dds)$passage,   levels = c("P7", "P14"))
print(as.data.frame(colData(dds)[, c("line", "passage", "condition")]))

keep <- rowSums(counts(dds) >= 10) >= 4
cat(sprintf("발현 필터: %d / %d 유전자\n", sum(keep), nrow(dds)))
dds <- dds[keep, ]
design(dds) <- ~ passage + condition
dds <- DESeq(dds)

vsd <- vst(dds, blind = TRUE); vm <- assay(vsd)
pca <- prcomp(t(vm[order(rowVars(vm), decreasing = TRUE)[1:2000], ]))
pv  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
cat(sprintf("\n=== PCA PC1-3 : %.1f%% / %.1f%% / %.1f%% ===\n", pv[1], pv[2], pv[3]))

png(file.path(FIG, "QC_PCA.png"), width = 6.5, height = 5, units = "in", res = 300, type = "cairo")
print(ggplot(data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2],
                        condition = colData(dds)$condition,
                        passage = colData(dds)$passage, lab = colnames(dds)),
             aes(PC1, PC2, color = condition, shape = passage, label = lab)) +
  geom_point(size = 4) + geom_text_repel(size = 3, show.legend = FALSE) +
  scale_color_manual(values = c(Control = "#4C72B0", HGPS = "#C44E52")) +
  labs(x = sprintf("PC1 (%.1f%%)", pv[1]), y = sprintf("PC2 (%.1f%%)", pv[2]),
       title = "GSE231761: HGPS vs control VSMC (2 lines each, P7/P14)") +
  theme_bw(base_size = 11))
dev.off()

# -------------------------------------------------------------- DE : 질환 ---
res  <- results(dds, contrast = c("condition", "HGPS", "Control"), alpha = 0.05)
resS <- lfcShrink(dds, coef = "condition_HGPS_vs_Control", type = "apeglm", res = res)
df <- as.data.frame(resS)
df$ensembl <- rownames(df)
df$gene <- sym[sub("\\..*$", "", df$ensembl)]
df <- df[order(df$padj, -abs(df$log2FoldChange)), ]
df$sig <- ifelse(!is.na(df$padj) & df$padj < .05 & df$log2FoldChange >=  1, "UP",
          ifelse(!is.na(df$padj) & df$padj < .05 & df$log2FoldChange <= -1, "DOWN", "NS"))

cat("\n================ DE : HGPS vs Control VSMC (passage 보정) ================\n")
cat(sprintf("검정 유전자          : %d\n", sum(!is.na(df$padj))))
cat(sprintf("padj<0.05            : %d\n", sum(df$padj < .05, na.rm = TRUE)))
cat(sprintf("padj<0.05 & |LFC|>=1 : %d  (UP %d / DOWN %d)\n",
            sum(df$sig != "NS"), sum(df$sig == "UP"), sum(df$sig == "DOWN")))
cat("\n--- 상위 30 ---\n")
print(head(df[!is.na(df$gene), c("gene","baseMean","log2FoldChange","padj","sig")], 30), digits = 3)

# -------------------------------------------------------- DE : passage 효과 --
rp <- results(dds, contrast = c("passage", "P14", "P7"), alpha = 0.05)
cat(sprintf("\n[참고] passage P14 vs P7 효과 : padj<0.05 = %d 유전자\n",
            sum(rp$padj < .05, na.rm = TRUE)))
dp <- as.data.frame(rp); dp$ensembl <- rownames(dp)
dp$gene <- sym[sub("\\..*$", "", dp$ensembl)]
write.table(dp[order(dp$padj), ], file.path(RES, "DE_passage_P14_vs_P7.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# --------------------------------------------------------------- 마커 점검 --
markers <- list(
  `Lamin/HGPS` = c("LMNA","LMNB1","LMNB2","ZMPSTE24"),
  `Senescence` = c("CDKN2A","CDKN1A","TP53","RB1","LMNB1"),
  `SASP`       = c("IL6","IL1A","IL1B","CXCL8","CCL2","MMP1","MMP3","SERPINE1","IGFBP3"),
  `VSMC 표지`  = c("MYH11","ACTA2","CNN1","TAGLN","MYOCD","LMOD1"),
  `Interferon` = c("ISG15","IFI6","IFIT1","MX1","OAS1","STAT1")
)
cat("\n================ 마커 점검 ================\n")
for (k in names(markers)) {
  s <- df[!is.na(df$gene) & df$gene %in% markers[[k]], c("gene","log2FoldChange","padj")]
  if (!nrow(s)) next
  s <- s[order(-abs(s$log2FoldChange)), ]
  cat(sprintf("\n[%s]\n", k)); print(s, digits = 3, row.names = FALSE)
}

write.table(df, file.path(RES, "DE_HGPS_vs_Control_VSMC.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

df$lab <- ifelse(df$sig != "NS" & rank(df$padj, ties.method="first") <= 30, df$gene, NA)
png(file.path(FIG, "volcano.png"), width = 7, height = 6, units = "in", res = 300, type = "cairo")
print(ggplot(df[!is.na(df$padj), ], aes(log2FoldChange, -log10(padj), color = sig)) +
  geom_point(alpha = .5, size = 1.1) +
  scale_color_manual(values = c(UP="#C44E52", DOWN="#4C72B0", NS="grey80")) +
  geom_vline(xintercept = c(-1,1), linetype = 2, color = "grey40") +
  geom_hline(yintercept = -log10(.05), linetype = 2, color = "grey40") +
  geom_text_repel(aes(label = lab), size = 2.6, max.overlaps = 40, show.legend = FALSE) +
  labs(title = "GSE231761: HGPS vs control VSMC (passage-adjusted)",
       subtitle = sprintf("UP %d / DOWN %d", sum(df$sig=="UP"), sum(df$sig=="DOWN")),
       x = "log2 fold change (HGPS / Control)", y = "-log10 adjusted p") +
  theme_bw(base_size = 11))
dev.off()

saveRDS(list(dds = dds, vsd = vsd, res = df, coldata = col), file.path(RES, "GSE231761_objects.rds"))
cat(sprintf("\n저장 완료 -> %s\n", RES))
