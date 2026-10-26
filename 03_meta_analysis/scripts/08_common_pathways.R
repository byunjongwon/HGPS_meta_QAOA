#!/usr/bin/env Rscript
# =============================================================================
# 사람 3개 데이터셋 Hallmark 경로 통합 비교 (수정판)
#
# 07번에서 GSE198983 의 GSEA 파일이 '유의한 16개'만 담고 있어 비교 대상이
# 50개 중 16개로 줄어 있었다. 여기서는 GSE198983 의 blocked DE 결과로부터
# Hallmark 50개 전체를 직접 계산해 공정하게 비교한다.
# =============================================================================
suppressPackageStartupMessages({
  library(fgsea); library(msigdbr); library(ggplot2); library(reshape2)
})
# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
.this <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
DEPOSIT <- Sys.getenv("PROGERIA_DEPOSIT", normalizePath(file.path(dirname(.this), "..", "..")))
ROOT <- file.path(DEPOSIT, "02_per_dataset_DE")  # [deposit] GSE* dataset dirs (was ~/car_project/progeria_meta)
RES  <- file.path(DEPOSIT, "03_meta_analysis/results"); FIG <- file.path(RES, "figures")  # [deposit]
PNG <- function(f, w, h) png(f, width=w, height=h, units="in", res=300, type="cairo")
set.seed(1)

H <- msigdbr(species = "Homo sapiens", collection = "H")
sets <- split(H$gene_symbol, H$gs_name)
sets <- sets[lengths(sets) >= 15 & lengths(sets) <= 500]

# --- GSE198983 : blocked DE -> Hallmark 50개 전체
d3 <- read.delim(file.path(DEPOSIT, "01_rnaseq_processing/GSE198983_human_AoSMC/results/09_deseq2/DE_Progerin_vs_Control_blocked.tsv"))  # [deposit]
d3 <- d3[!is.na(d3$pvalue) & !is.na(d3$symbol) & d3$symbol != "", ]
d3$stat <- sign(d3$log2FoldChange) * -log10(pmax(d3$pvalue, 1e-300))
d3 <- d3[order(-abs(d3$stat)), ]; d3 <- d3[!duplicated(d3$symbol), ]
r3 <- sort(setNames(d3$stat, d3$symbol), decreasing = TRUE)
g3 <- fgsea(sets, r3, minSize = 15, maxSize = 500, eps = 0)
cat(sprintf("GSE198983 Hallmark 재계산: %d 경로 (유의 %d)\n", nrow(g3), sum(g3$padj < .05)))
o3 <- as.data.frame(g3); o3$leadingEdge <- vapply(o3$leadingEdge, paste, character(1), collapse=",")
write.table(o3, file.path(RES, "GSEA_GSE198983_Hallmark_full.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

a <- read.delim(file.path(ROOT, "GSE306264/results/GSEA_patient_Hallmark.tsv"))
b <- read.delim(file.path(ROOT, "GSE231761/results/GSEA_condition_Hallmark.tsv"))

m <- Reduce(function(x, y) merge(x, y, by = "pathway"), list(
  data.frame(pathway = a$pathway, NES_306264 = a$NES, padj_306264 = a$padj),
  data.frame(pathway = b$pathway, NES_231761 = b$NES, padj_231761 = b$padj),
  data.frame(pathway = o3$pathway, NES_198983 = o3$NES, padj_198983 = o3$padj)))
cat(sprintf("병합된 Hallmark 경로: %d\n", nrow(m)))

NC <- c("NES_306264","NES_231761","NES_198983")
PC <- c("padj_306264","padj_231761","padj_198983")
m$n_sig <- rowSums(m[, PC] < .05, na.rm = TRUE)
m$all_same_dir <- apply(sign(m[, NC]), 1, function(x) length(unique(x)) == 1)
m$mean_NES <- rowMeans(m[, NC])
m <- m[order(-m$n_sig, -abs(m$mean_NES)), ]
write.table(m, file.path(RES, "TABLE2_hallmark_across_datasets.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

cat("\n=== 3개 전부 유의 + 방향 일치 ===\n")
core <- m[m$n_sig == 3 & m$all_same_dir, ]
if (nrow(core)) print(core[, c("pathway", NC, "mean_NES")], digits=3, row.names=FALSE) else cat("  없음\n")

cat("\n=== 3개 방향 일치 + 2개 이상 유의 (핵심 공통 경로) ===\n")
com <- m[m$all_same_dir & m$n_sig >= 2, ]
print(com[, c("pathway", NC, "n_sig", "mean_NES")], digits=3, row.names=FALSE)

cat("\n=== 3개 방향 일치 (유의도 무관) 상위 12 ===\n")
allc <- m[m$all_same_dir, ]
print(head(allc[order(-abs(allc$mean_NES)), c("pathway", NC, "n_sig")], 12), digits=3, row.names=FALSE)

cat("\n=== 환자 2개는 같은데 isogenic 만 반대인 경로 (급성 vs 만성) ===\n")
opp <- m[sign(m$NES_306264) == sign(m$NES_231761) &
         sign(m$NES_306264) != sign(m$NES_198983) & m$n_sig >= 2, ]
print(head(opp[order(-abs(opp$NES_306264)), c("pathway", NC, "n_sig")], 12), digits=3, row.names=FALSE)

# ------------------------------------------------------------------ 그림
plot_set <- com
if (nrow(plot_set) > 0) {
  long <- melt(plot_set[, c("pathway", NC)], id.vars="pathway",
               variable.name="dataset", value.name="NES")
  long$dataset <- factor(sub("^NES_", "GSE", long$dataset),
                         levels=c("GSE306264","GSE231761","GSE198983"),
                         labels=c("GSE306264 (fibroblast)","GSE231761 (VSMC)","GSE198983 (isogenic VSMC)"))
  ord <- plot_set$pathway[order(plot_set$mean_NES)]
  long$pathway <- factor(sub("^HALLMARK_", "", long$pathway),
                         levels=sub("^HALLMARK_", "", ord))
  PNG(file.path(FIG, "common_hallmark_pathways.png"), 9, max(3.2, .55*nrow(plot_set)+1.6))
  print(ggplot(long, aes(NES, pathway, fill=dataset)) +
    geom_col(position=position_dodge(.78), width=.72) +
    geom_vline(xintercept=0, color="grey35") +
    scale_fill_manual(values=c("#C44E52","#DD8452","#55A868"), name=NULL) +
    labs(title="Hallmark pathways concordant across all three human datasets",
         subtitle="same direction in 3/3 datasets, FDR<0.05 in at least 2",
         x="Normalized Enrichment Score", y=NULL) +
    theme_bw(base_size=10) + theme(legend.position="bottom"))
  dev.off()
  cat("\n  -> common_hallmark_pathways.png 저장\n")
}

# 전체 Hallmark 히트맵 (부록용)
hm <- m[order(-m$mean_NES), ]
long2 <- melt(hm[, c("pathway", NC)], id.vars="pathway",
              variable.name="dataset", value.name="NES")
long2$dataset <- factor(sub("^NES_", "GSE", long2$dataset),
                        levels=c("GSE306264","GSE231761","GSE198983"))
long2$pathway <- factor(sub("^HALLMARK_", "", long2$pathway),
                        levels=rev(sub("^HALLMARK_", "", hm$pathway)))
sigdf <- melt(hm[, c("pathway", PC)], id.vars="pathway", value.name="padj")
long2$star <- ifelse(sigdf$padj < .05, "*", "")
PNG(file.path(FIG, "hallmark_heatmap_all.png"), 7.5, 10)
print(ggplot(long2, aes(dataset, pathway, fill=NES)) +
  geom_tile(color="white", linewidth=.3) +
  geom_text(aes(label=star), size=4, vjust=.78) +
  scale_fill_gradient2(low="#4C72B0", mid="white", high="#C44E52", midpoint=0) +
  labs(title="Hallmark GSEA across the three human datasets",
       subtitle="* FDR < 0.05", x=NULL, y=NULL) +
  theme_minimal(base_size=8) +
  theme(axis.text.x=element_text(angle=30, hjust=1), panel.grid=element_blank()))
dev.off()
cat("  -> hallmark_heatmap_all.png 저장\n")
cat(sprintf("\n저장 -> %s\n", RES))
