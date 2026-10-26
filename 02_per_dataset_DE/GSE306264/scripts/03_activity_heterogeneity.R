#!/usr/bin/env Rscript
# =============================================================================
# GSE306264 : 샘플 수준 분석
#   (A) 성별 확인 (XIST / Y 유전자) 및 성별 보정 재검정
#   (B) 유전자세트 서명 점수 (샘플별 z-score) — 노화/SASP/IFN/EMT/평활근 등
#   (C) 경로 활성 (PROGENy) / 전사인자 활성 (DoRothEA + decoupleR)
#   (D) 환자 이질성 — 20개 HGPS 라인의 군집화, 하위군 탐색, 이상 라인 식별
#       (원 논문의 핵심 주장이 "환자 간 발현 이질성"이므로 반드시 재현해야 한다)
# =============================================================================
suppressPackageStartupMessages({
  library(DESeq2); library(ggplot2); library(ggrepel); library(pheatmap)
  library(RColorBrewer); library(matrixStats); library(msigdbr)
})

# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
.this <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
DEPOSIT <- Sys.getenv("PROGERIA_DEPOSIT", normalizePath(file.path(dirname(.this), "..", "..", "..")))
BASE <- file.path(DEPOSIT, "02_per_dataset_DE", "GSE306264")  # [deposit] was ~/car_project/progeria_meta/GSE306264
RES  <- file.path(BASE, "results"); FIG <- file.path(RES, "figures")
set.seed(1)

o   <- readRDS(file.path(RES, "GSE306264_objects.rds"))
dds <- o$dds; vm <- assay(o$vsd); geno <- colData(dds)$genotype
ann <- data.frame(genotype = geno, row.names = colnames(vm))
cat(sprintf("세포주 %d개 (HGPS %d / WT %d), 유전자 %d\n",
            ncol(vm), sum(geno == "HGPS"), sum(geno == "WT"), nrow(vm)))

# =============================================================================
# (A) 성별
# =============================================================================
cat("\n\n############ (A) 성별 확인 ############\n")
sexg <- c("XIST", "RPS4Y1", "DDX3Y", "UTY", "KDM5D", "USP9Y", "EIF1AY")
have <- intersect(sexg, rownames(vm))
cat("사용 가능한 성별 마커:", paste(have, collapse = ", "), "\n\n")
sx <- as.data.frame(t(vm[have, , drop = FALSE]))
ycols <- setdiff(have, "XIST")
sx$Y_mean <- if (length(ycols)) rowMeans(sx[, ycols, drop = FALSE]) else NA
sx$sex <- if ("XIST" %in% have && length(ycols))
            ifelse(sx$Y_mean > sx$XIST, "male", "female") else NA
sx$genotype <- geno
print(round(sx[, c(have, "Y_mean")], 2))
cat("\n=== 추정 성별 x 유전형 ===\n"); print(table(sx$sex, sx$genotype))
ft <- fisher.test(table(sx$sex, sx$genotype))
cat(sprintf("성별-유전형 연관 Fisher p = %.3f  %s\n", ft$p.value,
            ifelse(ft$p.value < .05, "<- 교란 가능! 보정 필요", "<- 균형, 교란 없음")))
write.table(sx, file.path(RES, "sex_assignment.tsv"), sep = "\t", quote = FALSE)

# 성별을 covariate 로 넣은 재검정
if (!any(is.na(sx$sex)) && length(unique(sx$sex)) == 2) {
  d2 <- dds; colData(d2)$sex <- factor(sx$sex)
  design(d2) <- ~ sex + genotype
  d2 <- DESeq(d2, quiet = TRUE)
  r2 <- results(d2, contrast = c("genotype", "HGPS", "WT"), alpha = .05)
  r1 <- results(dds, contrast = c("genotype", "HGPS", "WT"), alpha = .05)
  cat(sprintf("\n성별 미보정 padj<0.05 : %d\n", sum(r1$padj < .05, na.rm = TRUE)))
  cat(sprintf("성별 보정   padj<0.05 : %d\n", sum(r2$padj < .05, na.rm = TRUE)))
  cat(sprintf("log2FC 상관 (보정 전후) : %.4f\n",
              cor(r1$log2FoldChange, r2$log2FoldChange, use = "complete.obs")))
  dsex <- as.data.frame(r2); dsex$gene <- rownames(dsex)
  write.table(dsex[order(dsex$padj), ], file.path(RES, "DE_HGPS_vs_WT_sexAdjusted.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# =============================================================================
# (B) 서명 점수 (샘플별)
# =============================================================================
cat("\n\n############ (B) 유전자세트 서명 점수 ############\n")
H <- msigdbr(species = "Homo sapiens", collection = "H")
hs <- split(H$gene_symbol, H$gs_name)
custom <- list(
  SASP_core        = c("IL6","IL1A","IL1B","CXCL8","CCL2","MMP1","MMP3","MMP10",
                       "IGFBP3","IGFBP7","SERPINE1","TNF","CXCL1","CXCL2"),
  Senescence_CDKN  = c("CDKN1A","CDKN2A","CDKN2B","TP53","GLB1","TP53BP1"),
  Lamin_NE         = c("LMNA","LMNB1","LMNB2","LBR","EMD","TMPO","SUN1","SUN2","SYNE1","SYNE2"),
  SMC_program      = c("MYH11","ACTA2","CNN1","TAGLN","MYOCD","LMOD1","ACTG2","KCNMB1"),
  Ribosome         = grep("^RP[LS][0-9]", rownames(vm), value = TRUE)
)
sets <- c(hs[c("HALLMARK_MTORC1_SIGNALING","HALLMARK_GLYCOLYSIS",
               "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_MYOGENESIS",
               "HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_TNFA_SIGNALING_VIA_NFKB",
               "HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_P53_PATHWAY",
               "HALLMARK_ANGIOGENESIS","HALLMARK_HYPOXIA")], custom)

z <- t(scale(t(vm)))
score <- t(vapply(sets, function(g) {
  g <- intersect(g, rownames(z)); if (length(g) < 5) return(rep(NA_real_, ncol(z)))
  colMeans(z[g, , drop = FALSE])
}, numeric(ncol(z))))
colnames(score) <- colnames(vm)

st <- data.frame(set = rownames(score),
                 WT   = rowMeans(score[, geno == "WT",   drop = FALSE]),
                 HGPS = rowMeans(score[, geno == "HGPS", drop = FALSE]))
st$diff <- st$HGPS - st$WT
st$p <- apply(score, 1, function(x) if (all(is.na(x))) NA else
              wilcox.test(x[geno == "HGPS"], x[geno == "WT"])$p.value)
st$padj <- p.adjust(st$p, "BH")
st <- st[order(st$padj), ]
print(st, digits = 3, row.names = FALSE)
write.table(st, file.path(RES, "signature_scores_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(score, file.path(RES, "signature_scores_per_line.tsv"), sep = "\t", quote = FALSE)

png(file.path(FIG, "signature_scores_heatmap.png"), width = 8, height = 5.5, units = "in", res = 300, type = "cairo")
pheatmap(score[!is.na(rowSums(score)), ], annotation_col = ann,
         col = colorRampPalette(rev(brewer.pal(11, "RdBu")))(255),
         main = "Per-line signature scores (mean z of VST)", fontsize_row = 8)
dev.off()

# =============================================================================
# (C) 경로 활성 (PROGENy) / 전사인자 활성 (DoRothEA)
# =============================================================================
cat("\n\n############ (C) 경로 / 전사인자 활성 ############\n")
pg <- tryCatch({
  library(progeny)
  progeny(vm, scale = TRUE, organism = "Human", top = 500, perm = 1)
}, error = function(e) { cat("PROGENy 실패:", conditionMessage(e), "\n"); NULL })

if (!is.null(pg)) {
  pt <- data.frame(pathway = colnames(pg),
                   WT = colMeans(pg[geno == "WT", , drop = FALSE]),
                   HGPS = colMeans(pg[geno == "HGPS", , drop = FALSE]))
  pt$diff <- pt$HGPS - pt$WT
  pt$p <- apply(pg, 2, function(x) wilcox.test(x[geno=="HGPS"], x[geno=="WT"])$p.value)
  pt$padj <- p.adjust(pt$p, "BH"); pt <- pt[order(pt$padj), ]
  cat("\n=== PROGENy 경로 활성 (HGPS - WT) ===\n"); print(pt, digits = 3, row.names = FALSE)
  write.table(pt, file.path(RES, "progeny_pathway_activity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  png(file.path(FIG, "progeny_activity.png"), width = 7, height = 5, units = "in", res = 300, type = "cairo")
  pheatmap(t(pg), annotation_col = ann, main = "PROGENy pathway activity",
           col = colorRampPalette(rev(brewer.pal(11, "RdBu")))(255))
  dev.off()
}

tf <- tryCatch({
  library(decoupleR); library(dorothea)
  net <- dorothea::dorothea_hs
  net <- net[net$confidence %in% c("A","B","C"), ]
  names(net)[names(net) == "tf"] <- "source"
  names(net)[names(net) == "target"] <- "target"
  names(net)[names(net) == "mor"] <- "mor"
  decoupleR::run_ulm(mat = vm, net = net, .source = "source", .target = "target",
                     .mor = "mor", minsize = 10)
}, error = function(e) { cat("DoRothEA 실패:", conditionMessage(e), "\n"); NULL })

if (!is.null(tf)) {
  act <- reshape2::acast(as.data.frame(tf), source ~ condition, value.var = "score")
  act <- act[complete.cases(act), , drop = FALSE]
  g2 <- geno[match(colnames(act), colnames(vm))]
  tt <- data.frame(TF = rownames(act),
                   diff = rowMeans(act[, g2=="HGPS", drop=FALSE]) - rowMeans(act[, g2=="WT", drop=FALSE]),
                   p = apply(act, 1, function(x) wilcox.test(x[g2=="HGPS"], x[g2=="WT"])$p.value))
  tt$padj <- p.adjust(tt$p, "BH"); tt <- tt[order(tt$padj), ]
  cat(sprintf("\n=== 전사인자 활성: 검정 %d개, padj<0.05 %d개 ===\n", nrow(tt), sum(tt$padj < .05)))
  print(head(tt[order(-abs(tt$diff)), ][tt$padj[order(-abs(tt$diff))] < .1, ], 25),
        digits = 3, row.names = FALSE)
  write.table(tt, file.path(RES, "dorothea_TF_activity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  top <- head(tt$TF[order(tt$padj)], 30)
  png(file.path(FIG, "TF_activity_heatmap.png"), width = 8, height = 7, units = "in", res = 300, type = "cairo")
  pheatmap(act[top, ], annotation_col = ann, main = "DoRothEA TF activity (top 30)",
           col = colorRampPalette(rev(brewer.pal(11, "RdBu")))(255), fontsize_row = 8)
  dev.off()
}

# =============================================================================
# (D) 환자 이질성
# =============================================================================
cat("\n\n############ (D) 환자 간 이질성 ############\n")
hv <- vm[order(rowVars(vm), decreasing = TRUE)[1:3000], ]
hg <- hv[, geno == "HGPS"]; wt <- hv[, geno == "WT"]

cat(sprintf("HGPS 라인 간 평균 상관 : %.4f\n",
            mean(cor(hg, method = "spearman")[upper.tri(diag(ncol(hg)))])))
cat(sprintf("WT   라인 간 평균 상관 : %.4f\n",
            mean(cor(wt, method = "spearman")[upper.tri(diag(ncol(wt)))])))
cat("  -> HGPS 쪽이 더 낮으면 환자 간 이질성이 정상 대조보다 크다는 뜻\n")

# 각 HGPS 라인이 'HGPS 서명'을 얼마나 강하게 보이는가
de <- read.delim(file.path(RES, "DE_HGPS_vs_WT_patient.tsv"))
sigu <- de$gene[de$sig == "UP"]; sigd <- de$gene[de$sig == "DOWN"]
sigu <- intersect(sigu, rownames(z)); sigd <- intersect(sigd, rownames(z))
hgps_score <- colMeans(z[sigu, , drop = FALSE]) - colMeans(z[sigd, , drop = FALSE])
hs_df <- data.frame(line = names(hgps_score), score = hgps_score, genotype = geno)
hs_df <- hs_df[order(-hs_df$score), ]
cat("\n=== 라인별 HGPS 서명 점수 (UP평균 - DOWN평균) ===\n")
print(hs_df, digits = 3, row.names = FALSE)
write.table(hs_df, file.path(RES, "per_line_HGPS_score.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

ov <- sum(hs_df$score[hs_df$genotype == "HGPS"] < max(hs_df$score[hs_df$genotype == "WT"]))
cat(sprintf("\n정상 대조 최고점보다 낮은 HGPS 라인 : %d / %d\n", ov, sum(geno == "HGPS")))

# HGPS 라인 군집화 -> 하위군
hc <- hclust(as.dist(1 - cor(hg, method = "spearman")), method = "ward.D2")
k2 <- cutree(hc, k = 2); k3 <- cutree(hc, k = 3)
cat("\n=== HGPS 라인 군집 (k=2) ===\n"); print(split(names(k2), k2))
cat("\n=== HGPS 라인 군집 (k=3) ===\n"); print(split(names(k3), k3))
sub_ann <- data.frame(subgroup_k2 = factor(k2), subgroup_k3 = factor(k3),
                      HGPS_score = hgps_score[names(k2)], row.names = names(k2))
write.table(sub_ann, file.path(RES, "HGPS_subgroups.tsv"), sep = "\t", quote = FALSE)

png(file.path(FIG, "HGPS_line_clustering.png"), width = 7, height = 6, units = "in", res = 300, type = "cairo")
pheatmap(cor(hg, method = "spearman"), annotation_col = sub_ann["subgroup_k2"],
         annotation_row = sub_ann["HGPS_score"],
         main = "HGPS patient line similarity (top 3000 variable genes)",
         col = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(255), fontsize = 8)
dev.off()

png(file.path(FIG, "per_line_HGPS_score.png"), width = 7, height = 5, units = "in", res = 300, type = "cairo")
hs_df$line <- factor(hs_df$line, levels = hs_df$line)
print(ggplot(hs_df, aes(line, score, fill = genotype)) + geom_col() +
  scale_fill_manual(values = c(WT = "#4C72B0", HGPS = "#C44E52")) +
  labs(title = "HGPS transcriptional signature score per cell line",
       subtitle = "mean z(UP genes) - mean z(DOWN genes)", x = NULL, y = "signature score") +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 60, hjust = 1)))
dev.off()

# 하위군 간 차이 유전자
if (length(unique(k2)) == 2) {
  g1 <- names(k2)[k2 == 1]; g2n <- names(k2)[k2 == 2]
  d <- rowMeans(vm[, g1, drop=FALSE]) - rowMeans(vm[, g2n, drop=FALSE])
  p <- apply(vm, 1, function(x) tryCatch(wilcox.test(x[g1], x[g2n])$p.value, error=function(e) NA))
  sg <- data.frame(gene = rownames(vm), diff = d, p = p, padj = p.adjust(p, "BH"))
  sg <- sg[order(sg$padj, -abs(sg$diff)), ]
  cat(sprintf("\n=== HGPS 하위군1(n=%d) vs 하위군2(n=%d) 차이 유전자 padj<0.1 : %d ===\n",
              length(g1), length(g2n), sum(sg$padj < .1, na.rm = TRUE)))
  print(head(sg, 20), digits = 3, row.names = FALSE)
  write.table(sg, file.path(RES, "HGPS_subgroup_k2_DE.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
cat(sprintf("\n저장 완료 -> %s\n", RES))
