#!/usr/bin/env Rscript
# =============================================================================
# 4개 데이터셋(사람 3 + 마우스 1) 전부에서 공통인 결과만 추린다.
#   (A) 경로 : Hallmark — 마우스도 포함한 4-way 비교
#   (B) 유전자 : core16 / strict7
#   (C) LAD 대리지표 : 마우스에서도 재현되는가
# 3개(사람)에서만 확인된 것과 4개 전부에서 확인된 것을 엄격히 구분한다.
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

# =============================================================================
# (A) Hallmark 4-way
# =============================================================================
cat("############ (A) Hallmark 4-way 비교 ############\n")
Hh <- msigdbr(species="Homo sapiens", collection="H")
sh <- split(Hh$gene_symbol, Hh$gs_name); sh <- sh[lengths(sh) >= 15 & lengths(sh) <= 500]
Hm <- msigdbr(species="Mus musculus", collection="H")
sm <- split(Hm$gene_symbol, Hm$gs_name); sm <- sm[lengths(sm) >= 15 & lengths(sm) <= 500]

# 마우스 Hallmark 를 우리 기준(WT_24mo 제외본)으로 직접 계산
mo <- read.delim(file.path(DEPOSIT, "01_rnaseq_processing/GSE165409_mouse_aorta/results/08_deseq2/DE_LMNA_vs_WTyng.tsv"))  # [deposit]
mo <- mo[!is.na(mo$pvalue) & !is.na(mo$symbol) & mo$symbol != "", ]
mo$stat <- sign(mo$log2FoldChange) * -log10(pmax(mo$pvalue, 1e-300))
mo <- mo[order(-abs(mo$stat)), ]; mo <- mo[!duplicated(mo$symbol), ]
rm_ <- sort(setNames(mo$stat, mo$symbol), decreasing=TRUE)
gm <- fgsea(sm, rm_, minSize=15, maxSize=500, eps=0)
cat(sprintf("마우스 Hallmark: %d 경로 (유의 %d)\n", nrow(gm), sum(gm$padj < .05)))
om <- as.data.frame(gm); om$leadingEdge <- vapply(om$leadingEdge, paste, character(1), collapse=",")
write.table(om, file.path(RES, "GSEA_GSE165409_Hallmark_full.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

a <- read.delim(file.path(ROOT, "GSE306264/results/GSEA_patient_Hallmark.tsv"))
b <- read.delim(file.path(ROOT, "GSE231761/results/GSEA_condition_Hallmark.tsv"))
c3 <- read.delim(file.path(RES, "GSEA_GSE198983_Hallmark_full.tsv"))

M <- Reduce(function(x, y) merge(x, y, by="pathway"), list(
  data.frame(pathway=a$pathway,  NES_306264=a$NES,  padj_306264=a$padj),
  data.frame(pathway=b$pathway,  NES_231761=b$NES,  padj_231761=b$padj),
  data.frame(pathway=c3$pathway, NES_198983=c3$NES, padj_198983=c3$padj),
  data.frame(pathway=om$pathway, NES_165409=om$NES, padj_165409=om$padj)))
NC <- c("NES_306264","NES_231761","NES_198983","NES_165409")
PC <- c("padj_306264","padj_231761","padj_198983","padj_165409")
M$n_sig <- rowSums(M[, PC] < .05, na.rm=TRUE)
M$same4 <- apply(sign(M[, NC]), 1, function(x) length(unique(x)) == 1)
M$same3_human <- apply(sign(M[, NC[1:3]]), 1, function(x) length(unique(x)) == 1)
M$mean_NES <- rowMeans(M[, NC])
M <- M[order(-M$n_sig, -abs(M$mean_NES)), ]
write.table(M, file.path(RES, "TABLE3_hallmark_4way.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf("병합된 경로: %d\n", nrow(M)))

cat("\n=== [4-way] 4개 전부 방향 일치 + 2개 이상 유의 ===\n")
four <- M[M$same4 & M$n_sig >= 2, ]
if (nrow(four)) print(four[, c("pathway", NC, "n_sig")], digits=3, row.names=FALSE) else cat("  없음\n")

cat("\n=== [4-way] 4개 전부 방향 일치 (유의도 무관) 상위 12 ===\n")
f2 <- M[M$same4, ]; f2 <- f2[order(-abs(f2$mean_NES)), ]
print(head(f2[, c("pathway", NC, "n_sig")], 12), digits=3, row.names=FALSE)

cat("\n=== [참고] 사람 3개만 일치했던 4개 경로가 마우스에서는? ===\n")
key <- c("HALLMARK_MTORC1_SIGNALING","HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
         "HALLMARK_HEDGEHOG_SIGNALING","HALLMARK_UNFOLDED_PROTEIN_RESPONSE")
print(M[M$pathway %in% key, c("pathway", NC, "n_sig", "same4")], digits=3, row.names=FALSE)

# =============================================================================
# (B) 유전자 4-way
# =============================================================================
cat("\n\n############ (B) 유전자 4-way ############\n")
core <- read.delim(file.path(RES, "candidates_k3_final.tsv"))
strict <- read.delim(file.path(RES, "candidates_k3_strict.tsv"))
cat(sprintf("사람 3개 일치 + |LFC|>=1 + 마우스 방향 일치 : %d개\n", nrow(core)))
cat(sprintf("  그중 마우스에서도 FDR<0.05                 : %d개\n", nrow(strict)))
cat("\n[4-way 최고 신뢰 집합 — 4개 데이터셋 전부에서 유의하고 방향 일치]\n")
print(strict[, c("gene","meta_lfc","meta_padj","I2","lfc_306264","lfc_231761",
                 "lfc_198983","lfc_mouse","padj_mouse","chronicity")],
      digits=3, row.names=FALSE)

# =============================================================================
# (C) LAD 4-way
# =============================================================================
cat("\n\n############ (C) LAD 대리지표 4-way ############\n")
fc <- read.delim(file.path(DEPOSIT, "01_rnaseq_processing/GSE165409_mouse_aorta/results/05_counts/counts_raw.txt"),  # [deposit]
                 skip=1, check.names=FALSE)[, 1:6]
pos <- data.frame(gene_id = sub("\\..*$", "", fc$Geneid),
                  chr = vapply(strsplit(as.character(fc$Chr), ";"), function(x) x[1], character(1)),
                  start = vapply(strsplit(as.character(fc$Start), ";"), function(x) as.numeric(x[1]), numeric(1)),
                  end = vapply(strsplit(as.character(fc$End), ";"), function(x) max(as.numeric(x)), numeric(1)),
                  stringsAsFactors = FALSE)
mainc <- paste0("chr", c(1:19, "X"))
pos <- pos[pos$chr %in% mainc, ]; pos$mid <- (pos$start + pos$end)/2
mo2 <- read.delim(file.path(DEPOSIT, "01_rnaseq_processing/GSE165409_mouse_aorta/results/08_deseq2/DE_LMNA_vs_WTyng.tsv"))  # [deposit]
mo2$gene_id <- sub("\\..*$", "", mo2$gene_id)
mo2$sig <- ifelse(!is.na(mo2$padj) & mo2$padj < .05 & mo2$log2FoldChange >=  1, "UP",
           ifelse(!is.na(mo2$padj) & mo2$padj < .05 & mo2$log2FoldChange <= -1, "DOWN", "NS"))
d <- merge(mo2[, c("gene_id","symbol","log2FoldChange","padj","sig")], pos, by="gene_id")
d <- d[!is.na(d$padj), ]; d$chr <- factor(d$chr, levels=mainc)
cat(sprintf("마우스 좌표 확보 유전자 %d (UP %d / DOWN %d)\n",
            nrow(d), sum(d$sig=="UP"), sum(d$sig=="DOWN")))
dens <- vapply(seq_len(nrow(d)), function(i){
  s <- d[d$chr == d$chr[i], ]; sum(abs(s$mid - d$mid[i]) <= 2e6) - 1 }, numeric(1))
d$gene_density <- dens
d$dens_decile <- cut(rank(d$gene_density, ties.method="first"),
                     breaks=quantile(rank(d$gene_density, ties.method="first"), probs=seq(0,1,.1)),
                     include.lowest=TRUE, labels=1:10)
lo <- d$dens_decile %in% c("1","2"); hi <- d$dens_decile %in% c("9","10")
wt <- wilcox.test(d$log2FoldChange[lo], d$log2FoldChange[hi])
cat(sprintf("희소(1-2) 평균 log2FC = %.4f  vs  조밀(9-10) = %.4f,  p = %.3g\n",
            mean(d$log2FoldChange[lo]), mean(d$log2FoldChange[hi]), wt$p.value))
cat(sprintf("UP 비율: 희소 %.4f vs 조밀 %.4f\n", mean(d$sig[lo]=="UP"), mean(d$sig[hi]=="UP")))
write.table(d, file.path(RES, "GSE165409_positional_gene_table.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

mk <- function(x, lab) data.frame(dataset=lab, decile=as.integer(as.character(x$dens_decile)),
                                  up = x$sig=="UP", lfc = x$log2FoldChange)
ld <- rbind(
  mk(read.delim(file.path(ROOT,"GSE306264/results/positional_gene_table.tsv")), "GSE306264 (human fibroblast)"),
  mk(read.delim(file.path(ROOT,"GSE231761/results/positional_gene_table.tsv")), "GSE231761 (human VSMC)"),
  mk(read.delim(file.path(DEPOSIT, "01_rnaseq_processing/GSE198983_human_AoSMC/results/19_positional/positional_gene_table.tsv")), "GSE198983 (human isogenic VSMC)"),  # [deposit]
  mk(d, "GSE165409 (mouse aorta)"))
agg <- aggregate(cbind(up, lfc) ~ dataset + decile, ld, mean)
cols <- c("GSE306264 (human fibroblast)"="#C44E52", "GSE231761 (human VSMC)"="#4C72B0",
          "GSE198983 (human isogenic VSMC)"="#55A868", "GSE165409 (mouse aorta)"="#8172B3")
PNG(file.path(FIG, "LAD_proxy_4way.png"), 10.5, 4.6)
gridExtra::grid.arrange(
  ggplot(agg, aes(decile, 100*up, color=dataset)) + geom_line(linewidth=.9) + geom_point(size=2.2) +
    scale_x_continuous(breaks=1:10) + scale_color_manual(values=cols, name=NULL) +
    labs(title="Upregulated-gene fraction vs local gene density",
         subtitle="decile 1 = gene-poor (LAD-like)", x="gene density decile", y="% genes upregulated") +
    theme_bw(base_size=10) + theme(legend.position="bottom") + guides(color=guide_legend(nrow=2)),
  ggplot(agg, aes(decile, lfc, color=dataset)) + geom_line(linewidth=.9) + geom_point(size=2.2) +
    geom_hline(yintercept=0, linetype=2, color="grey50") +
    scale_x_continuous(breaks=1:10) + scale_color_manual(values=cols, name=NULL) +
    labs(title="Mean log2FC vs local gene density",
         subtitle="reproduced in three human datasets and in mouse aorta",
         x="gene density decile", y="mean log2FC") +
    theme_bw(base_size=10) + theme(legend.position="bottom") + guides(color=guide_legend(nrow=2)),
  ncol=2)
dev.off()
cat("  -> LAD_proxy_4way.png 저장\n")

# 4-way Hallmark 히트맵
hm <- M[order(-M$mean_NES), ]
l2 <- melt(hm[, c("pathway", NC)], id.vars="pathway", variable.name="dataset", value.name="NES")
sg <- melt(hm[, c("pathway", PC)], id.vars="pathway", value.name="padj")
l2$star <- ifelse(sg$padj < .05, "*", "")
l2$dataset <- factor(sub("^NES_", "GSE", l2$dataset),
                     levels=c("GSE306264","GSE231761","GSE198983","GSE165409"))
l2$pathway <- factor(sub("^HALLMARK_", "", l2$pathway),
                     levels=rev(sub("^HALLMARK_", "", hm$pathway)))
PNG(file.path(FIG, "hallmark_heatmap_4way.png"), 7.8, 10)
print(ggplot(l2, aes(dataset, pathway, fill=NES)) +
  geom_tile(color="white", linewidth=.3) + geom_text(aes(label=star), size=4, vjust=.78) +
  scale_fill_gradient2(low="#4C72B0", mid="white", high="#C44E52", midpoint=0) +
  labs(title="Hallmark GSEA across all four datasets", subtitle="* FDR < 0.05", x=NULL, y=NULL) +
  theme_minimal(base_size=8) +
  theme(axis.text.x=element_text(angle=30, hjust=1), panel.grid=element_blank()))
dev.off()
cat("  -> hallmark_heatmap_4way.png 저장\n")

if (nrow(four) > 0) {
  lg <- melt(four[, c("pathway", NC)], id.vars="pathway", variable.name="dataset", value.name="NES")
  lg$dataset <- factor(sub("^NES_", "GSE", lg$dataset),
                       levels=c("GSE306264","GSE231761","GSE198983","GSE165409"))
  ordp <- four$pathway[order(four$mean_NES)]
  lg$pathway <- factor(sub("^HALLMARK_", "", lg$pathway), levels=sub("^HALLMARK_", "", ordp))
  PNG(file.path(FIG, "common_hallmark_4way.png"), 9.5, max(3.4, .62*nrow(four)+1.8))
  print(ggplot(lg, aes(NES, pathway, fill=dataset)) +
    geom_col(position=position_dodge(.8), width=.74) + geom_vline(xintercept=0, color="grey35") +
    scale_fill_manual(values=c("#C44E52","#4C72B0","#55A868","#8172B3"), name=NULL) +
    labs(title="Hallmark pathways concordant across all four datasets",
         subtitle="same direction in 4/4 (three human + mouse aorta), FDR<0.05 in at least two",
         x="Normalized Enrichment Score", y=NULL) +
    theme_bw(base_size=10) + theme(legend.position="bottom"))
  dev.off()
  cat("  -> common_hallmark_4way.png 저장\n")
}
cat(sprintf("\n저장 -> %s\n", RES))
