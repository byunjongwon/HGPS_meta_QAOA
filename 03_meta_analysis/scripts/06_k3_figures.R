#!/usr/bin/env Rscript
# =============================================================================
# k=3 확정 후 그림 갱신 + GSE198983 위치기반(LAD) 분석 추가
#
#  [1] GSE198983 위치 기반 / LAD 대리지표  -> progerin_rnaseq/results/19_positional/
#  [2] LAD proxy 3개 데이터셋 통합 비교 (갱신)
#  [3] 필터 깔때기 (k=3)
#  [4] forest plot (4개 데이터셋 열)
#  [5] UpSet (5개 집합)
#  [6] 데이터셋 3쌍 일치도 산점도
#  [7] 이질성 I2 분포
#  [8] 데이터셋 요약 비교 (5개)
# =============================================================================
suppressPackageStartupMessages({
  library(DESeq2); library(ggplot2); library(ggrepel); library(ComplexHeatmap)
  library(RColorBrewer); library(reshape2); library(gridExtra)
})
# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
.this <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
DEPOSIT <- Sys.getenv("PROGERIA_DEPOSIT", normalizePath(file.path(dirname(.this), "..", "..")))
ROOT <- file.path(DEPOSIT, "02_per_dataset_DE")  # [deposit] GSE* dataset dirs (was ~/car_project/progeria_meta)
PR   <- file.path(DEPOSIT, "01_rnaseq_processing/GSE198983_human_AoSMC")  # [deposit] was ~/progerin_rnaseq
RES  <- file.path(DEPOSIT, "03_meta_analysis/results"); MFIG <- file.path(RES, "figures")  # [deposit]
P19  <- file.path(PR, "results/19_positional"); dir.create(P19, showWarnings = FALSE, recursive = TRUE)
set.seed(1)
PNG <- function(f, w, h) png(f, width = w, height = h, units = "in", res = 300, type = "cairo")

M <- read.delim(file.path(RES, "meta_k3_final.tsv"))

# =============================================================================
# [1] GSE198983 위치 기반 / LAD
# =============================================================================
cat("############ [1] GSE198983 LAD 분석 ############\n")
fc <- read.delim(file.path(PR, "results/07_counts/counts_raw.txt"), skip = 1, check.names = FALSE)[, 1:6]
pos <- data.frame(
  # featureCounts 는 GENCODE 버전 접미사(ENSG...\.2)를 갖고, DE 표는 버전이 제거되어 있다
  gene_id = sub("\\..*$", "", fc$Geneid),
  chr   = vapply(strsplit(as.character(fc$Chr),   ";"), function(x) x[1], character(1)),
  start = vapply(strsplit(as.character(fc$Start), ";"), function(x) as.numeric(x[1]), numeric(1)),
  end   = vapply(strsplit(as.character(fc$End),   ";"), function(x) max(as.numeric(x)), numeric(1)),
  stringsAsFactors = FALSE)
main <- paste0("chr", c(1:22, "X"))
pos <- pos[pos$chr %in% main, ]; pos$mid <- (pos$start + pos$end)/2

de3 <- read.delim(file.path(PR, "results/09_deseq2/DE_Progerin_vs_Control_blocked.tsv"))
de3$sig <- ifelse(!is.na(de3$padj) & de3$padj < .05 & de3$log2FoldChange >=  1, "UP",
           ifelse(!is.na(de3$padj) & de3$padj < .05 & de3$log2FoldChange <= -1, "DOWN", "NS"))
d <- merge(de3[, c("gene_id","symbol","log2FoldChange","padj","sig")], pos, by = "gene_id")
d <- d[!is.na(d$padj), ]; d$chr <- factor(d$chr, levels = main)
cat(sprintf("좌표 있는 검정 유전자 %d (UP %d / DOWN %d)\n", nrow(d), sum(d$sig=="UP"), sum(d$sig=="DOWN")))

dens <- vapply(seq_len(nrow(d)), function(i) {
  s <- d[d$chr == d$chr[i], ]; sum(abs(s$mid - d$mid[i]) <= 2e6) - 1 }, numeric(1))
d$gene_density <- dens
d$dens_decile <- cut(rank(d$gene_density, ties.method="first"),
                     breaks = quantile(rank(d$gene_density, ties.method="first"), probs = seq(0,1,.1)),
                     include.lowest = TRUE, labels = 1:10)
lo <- d$dens_decile %in% c("1","2"); hi <- d$dens_decile %in% c("9","10")
wt <- wilcox.test(d$log2FoldChange[lo], d$log2FoldChange[hi])
cat(sprintf("희소(1-2) 평균 log2FC = %.4f  vs  조밀(9-10) = %.4f,  p = %.3g\n",
            mean(d$log2FoldChange[lo]), mean(d$log2FoldChange[hi]), wt$p.value))
cat(sprintf("UP 비율: 희소 %.4f vs 조밀 %.4f\n", mean(d$sig[lo]=="UP"), mean(d$sig[hi]=="UP")))
write.table(d, file.path(P19, "positional_gene_table.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

PNG(file.path(P19, "LAD_proxy_density_vs_lfc.png"), 7, 5)
print(ggplot(d, aes(dens_decile, log2FoldChange)) + geom_boxplot(outlier.size=.3, fill="#B0C4DE") +
  geom_hline(yintercept=0, linetype=2, color="red") +
  labs(title="Expression change vs local gene density (LAD proxy)",
       subtitle=sprintf("GSE198983 isogenic VSMC; sparse(1-2) %.3f vs dense(9-10) %.3f, p=%.2g",
                        mean(d$log2FoldChange[lo]), mean(d$log2FoldChange[hi]), wt$p.value),
       x="local gene density decile (+-2 Mb)", y="log2FC (progerin / control)") +
  coord_cartesian(ylim=quantile(d$log2FoldChange, c(.005,.995))) + theme_bw(base_size=11))
dev.off()

# =============================================================================
# [2] LAD proxy 3개 데이터셋 통합
# =============================================================================
cat("\n[2] LAD proxy 3개 통합\n")
mk <- function(x, lab, lfc = "log2FoldChange")
  data.frame(dataset = lab, decile = as.integer(as.character(x$dens_decile)),
             up = x$sig == "UP", lfc = x[[lfc]])
ld <- rbind(
  mk(read.delim(file.path(ROOT,"GSE306264/results/positional_gene_table.tsv")), "GSE306264 (fibroblast)"),
  mk(read.delim(file.path(ROOT,"GSE231761/results/positional_gene_table.tsv")), "GSE231761 (VSMC)"),
  mk(d, "GSE198983 (isogenic VSMC)"))
agg <- aggregate(cbind(up, lfc) ~ dataset + decile, ld, mean)
cols <- c(`GSE306264 (fibroblast)`="#C44E52", `GSE231761 (VSMC)`="#4C72B0",
          `GSE198983 (isogenic VSMC)`="#55A868")
PNG(file.path(MFIG, "LAD_proxy_across_datasets.png"), 10, 4.5)
grid.arrange(
  ggplot(agg, aes(decile, 100*up, color = dataset)) + geom_line(linewidth=.9) + geom_point(size=2.2) +
    scale_x_continuous(breaks=1:10) + scale_color_manual(values=cols, name=NULL) +
    labs(title="Upregulated-gene fraction vs local gene density",
         subtitle="decile 1 = gene-poor (LAD-like)", x="gene density decile", y="% genes upregulated") +
    theme_bw(base_size=10) + theme(legend.position="bottom"),
  ggplot(agg, aes(decile, lfc, color = dataset)) + geom_line(linewidth=.9) + geom_point(size=2.2) +
    geom_hline(yintercept=0, linetype=2, color="grey50") +
    scale_x_continuous(breaks=1:10) + scale_color_manual(values=cols, name=NULL) +
    labs(title="Mean log2FC vs local gene density",
         subtitle="reproduced independently in three human datasets",
         x="gene density decile", y="mean log2FC") +
    theme_bw(base_size=10) + theme(legend.position="bottom"), ncol=2)
dev.off()

# =============================================================================
# [3] 깔때기 (k=3)
# =============================================================================
cat("[3] 깔때기\n")
base <- M$meta_padj < .05 & M$all_same_dir
fn <- data.frame(
  step = c("common genes (3 human sets)","meta FDR<0.05","+ all 3 concordant","+ I2 < 75%",
           "+ |meta log2FC| >= 0.5","+ |meta log2FC| >= 1",
           "+ mouse same direction","+ mouse significant & same direction"),
  n = c(nrow(M), sum(M$meta_padj<.05,na.rm=TRUE), sum(base,na.rm=TRUE), sum(base & M$I2<75,na.rm=TRUE),
        sum(base & M$I2<75 & abs(M$meta_lfc)>=.5,na.rm=TRUE),
        sum(base & M$I2<75 & abs(M$meta_lfc)>=1,na.rm=TRUE),
        sum(base & M$I2<75 & abs(M$meta_lfc)>=1 & M$mouse_dir_only,na.rm=TRUE),
        sum(base & M$I2<75 & abs(M$meta_lfc)>=1 & M$mouse_same_dir,na.rm=TRUE)))
fn$step <- factor(fn$step, levels = rev(fn$step)); fn$lab <- format(fn$n, big.mark=",")
PNG(file.path(MFIG, "filter_funnel.png"), 9.5, 5.5)
print(ggplot(fn, aes(n, step)) +
  annotate("rect", xmin=20, xmax=54, ymin=-Inf, ymax=Inf, alpha=.16, fill="#C44E52") +
  annotate("segment", x=c(20,54), xend=c(20,54), y=-Inf, yend=Inf, linetype=2, color="#C44E52", linewidth=.4) +
  geom_col(fill="#4C72B0", width=.68) + geom_text(aes(label=lab), hjust=-0.18, size=3.4) +
  scale_x_log10(expand=expansion(mult=c(0,.20))) +
  scale_y_discrete(expand=expansion(add=c(0.6,1.15))) +
  annotate("label", x=33, y=8.85, label="IQM hardware range\n20 (Garnet) - 54 (Emerald) qubits",
           size=2.9, color="#C44E52", fill="white", lineheight=.95) +
  labs(title="Candidate gene funnel for the QUBO panel (k = 3 random-effects meta-analysis)",
       subtitle="GSE306264 + GSE231761 + GSE198983, validated against GSE165409 (mouse aorta)",
       x="number of genes (log scale)", y=NULL) + theme_bw(base_size=11))
dev.off()

# =============================================================================
# [4] forest plot
# =============================================================================
cat("[4] forest plot\n")
cand <- read.delim(file.path(RES, "candidates_k3_final.tsv"))
tp <- head(cand[order(-abs(cand$meta_lfc)), ], 16)
fr <- rbind(
  data.frame(gene=tp$gene, dataset="GSE306264 (fibroblast)",     lfc=tp$lfc_306264),
  data.frame(gene=tp$gene, dataset="GSE231761 (VSMC)",           lfc=tp$lfc_231761),
  data.frame(gene=tp$gene, dataset="GSE198983 (isogenic VSMC)",  lfc=tp$lfc_198983),
  data.frame(gene=tp$gene, dataset="GSE165409 (mouse aorta)",    lfc=tp$lfc_mouse),
  data.frame(gene=tp$gene, dataset="META (random effects)",      lfc=tp$meta_lfc))
fr$gene <- factor(fr$gene, levels = rev(tp$gene))
fr$dataset <- factor(fr$dataset, levels=c("GSE306264 (fibroblast)","GSE231761 (VSMC)",
                                          "GSE198983 (isogenic VSMC)","GSE165409 (mouse aorta)",
                                          "META (random effects)"))
lab <- data.frame(gene=factor(tp$gene, levels=rev(tp$gene)),
                  txt=sprintf("I2=%.0f%%", tp$I2), x=max(abs(fr$lfc), na.rm=TRUE)*1.02)
PNG(file.path(MFIG, "forest_top_candidates.png"), 9.5, 7)
print(ggplot(fr, aes(lfc, gene, color=dataset, shape=dataset)) +
  geom_vline(xintercept=0, linetype=2, color="grey50") +
  geom_point(size=2.6, position=position_dodge(width=.68)) +
  geom_text(data=lab, aes(x=x, y=gene, label=txt), inherit.aes=FALSE, size=2.6, hjust=0, color="grey35") +
  scale_color_manual(values=c("#C44E52","#DD8452","#55A868","#8172B3","black"), name=NULL) +
  scale_shape_manual(values=c(16,17,15,18,4), name=NULL) +
  scale_x_continuous(expand=expansion(mult=c(.05,.16))) +
  labs(title="Effect size of the final candidate genes across all datasets",
       subtitle="k=3 meta FDR<0.05, all concordant, I2<75%, |log2FC|>=1, mouse same direction (n = 16)",
       x="log2 fold change (progeria / control)", y=NULL) +
  theme_bw(base_size=11) + theme(legend.position="bottom"))
dev.off()

# =============================================================================
# [5] UpSet (5개 집합)
# =============================================================================
cat("[5] UpSet\n")
d306 <- read.delim(file.path(ROOT,"GSE306264/results/DE_HGPS_vs_WT_patient.tsv"))
d231 <- read.delim(file.path(ROOT,"GSE231761/results/DE_HGPS_vs_Control_VSMC.tsv"))
mmz  <- read.delim(file.path(RES,"GSE165409_mouse_to_human.tsv"))
iso  <- read.delim(file.path(ROOT,"GSE306264/results/DE_isogenic_Dox6d_vs_noDox.tsv"))
clean <- function(x) { x <- as.character(x); unique(x[!is.na(x) & nzchar(trimws(x))]) }
lst <- list(
  `GSE306264 (fibroblast)`    = clean(d306$gene[d306$sig != "NS"]),
  `GSE231761 (VSMC)`          = clean(d231$gene[d231$sig != "NS"]),
  `GSE198983 (isogenic VSMC)` = clean(de3$symbol[de3$sig != "NS"]),
  `GSE165409 (mouse)`         = clean(mmz$gene[!is.na(mmz$padj) & mmz$padj<.05 & abs(mmz$lfc)>=1]),
  `TetON acute`               = clean(iso$gene[!is.na(iso$padj) & iso$padj<.05 & abs(iso$log2FoldChange)>=1]))
cat("   집합 크기:"); print(vapply(lst, length, integer(1)))
cm <- make_comb_mat(lst); cm <- cm[comb_size(cm) >= 5]
PNG(file.path(MFIG, "upset_dataset_overlap.png"), 11, 5.5)
draw(UpSet(cm, comb_order = order(-comb_size(cm)),
           top_annotation = upset_top_annotation(cm, add_numbers = TRUE),
           right_annotation = upset_right_annotation(cm, add_numbers = TRUE),
           column_title = "Overlap of significant DEGs (FDR<0.05, |log2FC|>=1)"))
dev.off()

# =============================================================================
# [6] 데이터셋 쌍별 일치도
# =============================================================================
cat("[6] 쌍별 일치도\n")
sub <- M[M$meta_padj < .05, ]
mkp <- function(x, y, xl, yl) {
  r <- cor(sub[[x]], sub[[y]], method="spearman")
  ggplot(sub, aes(.data[[x]], .data[[y]], color = all_same_dir)) +
    geom_point(alpha=.55, size=1.1) +
    geom_hline(yintercept=0, linetype=2) + geom_vline(xintercept=0, linetype=2) +
    geom_abline(slope=1, intercept=0, color="grey55") +
    scale_color_manual(values=c(`TRUE`="#C44E52",`FALSE`="grey75"),
                       labels=c("discordant","all concordant"), name=NULL) +
    labs(subtitle=sprintf("Spearman = %.3f", r), x=xl, y=yl) +
    theme_bw(base_size=9) + theme(legend.position="none")
}
PNG(file.path(MFIG, "concordance_pairwise.png"), 11, 4)
grid.arrange(
  mkp("lfc_306264","lfc_231761","GSE306264 (fibroblast)","GSE231761 (VSMC)"),
  mkp("lfc_306264","lfc_198983","GSE306264 (fibroblast)","GSE198983 (isogenic)"),
  mkp("lfc_231761","lfc_198983","GSE231761 (VSMC)","GSE198983 (isogenic)"),
  ncol=3, top = grid::textGrob("Pairwise concordance of the three human datasets (meta FDR<0.05 genes)",
                               gp = grid::gpar(fontsize=12)))
dev.off()

# =============================================================================
# [7] 이질성 I2
# =============================================================================
cat("[7] I2 분포\n")
PNG(file.path(MFIG, "heterogeneity_I2.png"), 8, 4.5)
h <- M[!is.na(M$I2), ]
h$grp <- ifelse(h$meta_padj < .05, "meta FDR<0.05", "not significant")
print(ggplot(h, aes(I2, fill = grp)) +
  geom_histogram(bins=40, position="identity", alpha=.65) +
  geom_vline(xintercept=75, linetype=2, color="#C44E52") +
  annotate("text", x=78, y=Inf, vjust=1.5, hjust=0, label="I2 = 75%\n(high heterogeneity)",
           size=3, color="#C44E52") +
  scale_fill_manual(values=c("#C44E52","grey75"), name=NULL) +
  labs(title="Between-dataset heterogeneity (random-effects I2)",
       subtitle=sprintf("median I2 = %.1f%%; %d of %d genes exceed 75%%",
                        median(M$I2, na.rm=TRUE), sum(M$I2>75, na.rm=TRUE), nrow(M)),
       x="I2 (%)", y="number of genes") + theme_bw(base_size=11))
dev.off()

# =============================================================================
# [8] 데이터셋 요약 (5개)
# =============================================================================
cat("[8] 데이터셋 요약\n")
# [release] GSE137083 (screened, rejected; not part of the published analysis) is not
#           distributed with this repository, so it is drawn only if its DE table exists.
f137 <- file.path(ROOT,"GSE137083/results/DE_HGcon_vs_NFcon_all.tsv")
has137 <- file.exists(f137)
if (has137) d137 <- read.delim(f137)
smry <- data.frame(
  dataset = factor(c("GSE137083\n(rejected)","GSE306264","GSE231761","GSE198983","GSE165409\n(mouse)"),
                   levels=c("GSE137083\n(rejected)","GSE306264","GSE231761","GSE198983","GSE165409\n(mouse)")),
  used = c("rejected","meta","meta","meta","validation"),
  n_donor = c(2, 30, 4, 8, 12),
  n_deg = c(if (has137) sum(d137$sig!="NS") else NA, sum(d306$sig!="NS"), sum(d231$sig!="NS",na.rm=TRUE),
            sum(de3$sig!="NS",na.rm=TRUE),
            sum(!is.na(mmz$padj)&mmz$padj<.05&abs(mmz$lfc)>=1)),
  n_big = c(if (has137) sum(abs(d137$log2FoldChange)>5,na.rm=TRUE) else NA, sum(abs(d306$log2FoldChange)>5,na.rm=TRUE),
            sum(abs(d231$log2FoldChange)>5,na.rm=TRUE), sum(abs(de3$log2FoldChange)>5,na.rm=TRUE),
            sum(abs(mmz$lfc)>5,na.rm=TRUE)))
if (!has137) { smry <- smry[-1, ]; smry$dataset <- droplevels(smry$dataset) }
sm <- melt(smry, id.vars=c("dataset","used"))
levels(sm$variable) <- c("independent donors / batches / animals","significant DEGs",
                         "genes with |log2FC| > 5\n(artifact indicator)")
PNG(file.path(MFIG, "dataset_summary_comparison.png"), 11, 4.4)
print(ggplot(sm, aes(dataset, value, fill=used)) +
  geom_col(width=.68) + geom_text(aes(label=value), vjust=-0.35, size=3) +
  facet_wrap(~variable, scales="free_y") +
  scale_fill_manual(values=c(rejected="grey60", meta="#C44E52", validation="#55A868"), name=NULL) +
  scale_y_continuous(expand=expansion(mult=c(0,.18))) +
  labs(title="Dataset comparison across the progeria meta-analysis",
       subtitle=if (has137) "GSE137083: only 2 donors yet 4,081 DEGs and 540 genes with |log2FC|>5 = cell-line identity artifact" else NULL,
       x=NULL, y=NULL) + theme_bw(base_size=10) + theme(legend.position="bottom"))
dev.off()

# 마우스 보존성 갱신
PNG(file.path(MFIG, "mouse_conservation.png"), 7, 6)
mv <- M[!is.na(M$lfc_mouse) & M$meta_padj < .05, ]
mv$grp <- ifelse(mv$mouse_same_dir, "conserved (mouse sig + same dir)",
          ifelse(mv$mouse_dir_only, "same direction only", "not conserved"))
print(ggplot(mv, aes(meta_lfc, lfc_mouse, color=grp)) + geom_point(alpha=.55, size=1.3) +
  geom_hline(yintercept=0, linetype=2) + geom_vline(xintercept=0, linetype=2) +
  scale_color_manual(values=c(`conserved (mouse sig + same dir)`="#C44E52",
                              `same direction only`="#DD8452", `not conserved`="grey75"), name=NULL) +
  geom_text_repel(data=mv[mv$mouse_same_dir,], aes(label=gene), size=2.6, max.overlaps=40, show.legend=FALSE) +
  labs(title="Cross-species conservation of the k=3 human meta-signature",
       subtitle=sprintf("meta FDR<0.05 with 1:1 mouse ortholog (n=%d); Spearman = %.3f",
                        nrow(mv), cor(mv$meta_lfc, mv$lfc_mouse, method="spearman")),
       x="human meta log2FC (random effects)", y="mouse log2FC (GSE165409 aorta)") +
  theme_bw(base_size=11) + theme(legend.position="bottom"))
dev.off()

cat("\n완료\n")
