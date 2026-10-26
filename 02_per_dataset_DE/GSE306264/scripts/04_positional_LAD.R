#!/usr/bin/env Rscript
# =============================================================================
# GSE306264 : 위치 기반(genomic position) 분석
#
# HGPS 의 핵심 병태는 lamina-associated domain(LAD) 붕괴다. featureCounts 원본에
# Chr/Start/End 좌표가 들어 있으므로, 원시 데이터 없이도 다음을 검정할 수 있다.
#   (A) 염색체별 DEG 과대표현
#   (B) DEG 의 게놈상 군집화 (10 Mb 슬라이딩 윈도우 + 순열검정)
#   (C) LAD 대리지표(국소 유전자 밀도)와 log2FC 의 관계
#       LAD 는 유전자 희소·전사 억제 영역이므로, LAD 붕괴가 일어나면
#       '유전자 희소 영역'에서 발현 증가가 상대적으로 커야 한다.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(matrixStats) })

# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
.this <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
DEPOSIT <- Sys.getenv("PROGERIA_DEPOSIT", normalizePath(file.path(dirname(.this), "..", "..", "..")))
BASE <- file.path(DEPOSIT, "02_per_dataset_DE", "GSE306264")  # [deposit] was ~/car_project/progeria_meta/GSE306264
RES  <- file.path(BASE, "results"); FIG <- file.path(RES, "figures")
set.seed(1)

# ------------------------------------------------------------- 좌표 로드 ----
raw <- read.delim(gzfile(file.path(BASE, "raw_data", "GSE306264_featureCounts_all_gene_name.txt.gz")),
                  nrows = -1, check.names = FALSE)[, 1:6]
pos <- data.frame(
  gene  = raw$Geneid,
  chr   = vapply(strsplit(as.character(raw$Chr),   ";"), function(x) x[1], character(1)),
  start = vapply(strsplit(as.character(raw$Start), ";"), function(x) as.numeric(x[1]), numeric(1)),
  end   = vapply(strsplit(as.character(raw$End),   ";"), function(x) max(as.numeric(x)), numeric(1)),
  len   = raw$Length, stringsAsFactors = FALSE
)
main <- paste0("chr", c(1:22, "X"))
pos <- pos[pos$chr %in% main & !duplicated(pos$gene), ]
pos$mid <- (pos$start + pos$end) / 2
cat(sprintf("좌표 확보 유전자: %d (주요 염색체)\n", nrow(pos)))

de <- read.delim(file.path(RES, "DE_HGPS_vs_WT_patient.tsv"))
d  <- merge(de[, c("gene", "log2FoldChange", "padj", "sig", "baseMean")], pos, by = "gene")
d  <- d[!is.na(d$padj), ]
d$chr <- factor(d$chr, levels = main)
cat(sprintf("검정 유전자 중 좌표 있는 것: %d (UP %d / DOWN %d)\n",
            nrow(d), sum(d$sig == "UP"), sum(d$sig == "DOWN")))

# =============================================================================
# (A) 염색체별 DEG 과대표현
# =============================================================================
cat("\n############ (A) 염색체별 DEG 과대표현 ############\n")
tab <- do.call(rbind, lapply(levels(d$chr), function(cc) {
  ic <- d$chr == cc
  fu <- fisher.test(matrix(c(sum(ic & d$sig=="UP"),   sum(ic & d$sig!="UP"),
                             sum(!ic & d$sig=="UP"),  sum(!ic & d$sig!="UP")), 2))
  fd <- fisher.test(matrix(c(sum(ic & d$sig=="DOWN"), sum(ic & d$sig!="DOWN"),
                             sum(!ic & d$sig=="DOWN"),sum(!ic & d$sig!="DOWN")), 2))
  data.frame(chr = cc, n_gene = sum(ic), n_up = sum(ic & d$sig=="UP"),
             n_down = sum(ic & d$sig=="DOWN"),
             OR_up = fu$estimate, p_up = fu$p.value,
             OR_down = fd$estimate, p_down = fd$p.value,
             mean_lfc = mean(d$log2FoldChange[ic]))
}))
tab$padj_up <- p.adjust(tab$p_up, "BH"); tab$padj_down <- p.adjust(tab$p_down, "BH")
print(tab[order(tab$p_up), c("chr","n_gene","n_up","OR_up","padj_up","n_down","OR_down","padj_down","mean_lfc")],
      digits = 3, row.names = FALSE)
write.table(tab, file.path(RES, "positional_chromosome_enrichment.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 염색체별 평균 log2FC 가 0과 다른가 (전역 편향)
gm <- mean(d$log2FoldChange)
cat(sprintf("\n전체 평균 log2FC = %.4f\n", gm))
chrp <- sapply(levels(d$chr), function(cc)
  wilcox.test(d$log2FoldChange[d$chr == cc], d$log2FoldChange[d$chr != cc])$p.value)
sig_chr <- names(which(p.adjust(chrp, "BH") < .05))
cat("평균 log2FC 가 유의하게 다른 염색체:", ifelse(length(sig_chr), paste(sig_chr, collapse=", "), "없음"), "\n")

# =============================================================================
# (B) 게놈상 군집화 — 10 Mb 슬라이딩 윈도우
# =============================================================================
cat("\n############ (B) DEG 게놈 군집화 (10Mb 윈도우) ############\n")
W <- 10e6; STEP <- 5e6
win <- do.call(rbind, lapply(levels(d$chr), function(cc) {
  s <- d[d$chr == cc, ]; if (nrow(s) < 20) return(NULL)
  brks <- seq(0, max(s$mid) + W, by = STEP)
  do.call(rbind, lapply(brks, function(b) {
    ii <- s$mid >= b & s$mid < b + W
    if (sum(ii) < 10) return(NULL)
    data.frame(chr = cc, start = b, n = sum(ii),
               n_up = sum(s$sig[ii] == "UP"), mean_lfc = mean(s$log2FoldChange[ii]))
  }))
}))
win$frac_up <- win$n_up / win$n
cat(sprintf("분석 윈도우 %d개 (유전자 >=10)\n", nrow(win)))

# 순열검정: 유전자 위치를 섞었을 때 이만큼 높은 UP 비율 윈도우가 나오는가
obs <- max(win$frac_up)
perm <- replicate(1000, {
  ds <- d; ds$sig <- sample(ds$sig)
  mx <- 0
  for (cc in levels(ds$chr)) {
    s <- ds[ds$chr == cc, ]; if (nrow(s) < 20) next
    brks <- seq(0, max(s$mid) + W, by = STEP)
    for (b in brks) {
      ii <- s$mid >= b & s$mid < b + W
      if (sum(ii) < 10) next
      mx <- max(mx, sum(s$sig[ii] == "UP") / sum(ii))
    }
  }
  mx
})
cat(sprintf("관측 최대 UP비율 = %.3f, 순열 최대 UP비율 평균 = %.3f (95%%분위 %.3f)\n",
            obs, mean(perm), quantile(perm, .95)))
cat(sprintf("순열 p = %.4f  %s\n", mean(perm >= obs),
            ifelse(mean(perm >= obs) < .05, "<- DEG 가 게놈상 군집을 이룸", "<- 군집화 근거 약함")))

hot <- win[order(-win$frac_up), ][1:15, ]
cat("\n=== UP 유전자가 밀집한 상위 15개 구간 ===\n")
print(data.frame(chr = hot$chr, Mb = sprintf("%.0f-%.0f", hot$start/1e6, (hot$start+W)/1e6),
                 n = hot$n, n_up = hot$n_up, frac_up = round(hot$frac_up, 3),
                 mean_lfc = round(hot$mean_lfc, 2)), row.names = FALSE)
write.table(win, file.path(RES, "positional_windows.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

png(file.path(FIG, "positional_manhattan_lfc.png"), width = 12, height = 4.5, units = "in", res = 300, type = "cairo")
print(ggplot(d, aes(mid/1e6, log2FoldChange, color = sig)) +
  geom_point(size = .5, alpha = .5) +
  facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
  scale_color_manual(values = c(UP="#C44E52", DOWN="#4C72B0", NS="grey85")) +
  labs(title = "Chromosomal distribution of expression change (HGPS vs control)",
       x = "position (Mb)", y = "log2FC") +
  theme_bw(base_size = 8) +
  theme(panel.spacing = unit(0.02,"lines"), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), strip.text = element_text(size = 5)))
dev.off()

# =============================================================================
# (C) LAD 대리지표 : 국소 유전자 밀도
# =============================================================================
cat("\n############ (C) LAD 대리지표(국소 유전자 밀도) 와 발현 변화 ############\n")
# 각 유전자 주변 ±2Mb 안의 (좌표를 가진) 유전자 수
dens <- vapply(seq_len(nrow(d)), function(i) {
  s <- d[d$chr == d$chr[i], ]
  sum(abs(s$mid - d$mid[i]) <= 2e6) - 1
}, numeric(1))
d$gene_density <- dens
d$dens_decile  <- cut(rank(d$gene_density, ties.method = "first"),
                      breaks = quantile(rank(d$gene_density, ties.method="first"), probs = seq(0,1,.1)),
                      include.lowest = TRUE, labels = 1:10)

ds <- aggregate(cbind(log2FoldChange) ~ dens_decile, d, mean)
ds$n      <- as.numeric(table(d$dens_decile))
ds$frac_up   <- as.numeric(tapply(d$sig == "UP",   d$dens_decile, mean))
ds$frac_down <- as.numeric(tapply(d$sig == "DOWN", d$dens_decile, mean))
ds$med_density <- as.numeric(tapply(d$gene_density, d$dens_decile, median))
cat("\n=== 유전자 밀도 십분위별 (1=가장 희소=LAD 유사, 10=가장 조밀) ===\n")
print(ds, digits = 3, row.names = FALSE)

ct <- cor.test(d$gene_density, d$log2FoldChange, method = "spearman")
cat(sprintf("\n유전자밀도 vs log2FC Spearman rho = %.4f, p = %.3g\n", ct$estimate, ct$p.value))
lo <- d$dens_decile %in% c("1","2"); hi <- d$dens_decile %in% c("9","10")
wt <- wilcox.test(d$log2FoldChange[lo], d$log2FoldChange[hi])
cat(sprintf("희소영역(1-2분위) 평균 log2FC = %.4f  vs  조밀영역(9-10분위) = %.4f,  p = %.3g\n",
            mean(d$log2FoldChange[lo]), mean(d$log2FoldChange[hi]), wt$p.value))
cat(sprintf("UP 비율:  희소 %.3f  vs  조밀 %.3f\n", mean(d$sig[lo]=="UP"), mean(d$sig[hi]=="UP")))
cat("  해석: 희소(LAD 유사) 영역에서 발현 증가가 더 크면 LAD 붕괴 가설과 부합한다\n")

write.table(d, file.path(RES, "positional_gene_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ds, file.path(RES, "positional_density_deciles.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

png(file.path(FIG, "LAD_proxy_density_vs_lfc.png"), width = 7, height = 5, units = "in", res = 300, type = "cairo")
print(ggplot(d, aes(dens_decile, log2FoldChange)) +
  geom_boxplot(outlier.size = .3, fill = "#B0C4DE") +
  geom_hline(yintercept = 0, linetype = 2, color = "red") +
  labs(title = "Expression change vs local gene density (LAD proxy)",
       subtitle = sprintf("Spearman rho = %.3f (p = %.2g);  decile 1 = gene-poor (LAD-like)",
                          ct$estimate, ct$p.value),
       x = "local gene density decile (±2 Mb)", y = "log2FC (HGPS / control)") +
  coord_cartesian(ylim = quantile(d$log2FoldChange, c(.005, .995))) +
  theme_bw(base_size = 11))
dev.off()
cat(sprintf("\n저장 완료 -> %s\n", RES))
