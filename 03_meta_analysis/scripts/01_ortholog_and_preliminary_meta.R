#!/usr/bin/env Rscript
# =============================================================================
# 조로증 메타분석 — 예비 단계 (GSE198983 도착 전)
#
#  [1] GSE165409 (마우스 대동맥, WT_24mo 제외본) 을 사람 유전자로 ortholog 매핑
#  [2] 사람 2개 데이터셋(GSE306264, GSE231761)으로 예비 메타분석 (k=2)
#  [3] 마우스 보존성까지 적용했을 때 최종 유전자가 몇 개 남는지 -> QUBO 큐빗 수 예측
#
#  주의: k=2 에서는 랜덤효과 tau^2 추정이 불안정하므로 고정효과(FE) 를 주로 쓰고,
#        RE 는 참고로만 계산한다. GSE198983 이 들어오면 k=3 RE 로 전환한다.
# =============================================================================
suppressPackageStartupMessages({
  library(metafor); library(babelgene); library(RobustRankAggreg); library(ggplot2)
})

# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
.this <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
DEPOSIT <- Sys.getenv("PROGERIA_DEPOSIT", normalizePath(file.path(dirname(.this), "..", "..")))
ROOT <- file.path(DEPOSIT, "02_per_dataset_DE")  # [deposit] GSE* dataset dirs (was ~/car_project/progeria_meta)
RES  <- file.path(DEPOSIT, "03_meta_analysis", "results")  # [deposit]
FIG  <- file.path(RES, "figures"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
set.seed(1)

# =============================================================================
# [1] 입력 로드
# =============================================================================
cat("############ [1] 입력 데이터 ############\n")

h1 <- read.delim(file.path(ROOT, "GSE306264/results/DE_HGPS_vs_WT_patient.tsv"))
h1 <- h1[!is.na(h1$padj), c("gene", "log2FoldChange", "lfcSE", "pvalue", "padj", "baseMean")]
names(h1) <- c("gene", "lfc", "se", "p", "padj", "baseMean")
h1 <- h1[!duplicated(h1$gene), ]
cat(sprintf("GSE306264 (섬유아세포, 20 vs 10) : %d 유전자, 유의(FDR<.05) %d\n",
            nrow(h1), sum(h1$padj < .05)))

h2 <- read.delim(file.path(ROOT, "GSE231761/results/DE_HGPS_vs_Control_VSMC.tsv"))
h2 <- h2[!is.na(h2$padj) & !is.na(h2$gene) & h2$gene != "", ]
h2 <- h2[order(h2$padj), c("gene", "log2FoldChange", "lfcSE", "pvalue", "padj", "baseMean")]
names(h2) <- c("gene", "lfc", "se", "p", "padj", "baseMean")
h2 <- h2[!duplicated(h2$gene), ]
cat(sprintf("GSE231761 (VSMC, 2 vs 2)         : %d 유전자, 유의(FDR<.05) %d\n",
            nrow(h2), sum(h2$padj < .05)))

iso <- read.delim(file.path(ROOT, "GSE306264/results/DE_isogenic_Dox6d_vs_noDox.tsv"))
iso <- iso[!is.na(iso$padj), c("gene", "log2FoldChange", "padj")]
names(iso) <- c("gene", "iso_lfc", "iso_padj")
iso <- iso[!duplicated(iso$gene), ]
cat(sprintf("isogenic TetON (급성 6일)         : %d 유전자, 유의 %d\n",
            nrow(iso), sum(iso$iso_padj < .05)))

# =============================================================================
# [2] GSE165409 마우스 -> 사람 ortholog 매핑
# =============================================================================
cat("\n\n############ [2] 마우스 ortholog 매핑 ############\n")
m <- read.delim(file.path(DEPOSIT, "01_rnaseq_processing/GSE165409_mouse_aorta/results/08_deseq2/DE_LMNA_vs_WTyng.tsv"))  # [deposit]
m <- m[!is.na(m$padj) & !is.na(m$symbol) & m$symbol != "", ]
m <- m[, c("symbol", "log2FoldChange", "lfcSE", "pvalue", "padj", "baseMean")]
names(m) <- c("mouse_symbol", "lfc", "se", "p", "padj", "baseMean")
m <- m[order(m$padj), ]; m <- m[!duplicated(m$mouse_symbol), ]
cat(sprintf("GSE165409 (마우스 대동맥, LMNA 2mo vs WT 2mo) : %d 유전자, 유의 %d\n",
            nrow(m), sum(m$padj < .05)))

# babelgene: 1:1 ortholog 만 사용 (모호한 1:多 는 제외)
orth <- orthologs(genes = m$mouse_symbol, species = "mouse", human = FALSE)
cat(sprintf("babelgene 매핑 후보: %d 행\n", nrow(orth)))
# 마우스/사람 양쪽에서 1:1 인 것만
tb_m <- table(orth$symbol); tb_h <- table(orth$human_symbol)
o11 <- orth[orth$symbol %in% names(tb_m)[tb_m == 1] &
            orth$human_symbol %in% names(tb_h)[tb_h == 1], c("symbol", "human_symbol")]
cat(sprintf("1:1 ortholog 만 남김: %d 쌍\n", nrow(o11)))

mm <- merge(m, o11, by.x = "mouse_symbol", by.y = "symbol")
mm <- mm[!duplicated(mm$human_symbol), ]
names(mm)[names(mm) == "human_symbol"] <- "gene"
cat(sprintf("마우스 DE 결과 중 사람 유전자로 매핑된 것: %d / %d (%.1f%%)\n",
            nrow(mm), nrow(m), 100 * nrow(mm) / nrow(m)))
cat(sprintf("  그중 마우스에서 유의(FDR<.05): %d\n", sum(mm$padj < .05)))
write.table(mm, file.path(RES, "GSE165409_mouse_to_human.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# =============================================================================
# [3] 예비 메타분석 (사람 k=2)
# =============================================================================
cat("\n\n############ [3] 예비 메타분석 (사람 k=2) ############\n")
common <- intersect(h1$gene, h2$gene)
cat(sprintf("두 사람 데이터셋 공통 유전자: %d\n", length(common)))

a <- h1[match(common, h1$gene), ]; b <- h2[match(common, h2$gene), ]
ok <- is.finite(a$lfc) & is.finite(a$se) & is.finite(b$lfc) & is.finite(b$se) &
      a$se > 0 & b$se > 0
common <- common[ok]; a <- a[ok, ]; b <- b[ok, ]
cat(sprintf("유효 SE 를 가진 유전자: %d\n", length(common)))

cat("메타분석 실행 중 (고정효과 + 랜덤효과)...\n")
meta <- t(vapply(seq_along(common), function(i) {
  y <- c(a$lfc[i], b$lfc[i]); v <- c(a$se[i], b$se[i])^2
  fe <- tryCatch(rma.uni(yi = y, vi = v, method = "FE"), error = function(e) NULL)
  re <- tryCatch(rma.uni(yi = y, vi = v, method = "DL"), error = function(e) NULL)
  if (is.null(fe)) return(c(NA, NA, NA, NA, NA))
  c(fe$beta[1], fe$se[1], fe$pval,
    if (is.null(re)) NA else re$I2,
    if (is.null(re)) NA else re$pval)
}, numeric(5)))
colnames(meta) <- c("meta_lfc", "meta_se", "meta_p_FE", "I2", "meta_p_RE")

M <- data.frame(gene = common, meta, stringsAsFactors = FALSE)
M$meta_padj <- p.adjust(M$meta_p_FE, "BH")
M$lfc_306264 <- a$lfc; M$padj_306264 <- a$padj
M$lfc_231761 <- b$lfc; M$padj_231761 <- b$padj
M$same_dir <- sign(a$lfc) == sign(b$lfc)

cat(sprintf("\n메타 FDR<0.05                       : %d\n", sum(M$meta_padj < .05, na.rm = TRUE)))
cat(sprintf("메타 FDR<0.05 & 방향 일치           : %d\n",
            sum(M$meta_padj < .05 & M$same_dir, na.rm = TRUE)))
cat(sprintf("메타 FDR<0.05 & 방향 일치 & |LFC|>=1: %d\n",
            sum(M$meta_padj < .05 & M$same_dir & abs(M$meta_lfc) >= 1, na.rm = TRUE)))
cat(sprintf("  이질성 I2 중앙값: %.1f%%\n", median(M$I2, na.rm = TRUE)))

# 두 데이터셋 각각에서 유의 + 방향 일치 (더 엄격)
M$both_sig <- M$padj_306264 < .05 & M$padj_231761 < .05 & M$same_dir
cat(sprintf("\n[엄격] 두 데이터셋 각각 FDR<0.05 + 방향 일치 : %d\n", sum(M$both_sig, na.rm = TRUE)))

# =============================================================================
# [4] 마우스 보존성 + isogenic 라벨 부착
# =============================================================================
cat("\n\n############ [4] 마우스 보존성 / isogenic 라벨 ############\n")
M <- merge(M, mm[, c("gene", "lfc", "padj")], by = "gene", all.x = TRUE, suffixes = c("", "_mouse"))
names(M)[names(M) == "lfc"] <- "lfc_mouse"; names(M)[names(M) == "padj"] <- "padj_mouse"
M <- merge(M, iso, by = "gene", all.x = TRUE)

M$mouse_tested    <- !is.na(M$padj_mouse)
M$mouse_sig       <- !is.na(M$padj_mouse) & M$padj_mouse < .05
M$mouse_same_dir  <- M$mouse_sig & (sign(M$lfc_mouse) == sign(M$meta_lfc))
M$iso_sig         <- !is.na(M$iso_padj) & M$iso_padj < .05
M$iso_same_dir    <- M$iso_sig & (sign(M$iso_lfc) == sign(M$meta_lfc))
M$chronicity      <- ifelse(!M$iso_sig, "chronic_only",
                     ifelse(M$iso_same_dir, "acute_and_chronic", "acute_opposite"))

cat(sprintf("메타 유전자 중 마우스에서도 검정된 것 : %d\n", sum(M$mouse_tested)))
cat(sprintf("                마우스에서도 유의     : %d\n", sum(M$mouse_sig)))
cat(sprintf("                마우스 방향까지 일치  : %d\n", sum(M$mouse_same_dir)))

M <- M[order(M$meta_padj, -abs(M$meta_lfc)), ]
write.table(M, file.path(RES, "preliminary_meta_k2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# =============================================================================
# [5] 필터 깔때기 -> QUBO 큐빗 수 예측
# =============================================================================
cat("\n\n############ [5] 필터 깔때기 (QUBO 규모 예측) ############\n")
f <- list(
  `공통 유전자 (사람 2개)`                        = rep(TRUE, nrow(M)),
  `+ 메타 FDR<0.05`                               = M$meta_padj < .05,
  `+ 두 데이터셋 방향 일치`                       = M$meta_padj < .05 & M$same_dir,
  `+ |메타 log2FC| >= 0.5`                        = M$meta_padj < .05 & M$same_dir & abs(M$meta_lfc) >= .5,
  `+ |메타 log2FC| >= 1`                          = M$meta_padj < .05 & M$same_dir & abs(M$meta_lfc) >= 1,
  `+ 마우스에서도 유의`                           = M$meta_padj < .05 & M$same_dir & abs(M$meta_lfc) >= 1 & M$mouse_sig,
  `+ 마우스 방향까지 일치 (최종)`                 = M$meta_padj < .05 & M$same_dir & abs(M$meta_lfc) >= 1 & M$mouse_same_dir
)
for (k in names(f)) cat(sprintf("  %-42s %5d\n", k, sum(f[[k]], na.rm = TRUE)))

final <- M[which(f[[length(f)]]), ]
cat(sprintf("\n=== 최종 후보 %d 유전자 ===\n", nrow(final)))
if (nrow(final) > 0) {
  print(head(final[order(-abs(final$meta_lfc)),
        c("gene","meta_lfc","meta_padj","I2","lfc_306264","lfc_231761","lfc_mouse","chronicity")],
        40), digits = 3, row.names = FALSE)
}
write.table(final, file.path(RES, "preliminary_final_candidates.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n--- QUBO 하드웨어 배정 예측 ---\n")
n <- nrow(final)
cat(sprintf("  최종 유전자 수 N = %d\n", n))
cat(sprintf("  IQM Garnet (20큐빗)  : %s\n", ifelse(n <= 20, "적합 ✅", "N 축소 필요")))
cat(sprintf("  IQM Emerald (54큐빗) : %s\n", ifelse(n <= 54, "적합 ✅", "상위 54개로 축소")))
cat(sprintf("  완전탐색 가능 여부   : %s (2^%d = %.3g)\n",
            ifelse(n <= 28, "가능 ✅", "상위 25개만"), min(n, 28), 2^min(n, 28)))

# 완화 시나리오 (마우스 조건을 뺐을 때)
n2 <- sum(f[["+ |메타 log2FC| >= 1"]], na.rm = TRUE)
n3 <- sum(f[["+ |메타 log2FC| >= 0.5"]], na.rm = TRUE)
cat(sprintf("\n  [참고] 마우스 조건 제외 시 : |LFC|>=1 -> %d,  |LFC|>=0.5 -> %d\n", n2, n3))
cat("         네트워크 QUBO 는 20~54개가 이상적이므로 이 범위에 맞게 임계값을 조정하면 된다\n")

# 스캐터: 두 데이터셋 log2FC 일치도
png(file.path(FIG, "concordance_306264_vs_231761.png"), width = 6, height = 5.5, units = "in", res = 300, type = "cairo")
sub <- M[M$meta_padj < .05, ]
print(ggplot(sub, aes(lfc_306264, lfc_231761, color = same_dir)) +
  geom_point(alpha = .4, size = 1) +
  geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 0, linetype = 2) +
  geom_abline(slope = 1, intercept = 0, color = "grey50") +
  scale_color_manual(values = c(`TRUE` = "#C44E52", `FALSE` = "grey70"),
                     labels = c("방향 불일치", "방향 일치"), name = NULL) +
  labs(title = "Concordance between the two human datasets",
       subtitle = sprintf("meta FDR<0.05 genes (n=%d), Spearman = %.3f", nrow(sub),
                          cor(sub$lfc_306264, sub$lfc_231761, method = "spearman")),
       x = "log2FC  GSE306264 (fibroblast)", y = "log2FC  GSE231761 (VSMC)") +
  theme_bw(base_size = 11))
dev.off()
cat(sprintf("\n두 사람 데이터셋 log2FC Spearman 상관 (메타 유의 유전자): %.3f\n",
            cor(sub$lfc_306264, sub$lfc_231761, method = "spearman")))

cat(sprintf("\n저장 완료 -> %s\n", RES))
