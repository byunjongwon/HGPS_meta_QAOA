#!/usr/bin/env Rscript
# =============================================================================
# 01번 수정판.
#
# 문제: 01번은 lfcShrink(apeglm) 결과의 log2FC / lfcSE 를 메타분석에 넣었다.
#       shrunken 추정치는 이미 0 쪽으로 당겨진 값이고 그 SE 는 사후분포 SD 라서,
#       이를 다시 역분산 가중 결합하면 '이중 축소'가 일어나 검정력이 떨어진다.
#       DESeq2 의 권고대로 메타분석 입력은 results() 의 MLE(비축소) log2FC + lfcSE 여야 한다.
#
# 이 스크립트는 저장된 dds 객체에서 비축소 추정치를 다시 뽑아 메타분석을 수행한다.
# (shrunken 값은 표시/랭킹용으로 계속 보존한다)
# =============================================================================
suppressPackageStartupMessages({
  library(DESeq2); library(metafor); library(ggplot2)
})

# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
.this <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
DEPOSIT <- Sys.getenv("PROGERIA_DEPOSIT", normalizePath(file.path(dirname(.this), "..", "..")))
ROOT <- file.path(DEPOSIT, "02_per_dataset_DE")  # [deposit] GSE* dataset dirs (was ~/car_project/progeria_meta)
RES  <- file.path(DEPOSIT, "03_meta_analysis", "results")  # [deposit]
FIG  <- file.path(RES, "figures"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
set.seed(1)

pull_mle <- function(rds, contrast, key = "gene") {
  o <- readRDS(rds)
  r <- results(o$dds, contrast = contrast, alpha = 0.05)   # MLE, 비축소
  d <- as.data.frame(r)
  d$id <- rownames(d)
  d <- d[!is.na(d$padj) & is.finite(d$lfcSE) & d$lfcSE > 0, ]
  d
}

cat("############ 비축소(MLE) 추정치 추출 ############\n")

# --- GSE306264 : gene symbol 이 rowname
d1 <- pull_mle(file.path(ROOT, "GSE306264/results/GSE306264_objects.rds"),
               c("genotype", "HGPS", "WT"))
d1$gene <- d1$id
d1 <- d1[!duplicated(d1$gene), ]
cat(sprintf("GSE306264 : %d 유전자 (유의 %d)\n", nrow(d1), sum(d1$padj < .05)))

# --- GSE231761 : Ensembl 이 rowname -> 저장된 res 에서 symbol 매핑을 가져온다
o2 <- readRDS(file.path(ROOT, "GSE231761/results/GSE231761_objects.rds"))
r2 <- results(o2$dds, contrast = c("condition", "HGPS", "Control"), alpha = 0.05)
d2 <- as.data.frame(r2); d2$ensembl <- rownames(d2)
map <- o2$res[, c("ensembl", "gene")]
d2 <- merge(d2, map, by = "ensembl")
d2 <- d2[!is.na(d2$padj) & is.finite(d2$lfcSE) & d2$lfcSE > 0 & !is.na(d2$gene) & d2$gene != "", ]
d2 <- d2[order(d2$padj), ]; d2 <- d2[!duplicated(d2$gene), ]
cat(sprintf("GSE231761 : %d 유전자 (유의 %d)\n", nrow(d2), sum(d2$padj < .05)))

# --- 마우스 (01번에서 만든 ortholog 매핑본; 이쪽은 원래 비축소 결과였다)
mm <- read.delim(file.path(RES, "GSE165409_mouse_to_human.tsv"))
iso <- read.delim(file.path(ROOT, "GSE306264/results/DE_isogenic_Dox6d_vs_noDox.tsv"))
iso <- iso[!is.na(iso$padj), c("gene", "log2FoldChange", "padj")]
names(iso) <- c("gene", "iso_lfc", "iso_padj"); iso <- iso[!duplicated(iso$gene), ]

# =============================================================================
# 메타분석 (사람 k=2, 고정효과)
# =============================================================================
cat("\n############ 메타분석 (비축소 입력) ############\n")
common <- intersect(d1$gene, d2$gene)
a <- d1[match(common, d1$gene), ]; b <- d2[match(common, d2$gene), ]
cat(sprintf("공통 유전자: %d\n", length(common)))

wa <- 1 / a$lfcSE^2; wb <- 1 / b$lfcSE^2
meta_lfc <- (wa * a$log2FoldChange + wb * b$log2FoldChange) / (wa + wb)
meta_se  <- sqrt(1 / (wa + wb))
meta_z   <- meta_lfc / meta_se
meta_p   <- 2 * pnorm(-abs(meta_z))
# Cochran Q -> I2
Q  <- wa * (a$log2FoldChange - meta_lfc)^2 + wb * (b$log2FoldChange - meta_lfc)^2
I2 <- pmax(0, 100 * (Q - 1) / Q)

M <- data.frame(
  gene = common, meta_lfc = meta_lfc, meta_se = meta_se, meta_p = meta_p,
  meta_padj = p.adjust(meta_p, "BH"), I2 = I2,
  lfc_306264 = a$log2FoldChange, padj_306264 = a$padj,
  lfc_231761 = b$log2FoldChange, padj_231761 = b$padj,
  same_dir = sign(a$log2FoldChange) == sign(b$log2FoldChange),
  stringsAsFactors = FALSE
)

cat(sprintf("\n[비교] 01번(shrunken 입력) 메타 FDR<0.05 : 90\n"))
cat(sprintf("[수정] 02번(비축소 입력) 메타 FDR<0.05   : %d\n", sum(M$meta_padj < .05)))

# 마우스 / isogenic 라벨
M <- merge(M, mm[, c("gene", "lfc", "padj")], by = "gene", all.x = TRUE)
names(M)[names(M) == "lfc"]  <- "lfc_mouse"
names(M)[names(M) == "padj"] <- "padj_mouse"
M <- merge(M, iso, by = "gene", all.x = TRUE)
M$mouse_sig       <- !is.na(M$padj_mouse) & M$padj_mouse < .05
M$mouse_same_dir  <- M$mouse_sig & (sign(M$lfc_mouse) == sign(M$meta_lfc))
M$mouse_dir_only  <- !is.na(M$lfc_mouse) & (sign(M$lfc_mouse) == sign(M$meta_lfc))
M$iso_sig         <- !is.na(M$iso_padj) & M$iso_padj < .05
M$chronicity      <- ifelse(!M$iso_sig, "chronic_only",
                     ifelse(sign(M$iso_lfc) == sign(M$meta_lfc), "acute_and_chronic", "acute_opposite"))
M <- M[order(M$meta_padj, -abs(M$meta_lfc)), ]
write.table(M, file.path(RES, "meta_k2_unshrunken.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# =============================================================================
# 필터 깔때기
# =============================================================================
cat("\n############ 필터 깔때기 ############\n")
base <- M$meta_padj < .05 & M$same_dir
steps <- list(
  `공통 유전자`                              = rep(TRUE, nrow(M)),
  `메타 FDR<0.05`                            = M$meta_padj < .05,
  `+ 사람 2개 방향 일치`                     = base,
  `+ |meta LFC|>=0.5`                        = base & abs(M$meta_lfc) >= .5,
  `+ |meta LFC|>=1`                          = base & abs(M$meta_lfc) >= 1,
  `+ 마우스 방향만 일치 (느슨)`              = base & abs(M$meta_lfc) >= 1 & M$mouse_dir_only,
  `+ 마우스 유의 + 방향 일치 (엄격)`         = base & abs(M$meta_lfc) >= 1 & M$mouse_same_dir
)
for (k in names(steps)) cat(sprintf("  %-38s %6d\n", k, sum(steps[[k]], na.rm = TRUE)))

cat("\n--- QUBO 규모 시나리오 ---\n")
scen <- list(
  `A: 사람만, |LFC|>=1`                 = base & abs(M$meta_lfc) >= 1,
  `B: 사람 + 마우스 방향 일치`          = base & abs(M$meta_lfc) >= 1 & M$mouse_dir_only,
  `C: 사람 + 마우스 유의·방향 일치`     = base & abs(M$meta_lfc) >= 1 & M$mouse_same_dir,
  `D: C + |LFC|>=0.5 로 완화`           = base & abs(M$meta_lfc) >= .5 & M$mouse_same_dir
)
for (k in names(scen)) {
  n <- sum(scen[[k]], na.rm = TRUE)
  hw <- if (n == 0) "-" else if (n <= 20) "Garnet(20) ✅ / Emerald ✅" else
        if (n <= 54) "Emerald(54) ✅" else sprintf("상위 54개로 축소 필요")
  ex <- if (n <= 28) sprintf("완전탐색 가능 (2^%d)", n) else "완전탐색 불가"
  cat(sprintf("  %-34s N=%3d   %-28s %s\n", k, n, hw, ex))
}

# 최종 후보 (시나리오 B: 20~54 목표)
final <- M[which(scen[["B: 사람 + 마우스 방향 일치"]]), ]
final <- final[order(-abs(final$meta_lfc)), ]
cat(sprintf("\n=== 시나리오 B 후보 %d 유전자 ===\n", nrow(final)))
print(head(final[, c("gene","meta_lfc","meta_padj","I2","lfc_306264","lfc_231761",
                     "lfc_mouse","mouse_sig","chronicity")], 60), digits = 3, row.names = FALSE)
write.table(final, file.path(RES, "candidates_scenarioB.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# 일치도 그림 (영문 라벨)
png(file.path(FIG, "concordance_306264_vs_231761.png"), width = 6, height = 5.5, units = "in", res = 300, type = "cairo")
sub <- M[M$meta_padj < .05, ]
print(ggplot(sub, aes(lfc_306264, lfc_231761, color = same_dir)) +
  geom_point(alpha = .45, size = 1.2) +
  geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 0, linetype = 2) +
  geom_abline(slope = 1, intercept = 0, color = "grey50") +
  scale_color_manual(values = c(`TRUE` = "#C44E52", `FALSE` = "grey70"),
                     labels = c("discordant", "concordant"), name = NULL) +
  labs(title = "Concordance between the two human datasets",
       subtitle = sprintf("meta FDR<0.05 (n=%d), Spearman = %.3f", nrow(sub),
                          cor(sub$lfc_306264, sub$lfc_231761, method = "spearman")),
       x = "log2FC  GSE306264 (fibroblast, 20 vs 10)",
       y = "log2FC  GSE231761 (VSMC, 2 vs 2)") +
  theme_bw(base_size = 11))
dev.off()
cat(sprintf("\n메타 유의 유전자에서 두 데이터셋 Spearman: %.3f\n",
            cor(sub$lfc_306264, sub$lfc_231761, method = "spearman")))
cat(sprintf("\n저장 완료 -> %s\n", RES))
