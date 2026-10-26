#!/usr/bin/env Rscript
# =============================================================================
# 사람 3개 데이터셋 메타분석 (k = 3)  — 최종본
#
#   GSE306264  섬유아세포   HGPS 20라인 vs 정상 10라인      (~ genotype, 라인 collapse)
#   GSE231761  VSMC         HGPS 2라인 vs 대조 2라인         (~ passage + condition)
#   GSE198983  VSMC         isogenic Flag-progerin vs Flag   (~ batch + group)   <- 신규
#
#  * GSE198983 은 '4 replicates' 가 기술반복이 아니라 **균형 2x2 배치구조** 였다.
#    (배치 A = rep 1,4 / 배치 B = rep 2,3, 각 배치에 Control 2 + Progerin 2)
#    따라서 collapse 가 아니라 blocking(~ batch + group) 이 정답이고, 그 결과를 쓴다.
#    - blocked  644 DEG   <- 채택
#    - simple   431
#    - collapse  65 (반복을 버려 과소검출)
#
#  * 입력은 전부 비축소(MLE) log2FC + lfcSE. k=3 이므로 랜덤효과(DL) 로 tau^2/I2 를 추정한다.
#  * GSE231761 의 CNV 위양성 4종은 입력에서 제외한다.
# =============================================================================
suppressPackageStartupMessages({
  library(DESeq2); library(metafor); library(ggplot2); library(ggrepel)
})
# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
.this <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
DEPOSIT <- Sys.getenv("PROGERIA_DEPOSIT", normalizePath(file.path(dirname(.this), "..", "..")))
ROOT <- file.path(DEPOSIT, "02_per_dataset_DE")  # [deposit] GSE* dataset dirs (was ~/car_project/progeria_meta)
RES  <- file.path(DEPOSIT, "03_meta_analysis", "results"); FIG <- file.path(RES, "figures")  # [deposit]
set.seed(1)

# ------------------------------------------------------------------ 입력 ----
cat("############ 입력 (비축소 MLE) ############\n")

o1 <- readRDS(file.path(ROOT, "GSE306264/results/GSE306264_objects.rds"))
d1 <- as.data.frame(results(o1$dds, contrast = c("genotype","HGPS","WT"), alpha = .05))
d1$gene <- rownames(d1)
d1 <- d1[!is.na(d1$padj) & is.finite(d1$lfcSE) & d1$lfcSE > 0, ]
d1 <- d1[!duplicated(d1$gene), ]
cat(sprintf("GSE306264 : %5d 유전자 (유의 %d)\n", nrow(d1), sum(d1$padj < .05)))

o2 <- readRDS(file.path(ROOT, "GSE231761/results/GSE231761_objects.rds"))
d2 <- as.data.frame(results(o2$dds, contrast = c("condition","HGPS","Control"), alpha = .05))
d2$ensembl <- rownames(d2)
d2 <- merge(d2, unique(o2$res[, c("ensembl","gene")]), by = "ensembl")
d2 <- d2[!is.na(d2$padj) & is.finite(d2$lfcSE) & d2$lfcSE > 0 & !is.na(d2$gene) & d2$gene != "", ]
d2 <- d2[order(d2$padj), ]; d2 <- d2[!duplicated(d2$gene), ]
# CNV 위양성 제외
cnv <- read.delim(file.path(ROOT, "GSE231761/results/CNV_polymorphism_candidates.tsv"), check.names = FALSE)
bad <- cnv$gene[cnv$splits_by_condition %in% c(TRUE, "TRUE")]
cat(sprintf("GSE231761 : %5d 유전자 (유의 %d)  | CNV 위양성 제외: %s\n",
            nrow(d2), sum(d2$padj < .05), paste(bad, collapse = ", ")))
d2 <- d2[!d2$gene %in% bad, ]

o3 <- readRDS(file.path(DEPOSIT, "01_rnaseq_processing/GSE198983_human_AoSMC/results/09_deseq2/dds.rds"))  # [deposit]
cat("  GSE198983 설계: "); print(design(o3))
r3 <- results(o3, contrast = c("group","Progerin","Control"), alpha = .05)
d3 <- as.data.frame(r3); d3$gene_id <- rownames(d3)
ann3 <- read.delim(file.path(DEPOSIT, "01_rnaseq_processing/GSE198983_human_AoSMC/results/09_deseq2/DE_Progerin_vs_Control_blocked.tsv"))  # [deposit]
d3 <- merge(d3, ann3[, c("gene_id","symbol")], by = "gene_id")
names(d3)[names(d3) == "symbol"] <- "gene"
d3 <- d3[!is.na(d3$padj) & is.finite(d3$lfcSE) & d3$lfcSE > 0 & !is.na(d3$gene) & d3$gene != "", ]
d3 <- d3[order(d3$padj), ]; d3 <- d3[!duplicated(d3$gene), ]
cat(sprintf("GSE198983 : %5d 유전자 (유의 %d)\n", nrow(d3), sum(d3$padj < .05)))

# ------------------------------------------------------------ 메타분석 -----
common <- Reduce(intersect, list(d1$gene, d2$gene, d3$gene))
cat(sprintf("\n3개 데이터셋 공통 유전자: %d\n", length(common)))
a <- d1[match(common, d1$gene), ]; b <- d2[match(common, d2$gene), ]; c3 <- d3[match(common, d3$gene), ]

cat("랜덤효과(DL) 메타분석 실행 중...\n")
mt <- t(vapply(seq_along(common), function(i) {
  y <- c(a$log2FoldChange[i], b$log2FoldChange[i], c3$log2FoldChange[i])
  v <- c(a$lfcSE[i], b$lfcSE[i], c3$lfcSE[i])^2
  fe <- tryCatch(rma.uni(yi = y, vi = v, method = "FE"), error = function(e) NULL)
  re <- tryCatch(rma.uni(yi = y, vi = v, method = "DL"), error = function(e) NULL)
  if (is.null(fe) || is.null(re)) return(rep(NA_real_, 6))
  c(re$beta[1], re$se[1], re$pval, re$I2, re$tau2, fe$pval)
}, numeric(6)))
colnames(mt) <- c("meta_lfc","meta_se","meta_p","I2","tau2","meta_p_FE")

M <- data.frame(gene = common, mt, stringsAsFactors = FALSE)
M$meta_padj <- p.adjust(M$meta_p, "BH")
M$lfc_306264 <- a$log2FoldChange; M$padj_306264 <- a$padj
M$lfc_231761 <- b$log2FoldChange; M$padj_231761 <- b$padj
M$lfc_198983 <- c3$log2FoldChange; M$padj_198983 <- c3$padj
sg <- cbind(sign(a$log2FoldChange), sign(b$log2FoldChange), sign(c3$log2FoldChange))
M$n_same_dir <- apply(sg, 1, function(x) max(sum(x > 0), sum(x < 0)))
M$all_same_dir <- M$n_same_dir == 3

cat(sprintf("\n메타 FDR<0.05                     : %d\n", sum(M$meta_padj < .05, na.rm = TRUE)))
cat(sprintf("메타 FDR<0.05 & 3개 방향 일치     : %d\n", sum(M$meta_padj < .05 & M$all_same_dir, na.rm = TRUE)))
cat(sprintf("이질성 I2 중앙값 %.1f%%,  I2>75%% 인 유전자 %d개 (%.1f%%)\n",
            median(M$I2, na.rm = TRUE), sum(M$I2 > 75, na.rm = TRUE), 100*mean(M$I2 > 75, na.rm = TRUE)))

# 마우스 / isogenic 라벨
mm <- read.delim(file.path(RES, "GSE165409_mouse_to_human.tsv"))
M <- merge(M, mm[, c("gene","lfc","padj")], by = "gene", all.x = TRUE)
names(M)[names(M)=="lfc"] <- "lfc_mouse"; names(M)[names(M)=="padj"] <- "padj_mouse"
M$mouse_sig      <- !is.na(M$padj_mouse) & M$padj_mouse < .05
M$mouse_same_dir <- M$mouse_sig & sign(M$lfc_mouse) == sign(M$meta_lfc)
M$mouse_dir_only <- !is.na(M$lfc_mouse) & sign(M$lfc_mouse) == sign(M$meta_lfc)
iso <- read.delim(file.path(ROOT, "GSE306264/results/DE_isogenic_Dox6d_vs_noDox.tsv"))
iso <- iso[!is.na(iso$padj), c("gene","log2FoldChange","padj")]
names(iso) <- c("gene","iso_lfc","iso_padj"); iso <- iso[!duplicated(iso$gene), ]
M <- merge(M, iso, by = "gene", all.x = TRUE)
M$chronicity <- ifelse(is.na(M$iso_padj) | M$iso_padj >= .05, "chronic_only",
                ifelse(sign(M$iso_lfc) == sign(M$meta_lfc), "acute_and_chronic", "acute_opposite"))
M <- M[order(M$meta_padj, -abs(M$meta_lfc)), ]
write.table(M, file.path(RES, "meta_k3_final.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------------------------------------ 깔때기 -------
cat("\n############ 필터 깔때기 (k=3) ############\n")
base <- M$meta_padj < .05 & M$all_same_dir
steps <- list(
  `공통 유전자 (사람 3개)`             = rep(TRUE, nrow(M)),
  `메타 FDR<0.05`                      = M$meta_padj < .05,
  `+ 3개 데이터셋 방향 일치`           = base,
  `+ I2 < 75 (이질성 낮음)`            = base & M$I2 < 75,
  `+ |meta log2FC| >= 0.5`             = base & M$I2 < 75 & abs(M$meta_lfc) >= .5,
  `+ |meta log2FC| >= 1`               = base & M$I2 < 75 & abs(M$meta_lfc) >= 1,
  `+ 마우스 방향 일치`                 = base & M$I2 < 75 & abs(M$meta_lfc) >= 1 & M$mouse_dir_only,
  `+ 마우스 유의 & 방향 일치`          = base & M$I2 < 75 & abs(M$meta_lfc) >= 1 & M$mouse_same_dir)
fn <- data.frame(step = names(steps), n = vapply(steps, function(x) sum(x, na.rm = TRUE), integer(1)))
for (i in seq_len(nrow(fn))) cat(sprintf("  %-38s %6d\n", fn$step[i], fn$n[i]))

for (nm in c("+ |meta log2FC| >= 0.5","+ |meta log2FC| >= 1","+ 마우스 방향 일치","+ 마우스 유의 & 방향 일치")) {
  n <- sum(steps[[nm]], na.rm = TRUE)
  hw <- if (n == 0) "-" else if (n <= 20) "Garnet(20) OK / Emerald OK" else
        if (n <= 54) "Emerald(54) OK" else "상위 54개로 축소"
  cat(sprintf("    [QUBO] %-28s N=%4d  %s\n", nm, n, hw))
}

FINAL <- M[which(steps[["+ 마우스 방향 일치"]]), ]
FINAL <- FINAL[order(-abs(FINAL$meta_lfc)), ]
cat(sprintf("\n=== 최종 후보 (마우스 방향 일치) %d 유전자 — 상위 40 ===\n", nrow(FINAL)))
print(head(FINAL[, c("gene","meta_lfc","meta_padj","I2","lfc_306264","lfc_231761",
                     "lfc_198983","lfc_mouse","mouse_sig","chronicity")], 40),
      digits = 3, row.names = FALSE)
write.table(FINAL, file.path(RES, "candidates_k3_final.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

STRICT <- M[which(steps[["+ 마우스 유의 & 방향 일치"]]), ]
STRICT <- STRICT[order(-abs(STRICT$meta_lfc)), ]
write.table(STRICT, file.path(RES, "candidates_k3_strict.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n엄격 집합 (마우스 유의): %d 유전자\n", nrow(STRICT)))
if (nrow(STRICT)) print(STRICT[, c("gene","meta_lfc","meta_padj","I2","lfc_mouse")], digits = 3, row.names = FALSE)

# 데이터셋 쌍별 상관
cat("\n=== 데이터셋 간 log2FC 일치도 (메타 유의 유전자) ===\n")
sub <- M[M$meta_padj < .05, ]
prs <- list(c("lfc_306264","lfc_231761"), c("lfc_306264","lfc_198983"), c("lfc_231761","lfc_198983"))
for (p in prs) cat(sprintf("  %-24s vs %-12s Spearman = %.3f\n", p[1], p[2],
                           cor(sub[[p[1]]], sub[[p[2]]], method = "spearman")))
cat(sprintf("\n저장 완료 -> %s\n", RES))
