#!/usr/bin/env Rscript
# =============================================================================
# GSE306264 경로 분석
#   (A) 환자 비교 (HGPS 20라인 vs WT 10라인)  : GSEA Hallmark / GO:BP / Reactome / KEGG
#   (B) isogenic  (TetON progerin 6일)        : 동일 세트
#   (C) 급성 vs 만성 경로 수준 비교
#
# GSEA 입력은 비축소(MLE) log2FC 를 SE 로 나눈 t 유사 통계량(Wald stat)을 쓴다.
# (shrunken LFC 는 랭킹을 왜곡할 수 있어 메타분석과 동일한 원칙 적용)
# =============================================================================
suppressPackageStartupMessages({
  library(DESeq2); library(fgsea); library(msigdbr); library(ggplot2)
  library(clusterProfiler); library(ReactomePA); library(org.Hs.eg.db); library(AnnotationDbi)
})

# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
.this <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
DEPOSIT <- Sys.getenv("PROGERIA_DEPOSIT", normalizePath(file.path(dirname(.this), "..", "..", "..")))
BASE <- file.path(DEPOSIT, "02_per_dataset_DE", "GSE306264")  # [deposit] was ~/car_project/progeria_meta/GSE306264
RES  <- file.path(BASE, "results"); FIG <- file.path(RES, "figures")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
set.seed(1)

o <- readRDS(file.path(RES, "GSE306264_objects.rds"))

# ------------------------------------------------------- 랭킹 통계량 준비 ----
rank_from <- function(dds, contrast, coefname) {
  r <- as.data.frame(results(dds, contrast = contrast, alpha = 0.05))
  r <- r[!is.na(r$stat), ]
  s <- setNames(r$stat, rownames(r))
  sort(s[!duplicated(names(s))], decreasing = TRUE)
}
stat_pat <- rank_from(o$dds, c("genotype", "HGPS", "WT"))
cat(sprintf("환자 랭킹 유전자: %d\n", length(stat_pat)))

# isogenic 은 dds 를 저장하지 않았으므로 결과 파일에서 재구성 (Wald stat 대용: lfc/lfcSE)
iso <- read.delim(file.path(RES, "DE_isogenic_Dox6d_vs_noDox.tsv"))
iso <- iso[!is.na(iso$padj) & is.finite(iso$lfcSE) & iso$lfcSE > 0, ]
stat_iso <- setNames(iso$log2FoldChange / iso$lfcSE, iso$gene)
stat_iso <- sort(stat_iso[!duplicated(names(stat_iso))], decreasing = TRUE)
cat(sprintf("isogenic 랭킹 유전자: %d\n", length(stat_iso)))

# ------------------------------------------------------------ GSEA 실행 -----
msig <- function(coll, sub = NULL) {
  m <- if (is.null(sub)) msigdbr(species = "Homo sapiens", collection = coll)
       else msigdbr(species = "Homo sapiens", collection = coll, subcollection = sub)
  s <- split(m$gene_symbol, m$gs_name)
  s[lengths(s) >= 15 & lengths(s) <= 500]
}
SETS <- list(Hallmark = msig("H"), GO_BP = msig("C5", "GO:BP"),
             Reactome = msig("C2", "CP:REACTOME"), KEGG = msig("C2", "CP:KEGG_MEDICUS"))
cat("\n유전자세트 수: "); print(vapply(SETS, length, integer(1)))

run_gsea <- function(stat, sets, tag) {
  fg <- fgsea(sets, stat, minSize = 15, maxSize = 500, eps = 0)
  fg <- fg[order(fg$padj, -abs(fg$NES)), ]
  out <- as.data.frame(fg)
  out$leadingEdge <- vapply(out$leadingEdge, paste, character(1), collapse = ",")
  write.table(out, file.path(RES, sprintf("GSEA_%s.tsv", tag)), sep = "\t",
              quote = FALSE, row.names = FALSE)
  fg
}

report <- function(fg, label, n = 12) {
  sig <- fg[!is.na(fg$padj) & fg$padj < 0.05, ]
  cat(sprintf("\n=== %s : 유의 %d / %d ===\n", label, nrow(sig), nrow(fg)))
  if (!nrow(sig)) return(invisible())
  up <- sig[sig$NES > 0, ][order(-sig[sig$NES > 0]$NES), ]
  dn <- sig[sig$NES < 0, ][order( sig[sig$NES < 0]$NES), ]
  cat("-- 활성화 (NES>0) --\n")
  print(head(as.data.frame(up[, c("pathway","NES","padj","size")]), n), digits = 3, row.names = FALSE)
  cat("-- 억제 (NES<0) --\n")
  print(head(as.data.frame(dn[, c("pathway","NES","padj","size")]), n), digits = 3, row.names = FALSE)
}

cat("\n\n################ (A) 환자 : HGPS vs WT ################\n")
GP <- list()
for (nm in names(SETS)) {
  GP[[nm]] <- run_gsea(stat_pat, SETS[[nm]], paste0("patient_", nm))
  report(GP[[nm]], paste("환자", nm))
}

cat("\n\n################ (B) isogenic : progerin 6일 유도 ################\n")
GI <- list()
for (nm in names(SETS)) {
  GI[[nm]] <- run_gsea(stat_iso, SETS[[nm]], paste0("isogenic_", nm))
  report(GI[[nm]], paste("isogenic", nm), n = 10)
}

# ------------------------------------------------ (C) 급성 vs 만성 경로 비교 --
cat("\n\n################ (C) 급성(isogenic) vs 만성(환자) 경로 비교 ################\n")
cmp <- merge(as.data.frame(GP$Hallmark)[, c("pathway","NES","padj")],
             as.data.frame(GI$Hallmark)[, c("pathway","NES","padj")],
             by = "pathway", suffixes = c("_chronic", "_acute"))
cmp$both_sig <- cmp$padj_chronic < .05 & cmp$padj_acute < .05
cmp$agree    <- sign(cmp$NES_chronic) == sign(cmp$NES_acute)
cat(sprintf("Hallmark 50개 중 둘 다 유의: %d, 그중 방향 일치: %d\n",
            sum(cmp$both_sig), sum(cmp$both_sig & cmp$agree)))
cat(sprintf("NES Spearman 상관 (전체 50개): %.3f\n",
            cor(cmp$NES_chronic, cmp$NES_acute, method = "spearman")))
cat("\n-- 만성에서만 유의한 경로 (급성에는 없음) --\n")
co <- cmp[cmp$padj_chronic < .05 & cmp$padj_acute >= .05, ]
print(head(co[order(-abs(co$NES_chronic)), c("pathway","NES_chronic","NES_acute")], 12),
      digits = 3, row.names = FALSE)
cat("\n-- 급성에서만 유의한 경로 --\n")
ao <- cmp[cmp$padj_acute < .05 & cmp$padj_chronic >= .05, ]
print(head(ao[order(-abs(ao$NES_acute)), c("pathway","NES_chronic","NES_acute")], 12),
      digits = 3, row.names = FALSE)
cat("\n-- 방향이 반대인 경로 (주의) --\n")
op <- cmp[cmp$both_sig & !cmp$agree, ]
if (nrow(op)) print(op[, c("pathway","NES_chronic","NES_acute")], digits = 3, row.names = FALSE) else cat("  없음\n")
write.table(cmp, file.path(RES, "hallmark_acute_vs_chronic.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

png(file.path(FIG, "hallmark_acute_vs_chronic.png"), width = 7, height = 6, units = "in", res = 300, type = "cairo")
print(ggplot(cmp, aes(NES_chronic, NES_acute, color = both_sig)) +
  geom_point(size = 2, alpha = .8) +
  geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 0, linetype = 2) +
  geom_abline(slope = 1, intercept = 0, color = "grey60") +
  ggrepel::geom_text_repel(data = cmp[cmp$both_sig, ],
      aes(label = sub("^HALLMARK_", "", pathway)), size = 2.4, max.overlaps = 25, show.legend = FALSE) +
  scale_color_manual(values = c(`TRUE` = "#C44E52", `FALSE` = "grey75"),
                     labels = c("not both", "both significant"), name = NULL) +
  labs(title = "Hallmark: chronic (patient) vs acute (6-day progerin induction)",
       x = "NES  chronic (HGPS patients)", y = "NES  acute (TetON progerin)") +
  theme_bw(base_size = 11))
dev.off()

# ------------------------------------------------------------ ORA (DEG) -----
cat("\n\n################ (D) ORA : 환자 DEG 과대표현 ################\n")
de <- read.delim(file.path(RES, "DE_HGPS_vs_WT_patient.tsv"))
up  <- de$gene[de$sig == "UP"]; dn <- de$gene[de$sig == "DOWN"]
univ <- de$gene[!is.na(de$padj)]
e2 <- function(g) na.omit(mapIds(org.Hs.eg.db, keys = g, column = "ENTREZID",
                                 keytype = "SYMBOL", multiVals = "first"))
for (lab in c("UP", "DOWN")) {
  g <- if (lab == "UP") up else dn
  cat(sprintf("\n--- %s 유전자 %d개 ---\n", lab, length(g)))
  go <- enrichGO(e2(g), OrgDb = org.Hs.eg.db, universe = e2(univ), ont = "BP",
                 pAdjustMethod = "BH", qvalueCutoff = .1, readable = TRUE)
  if (!is.null(go) && nrow(as.data.frame(go))) {
    print(head(as.data.frame(go)[, c("Description","GeneRatio","p.adjust","Count")], 12),
          digits = 3, row.names = FALSE)
    write.table(as.data.frame(go), file.path(RES, sprintf("ORA_GO_%s.tsv", lab)),
                sep = "\t", quote = FALSE, row.names = FALSE)
  } else cat("  유의 결과 없음\n")
  rp <- tryCatch(enrichPathway(e2(g), universe = e2(univ), readable = TRUE), error = function(e) NULL)
  if (!is.null(rp) && nrow(as.data.frame(rp))) {
    cat("  [Reactome]\n")
    print(head(as.data.frame(rp)[, c("Description","p.adjust","Count")], 8), digits = 3, row.names = FALSE)
    write.table(as.data.frame(rp), file.path(RES, sprintf("ORA_Reactome_%s.tsv", lab)),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

# ---------------------------------------------------------- Hallmark 그림 ---
sh <- as.data.frame(GP$Hallmark); sh <- sh[!is.na(sh$padj) & sh$padj < .05, ]
if (nrow(sh)) {
  sh <- sh[order(sh$NES), ]
  sh$pathway <- factor(sub("^HALLMARK_", "", sh$pathway), levels = sub("^HALLMARK_", "", sh$pathway))
  png(file.path(FIG, "GSEA_hallmark_patient.png"), width = 8, height = max(4, .3 * nrow(sh)), units = "in", res = 300, type = "cairo")
  print(ggplot(sh, aes(NES, pathway, fill = NES > 0)) + geom_col(width = .7) +
    scale_fill_manual(values = c(`TRUE`="#C44E52", `FALSE`="#4C72B0"),
                      labels = c("Down in HGPS","Up in HGPS"), name = NULL) +
    labs(title = "Hallmark GSEA: HGPS vs control fibroblasts",
         subtitle = "GSE306264, 20 patient vs 10 control lines, FDR < 0.05",
         x = "Normalized Enrichment Score", y = NULL) + theme_bw(base_size = 10))
  dev.off()
}
cat(sprintf("\n저장 완료 -> %s\n", RES))
