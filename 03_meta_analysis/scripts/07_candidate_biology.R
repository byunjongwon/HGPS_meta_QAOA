#!/usr/bin/env Rscript
# =============================================================================
# 최종 후보 유전자의 생물학적 해석 — 논문 Results 작성용
#   (A) 후보 집합별 GO/Reactome 과대표현
#   (B) 16개 핵심 후보의 기능 주석 + 방향/급만성 라벨
#   (C) 데이터셋 공통 경로 (3개 사람 데이터셋 GSEA 교집합)
#   (D) 논문에 넣을 핵심 수치 표
# =============================================================================
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi)
  library(ReactomePA); library(ggplot2)
})
# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
.this <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
DEPOSIT <- Sys.getenv("PROGERIA_DEPOSIT", normalizePath(file.path(dirname(.this), "..", "..")))
ROOT <- file.path(DEPOSIT, "02_per_dataset_DE")  # [deposit] GSE* dataset dirs (was ~/car_project/progeria_meta)
RES  <- file.path(DEPOSIT, "03_meta_analysis/results"); FIG <- file.path(RES, "figures")  # [deposit]
PNG <- function(f, w, h) png(f, width=w, height=h, units="in", res=300, type="cairo")

M    <- read.delim(file.path(RES, "meta_k3_final.tsv"))
cand <- read.delim(file.path(RES, "candidates_k3_final.tsv"))
base <- M$meta_padj < .05 & M$all_same_dir & M$I2 < 75
setB <- M[base & abs(M$meta_lfc) >= 1, ]
setC <- M[base & abs(M$meta_lfc) >= .5, ]
univ <- M$gene

e2 <- function(g) na.omit(mapIds(org.Hs.eg.db, keys=unique(g), column="ENTREZID",
                                 keytype="SYMBOL", multiVals="first"))

cat("############ (A) 후보 집합별 과대표현 ############\n")
for (nm in c("core16","panel34","panel97")) {
  g <- switch(nm, core16=cand$gene, panel34=setB$gene, panel97=setC$gene)
  cat(sprintf("\n=== %s (n=%d) ===\n", nm, length(g)))
  go <- tryCatch(enrichGO(e2(g), OrgDb=org.Hs.eg.db, universe=e2(univ), ont="BP",
                          pvalueCutoff=.05, qvalueCutoff=.2, readable=TRUE),
                 error=function(e) NULL)
  if (!is.null(go) && nrow(as.data.frame(go))) {
    d <- as.data.frame(go)[, c("Description","GeneRatio","p.adjust","Count","geneID")]
    print(head(d[, 1:4], 12), digits=3, row.names=FALSE)
    write.table(d, file.path(RES, sprintf("ORA_GO_%s.tsv", nm)), sep="\t",
                quote=FALSE, row.names=FALSE)
  } else cat("  유의 결과 없음\n")
  rp <- tryCatch(enrichPathway(e2(g), universe=e2(univ), readable=TRUE),
                 error=function(e) NULL)
  if (!is.null(rp) && nrow(as.data.frame(rp))) {
    cat("  [Reactome]\n")
    print(head(as.data.frame(rp)[, c("Description","p.adjust","Count")], 8),
          digits=3, row.names=FALSE)
  }
}

cat("\n\n############ (B) 핵심 16개 기능 주석 ############\n")
ann <- AnnotationDbi::select(org.Hs.eg.db, keys=cand$gene, keytype="SYMBOL",
                             columns=c("GENENAME","ENTREZID"))
ann <- ann[!duplicated(ann$SYMBOL), ]
t16 <- merge(cand[, c("gene","meta_lfc","meta_padj","I2","lfc_306264","lfc_231761",
                      "lfc_198983","lfc_mouse","mouse_sig","chronicity")],
             ann, by.x="gene", by.y="SYMBOL", all.x=TRUE)
t16 <- t16[order(-abs(t16$meta_lfc)), ]
t16$direction <- ifelse(t16$meta_lfc > 0, "UP", "DOWN")
print(t16[, c("gene","direction","meta_lfc","meta_padj","I2","mouse_sig","chronicity","GENENAME")],
      digits=3, row.names=FALSE)
write.table(t16, file.path(RES, "TABLE1_core16_annotated.tsv"), sep="\t",
            quote=FALSE, row.names=FALSE)

cat("\n\n############ (C) 3개 사람 데이터셋 공통 경로 ############\n")
gs <- list(
  GSE306264 = file.path(ROOT, "GSE306264/results/GSEA_patient_Hallmark.tsv"),
  GSE231761 = file.path(ROOT, "GSE231761/results/GSEA_condition_Hallmark.tsv"))
H <- lapply(gs, function(p) if (file.exists(p)) read.delim(p) else NULL)
H <- H[!vapply(H, is.null, logical(1))]
# GSE198983 Hallmark GSEA (자체 파이프라인 산출물)
p3 <- Sys.glob(file.path(DEPOSIT, "01_rnaseq_processing/GSE198983_human_AoSMC/results/10_enrichment/*Hallmark*.tsv"))  # [deposit]
p3 <- p3[!grepl("runningscore", p3)]
if (length(p3)) {
  h3 <- read.delim(p3[1])
  nmcol <- intersect(c("ID","pathway","Description"), names(h3))[1]
  h3$pathway <- h3[[nmcol]]
  H[["GSE198983"]] <- h3
}
cat("사용 가능한 Hallmark GSEA:", paste(names(H), collapse=", "), "\n")

if (length(H) >= 2) {
  nes <- Reduce(function(a, b) merge(a, b, by="pathway"),
                lapply(names(H), function(k) {
                  d <- H[[k]]
                  cc <- intersect(c("NES","nes"), names(d))[1]
                  pc <- intersect(c("padj","p.adjust","qvalue"), names(d))[1]
                  setNames(data.frame(d$pathway, d[[cc]], d[[pc]]),
                           c("pathway", paste0("NES_", k), paste0("padj_", k)))
                }))
  ncol_nes <- grep("^NES_", names(nes), value=TRUE)
  pcol <- grep("^padj_", names(nes), value=TRUE)
  nes$n_sig <- rowSums(nes[, pcol, drop=FALSE] < .05, na.rm=TRUE)
  nes$all_same_dir <- apply(sign(nes[, ncol_nes, drop=FALSE]), 1,
                            function(x) length(unique(x[!is.na(x)])) == 1)
  nes$mean_NES <- rowMeans(nes[, ncol_nes, drop=FALSE], na.rm=TRUE)
  com <- nes[nes$n_sig == length(ncol_nes) & nes$all_same_dir, ]
  com <- com[order(-abs(com$mean_NES)), ]
  cat(sprintf("\n%d개 데이터셋 전부에서 유의 + 방향 일치한 Hallmark 경로: %d개\n",
              length(ncol_nes), nrow(com)))
  print(com[, c("pathway", ncol_nes, "mean_NES")], digits=3, row.names=FALSE)
  write.table(nes, file.path(RES, "TABLE2_hallmark_across_datasets.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)

  if (nrow(com) > 0) {
    d <- com[, c("pathway", ncol_nes)]
    long <- reshape2::melt(d, id.vars="pathway", variable.name="dataset", value.name="NES")
    long$dataset <- sub("^NES_", "", long$dataset)
    long$pathway <- factor(sub("^HALLMARK_", "", long$pathway),
                           levels=rev(sub("^HALLMARK_", "", com$pathway)))
    PNG(file.path(FIG, "common_hallmark_pathways.png"), 8.5, max(3.5, .34*nrow(com)))
    print(ggplot(long, aes(NES, pathway, fill=dataset)) +
      geom_col(position=position_dodge(.78), width=.72) +
      geom_vline(xintercept=0, color="grey35") +
      scale_fill_manual(values=c("#C44E52","#DD8452","#55A868"), name=NULL) +
      labs(title="Hallmark pathways significant and concordant in all human datasets",
           subtitle="GSEA FDR < 0.05 in every dataset, same direction",
           x="Normalized Enrichment Score", y=NULL) +
      theme_bw(base_size=10) + theme(legend.position="bottom"))
    dev.off()
    cat("  -> common_hallmark_pathways.png 저장\n")
  }
}

cat("\n\n############ (D) 논문 핵심 수치 ############\n")
nums <- data.frame(
  item = c("사람 3개 데이터셋 공통 유전자",
           "메타 FDR<0.05", "3개 방향 일치", "I2<75% 통과",
           "|meta log2FC|>=0.5 (panel97)", "|meta log2FC|>=1 (panel34)",
           "마우스 방향 일치 (core16)", "마우스 유의+방향 일치 (strict7)",
           "이질성 I2 중앙값(%)",
           "쌍별 상관 306264-231761", "쌍별 상관 306264-198983",
           "쌍별 상관 231761-198983"),
  value = c(nrow(M),
            sum(M$meta_padj < .05, na.rm=TRUE),
            sum(M$meta_padj < .05 & M$all_same_dir, na.rm=TRUE),
            sum(base, na.rm=TRUE),
            nrow(setC), nrow(setB), nrow(cand),
            sum(base & abs(M$meta_lfc) >= 1 & M$mouse_same_dir, na.rm=TRUE),
            round(median(M$I2, na.rm=TRUE), 1),
            round(cor(M$lfc_306264[M$meta_padj<.05], M$lfc_231761[M$meta_padj<.05],
                      method="spearman"), 3),
            round(cor(M$lfc_306264[M$meta_padj<.05], M$lfc_198983[M$meta_padj<.05],
                      method="spearman"), 3),
            round(cor(M$lfc_231761[M$meta_padj<.05], M$lfc_198983[M$meta_padj<.05],
                      method="spearman"), 3)))
print(nums, row.names=FALSE)
write.table(nums, file.path(RES, "paper_key_numbers.tsv"), sep="\t",
            quote=FALSE, row.names=FALSE)
cat(sprintf("\n저장 -> %s\n", RES))
