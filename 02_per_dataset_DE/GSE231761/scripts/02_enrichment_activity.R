#!/usr/bin/env Rscript
# =============================================================================
# GSE231761 종합 분석 (1/2) — 경로 / 활성 / 서명
#   (A) 성별 확인
#   (B) GSEA : Hallmark / GO:BP / Reactome / KEGG  — 질환축 + passage축
#   (C) ORA
#   (D) 샘플별 서명 점수
#   (E) PROGENy 경로활성 / DoRothEA 전사인자 활성
# =============================================================================
suppressPackageStartupMessages({
  library(DESeq2); library(fgsea); library(msigdbr); library(ggplot2); library(ggrepel)
  library(clusterProfiler); library(ReactomePA); library(org.Hs.eg.db); library(AnnotationDbi)
  library(pheatmap); library(RColorBrewer); library(matrixStats)
})

# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
.this <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
DEPOSIT <- Sys.getenv("PROGERIA_DEPOSIT", normalizePath(file.path(dirname(.this), "..", "..", "..")))
BASE <- file.path(DEPOSIT, "02_per_dataset_DE", "GSE231761")  # [deposit] was ~/car_project/progeria_meta/GSE231761
RES  <- file.path(BASE, "results"); FIG <- file.path(RES, "figures")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
set.seed(1)

o   <- readRDS(file.path(RES, "GSE231761_objects.rds"))
dds <- o$dds; vsd <- o$vsd
map <- unique(o$res[!is.na(o$res$gene) & o$res$gene != "", c("ensembl", "gene")])
sym <- setNames(map$gene, map$ensembl)

cd  <- as.data.frame(colData(dds))
cond <- cd$condition; pas <- cd$passage
cat(sprintf("샘플 %d (HGPS %d / Control %d), 유전자 %d\n",
            ncol(dds), sum(cond=="HGPS"), sum(cond=="Control"), nrow(dds)))
print(cd[, c("line","passage","condition")])

# symbol 기준 VST 행렬
vm <- assay(vsd); rownames(vm) <- sym[rownames(vm)]
vm <- vm[!is.na(rownames(vm)), ]
vm <- vm[!duplicated(rownames(vm)), ]
cat(sprintf("symbol 매핑 VST 행렬: %d 유전자\n", nrow(vm)))
ann <- data.frame(condition = cond, passage = pas, row.names = colnames(vm))

# =============================================================================
# (A) 성별
# =============================================================================
cat("\n\n############ (A) 성별 ############\n")
sg <- intersect(c("XIST","RPS4Y1","DDX3Y","UTY","KDM5D","USP9Y","EIF1AY"), rownames(vm))
sx <- as.data.frame(t(vm[sg, , drop = FALSE]))
yc <- setdiff(sg, "XIST"); sx$Y_mean <- rowMeans(sx[, yc, drop = FALSE])
sx$sex <- ifelse(sx$Y_mean > sx$XIST, "male", "female")
sx$line <- cd$line; sx$condition <- cond
print(round(sx[, c(sg, "Y_mean")], 2))
cat("\n=== 세포주별 추정 성별 ===\n")
print(unique(sx[, c("line","condition","sex")]), row.names = FALSE)
tb <- table(sx$sex, sx$condition); print(tb)
if (nrow(tb) > 1 && all(dim(tb) == c(2,2)))
  cat(sprintf("Fisher p = %.3f %s\n", fisher.test(tb)$p.value,
              ifelse(fisher.test(tb)$p.value < .05, "<- 성별이 조건과 교락!", "<- 균형")))
write.table(sx, file.path(RES, "sex_assignment.tsv"), sep = "\t", quote = FALSE)

# =============================================================================
# (B) GSEA
# =============================================================================
cat("\n\n############ (B) GSEA ############\n")
msig <- function(coll, sub = NULL) {
  m <- if (is.null(sub)) msigdbr(species="Homo sapiens", collection=coll)
       else msigdbr(species="Homo sapiens", collection=coll, subcollection=sub)
  s <- split(m$gene_symbol, m$gs_name); s[lengths(s) >= 15 & lengths(s) <= 500]
}
SETS <- list(Hallmark = msig("H"), GO_BP = msig("C5","GO:BP"),
             Reactome = msig("C2","CP:REACTOME"), KEGG = msig("C2","CP:KEGG_MEDICUS"))

rank_of <- function(contrast) {
  r <- as.data.frame(results(dds, contrast = contrast, alpha = .05))
  r$gene <- sym[rownames(r)]
  r <- r[!is.na(r$stat) & !is.na(r$gene) & r$gene != "", ]
  r <- r[order(-abs(r$stat)), ]; r <- r[!duplicated(r$gene), ]
  sort(setNames(r$stat, r$gene), decreasing = TRUE)
}
st_cond <- rank_of(c("condition","HGPS","Control"))
st_pas  <- rank_of(c("passage","P14","P7"))
cat(sprintf("질환축 랭킹 %d / passage축 랭킹 %d 유전자\n", length(st_cond), length(st_pas)))

do_gsea <- function(stat, tag) {
  res <- list()
  for (nm in names(SETS)) {
    fg <- fgsea(SETS[[nm]], stat, minSize=15, maxSize=500, eps=0)
    fg <- fg[order(fg$padj, -abs(fg$NES)), ]
    out <- as.data.frame(fg); out$leadingEdge <- vapply(out$leadingEdge, paste, character(1), collapse=",")
    write.table(out, file.path(RES, sprintf("GSEA_%s_%s.tsv", tag, nm)), sep="\t", quote=FALSE, row.names=FALSE)
    sig <- fg[!is.na(fg$padj) & fg$padj < .05, ]
    cat(sprintf("\n=== [%s] %s : 유의 %d / %d ===\n", tag, nm, nrow(sig), nrow(fg)))
    if (nrow(sig)) {
      up <- sig[sig$NES > 0, ][order(-sig[sig$NES > 0]$NES), ]
      dn <- sig[sig$NES < 0, ][order( sig[sig$NES < 0]$NES), ]
      if (nrow(up)) { cat("-- 활성화 --\n"); print(head(as.data.frame(up[,c("pathway","NES","padj","size")]),10), digits=3, row.names=FALSE) }
      if (nrow(dn)) { cat("-- 억제 --\n");   print(head(as.data.frame(dn[,c("pathway","NES","padj","size")]),10), digits=3, row.names=FALSE) }
    }
    res[[nm]] <- fg
  }
  res
}
cat("\n\n---------------- 질환축 : HGPS vs Control ----------------")
GC <- do_gsea(st_cond, "condition")
cat("\n\n---------------- passage축 : P14 vs P7 ----------------")
GPa <- do_gsea(st_pas, "passage")

# 질환축 vs passage축 Hallmark 비교
cmp <- merge(as.data.frame(GC$Hallmark)[,c("pathway","NES","padj")],
             as.data.frame(GPa$Hallmark)[,c("pathway","NES","padj")],
             by="pathway", suffixes=c("_disease","_passage"))
cat(sprintf("\n\n=== Hallmark: 질환축 vs passage축 NES Spearman = %.3f ===\n",
            cor(cmp$NES_disease, cmp$NES_passage, method="spearman")))
cat("  (양수면 조로증 변화가 계대노화 변화와 같은 방향 = 가속노화 가설과 부합)\n")
write.table(cmp, file.path(RES, "hallmark_disease_vs_passage.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

png(file.path(FIG, "hallmark_disease_vs_passage.png"), width=7, height=6, units = "in", res = 300, type = "cairo")
cmp$lab <- ifelse(cmp$padj_disease < .05 | cmp$padj_passage < .05, sub("^HALLMARK_","",cmp$pathway), NA)
print(ggplot(cmp, aes(NES_disease, NES_passage)) +
  geom_point(aes(color = padj_disease < .05), size=2) +
  geom_hline(yintercept=0, linetype=2) + geom_vline(xintercept=0, linetype=2) +
  geom_abline(slope=1, intercept=0, color="grey60") +
  geom_text_repel(aes(label=lab), size=2.4, max.overlaps=25, show.legend=FALSE) +
  scale_color_manual(values=c(`TRUE`="#C44E52",`FALSE`="grey70"), name=NULL,
                     labels=c("ns","disease FDR<0.05")) +
  labs(title="Hallmark: HGPS disease effect vs replicative passage effect",
       subtitle=sprintf("GSE231761 VSMC;  Spearman = %.3f",
                        cor(cmp$NES_disease, cmp$NES_passage, method="spearman")),
       x="NES  HGPS vs Control", y="NES  P14 vs P7") + theme_bw(base_size=11))
dev.off()

# =============================================================================
# (C) ORA
# =============================================================================
cat("\n\n############ (C) ORA ############\n")
de <- read.delim(file.path(RES, "DE_HGPS_vs_Control_VSMC.tsv"))
de <- de[!is.na(de$gene) & de$gene != "", ]
univ <- de$gene[!is.na(de$padj)]
e2 <- function(g) na.omit(mapIds(org.Hs.eg.db, keys=unique(g), column="ENTREZID", keytype="SYMBOL", multiVals="first"))
for (lab in c("UP","DOWN")) {
  g <- de$gene[de$sig == lab]
  cat(sprintf("\n--- %s %d개 ---\n", lab, length(g)))
  go <- tryCatch(enrichGO(e2(g), OrgDb=org.Hs.eg.db, universe=e2(univ), ont="BP",
                          qvalueCutoff=.1, readable=TRUE), error=function(e) NULL)
  if (!is.null(go) && nrow(as.data.frame(go))) {
    print(head(as.data.frame(go)[,c("Description","GeneRatio","p.adjust","Count")],12), digits=3, row.names=FALSE)
    write.table(as.data.frame(go), file.path(RES, sprintf("ORA_GO_%s.tsv", lab)), sep="\t", quote=FALSE, row.names=FALSE)
  } else cat("  유의 결과 없음\n")
}

# =============================================================================
# (D) 서명 점수
# =============================================================================
cat("\n\n############ (D) 서명 점수 ############\n")
H <- msigdbr(species="Homo sapiens", collection="H"); hs <- split(H$gene_symbol, H$gs_name)
custom <- list(
  SASP_core       = c("IL6","IL1A","IL1B","CXCL8","CCL2","MMP1","MMP3","MMP10","IGFBP3","IGFBP7","SERPINE1","CXCL1","CXCL2"),
  Senescence_CDKN = c("CDKN1A","CDKN2A","CDKN2B","TP53","GLB1","TP53BP1"),
  Lamin_NE        = c("LMNA","LMNB1","LMNB2","LBR","EMD","TMPO","SUN1","SUN2","SYNE1","SYNE2"),
  VSMC_contractile= c("MYH11","ACTA2","CNN1","TAGLN","MYOCD","LMOD1","ACTG2"),
  VSMC_synthetic  = c("KLF4","SPP1","VCAN","COL1A1","FN1","MMP2","ELN"),
  Ribosome        = grep("^RP[LS][0-9]", rownames(vm), value=TRUE),
  Proliferation   = c("MKI67","PCNA","TOP2A","CCNB1","CDK1","MCM2","E2F1")
)
sets <- c(hs[c("HALLMARK_MTORC1_SIGNALING","HALLMARK_GLYCOLYSIS","HALLMARK_MYOGENESIS",
               "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_INTERFERON_ALPHA_RESPONSE",
               "HALLMARK_TNFA_SIGNALING_VIA_NFKB","HALLMARK_P53_PATHWAY","HALLMARK_G2M_CHECKPOINT",
               "HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_HYPOXIA")], custom)
z <- t(scale(t(vm)))
score <- t(vapply(sets, function(g){ g <- intersect(g, rownames(z))
  if (length(g) < 5) return(rep(NA_real_, ncol(z))); colMeans(z[g,,drop=FALSE]) }, numeric(ncol(z))))
colnames(score) <- colnames(vm)
st <- data.frame(set=rownames(score),
                 Ctrl_P7 = rowMeans(score[, cond=="Control" & pas=="P7", drop=FALSE]),
                 Ctrl_P14= rowMeans(score[, cond=="Control" & pas=="P14",drop=FALSE]),
                 HGPS_P7 = rowMeans(score[, cond=="HGPS"    & pas=="P7", drop=FALSE]),
                 HGPS_P14= rowMeans(score[, cond=="HGPS"    & pas=="P14",drop=FALSE]))
st$disease_effect <- rowMeans(score[, cond=="HGPS",drop=FALSE]) - rowMeans(score[, cond=="Control",drop=FALSE])
st$passage_effect <- rowMeans(score[, pas=="P14", drop=FALSE]) - rowMeans(score[, pas=="P7", drop=FALSE])
st <- st[order(-abs(st$disease_effect)), ]
print(st, digits=3, row.names=FALSE)
write.table(st, file.path(RES,"signature_scores_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(score, file.path(RES,"signature_scores_per_sample.tsv"), sep="\t", quote=FALSE)

png(file.path(FIG,"signature_scores_heatmap.png"), width=7.5, height=5.5, units = "in", res = 300, type = "cairo")
pheatmap(score[!is.na(rowSums(score)),], annotation_col=ann,
         col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(255),
         main="Signature scores (GSE231761 VSMC)", fontsize_row=8)
dev.off()

# =============================================================================
# (E) PROGENy / DoRothEA
# =============================================================================
cat("\n\n############ (E) 경로 / 전사인자 활성 ############\n")
pg <- tryCatch({ library(progeny); progeny(vm, scale=TRUE, organism="Human", top=500, perm=1) },
               error=function(e){ cat("PROGENy 실패:", conditionMessage(e), "\n"); NULL })
if (!is.null(pg)) {
  pt <- data.frame(pathway=colnames(pg),
                   disease = colMeans(pg[cond=="HGPS",,drop=FALSE]) - colMeans(pg[cond=="Control",,drop=FALSE]),
                   passage = colMeans(pg[pas=="P14",,drop=FALSE])   - colMeans(pg[pas=="P7",,drop=FALSE]))
  pt <- pt[order(-abs(pt$disease)), ]
  cat("\n=== PROGENy (질환 효과 / passage 효과) ===\n"); print(pt, digits=3, row.names=FALSE)
  write.table(pt, file.path(RES,"progeny_activity.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  png(file.path(FIG,"progeny_activity.png"), width=7, height=5, units = "in", res = 300, type = "cairo")
  pheatmap(t(pg), annotation_col=ann, main="PROGENy pathway activity",
           col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(255)); dev.off()
}

tf <- tryCatch({
  library(decoupleR); library(dorothea)
  net <- dorothea::dorothea_hs; net <- net[net$confidence %in% c("A","B","C"), ]
  names(net)[names(net)=="tf"] <- "source"
  decoupleR::run_ulm(mat=vm, net=net, .source="source", .target="target", .mor="mor", minsize=10)
}, error=function(e){ cat("DoRothEA 실패:", conditionMessage(e), "\n"); NULL })
if (!is.null(tf)) {
  act <- reshape2::acast(as.data.frame(tf), source ~ condition, value.var="score")
  act <- act[complete.cases(act), , drop=FALSE]
  cc <- cond[match(colnames(act), colnames(vm))]; pp <- pas[match(colnames(act), colnames(vm))]
  tt <- data.frame(TF=rownames(act),
                   disease = rowMeans(act[,cc=="HGPS",drop=FALSE]) - rowMeans(act[,cc=="Control",drop=FALSE]),
                   passage = rowMeans(act[,pp=="P14",drop=FALSE])  - rowMeans(act[,pp=="P7",drop=FALSE]))
  tt <- tt[order(-abs(tt$disease)), ]
  cat(sprintf("\n=== DoRothEA TF 활성 (%d개) — 질환 효과 상위 25 ===\n", nrow(tt)))
  print(head(tt, 25), digits=3, row.names=FALSE)
  write.table(tt, file.path(RES,"dorothea_TF_activity.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  png(file.path(FIG,"TF_activity_heatmap.png"), width=7.5, height=7, units = "in", res = 300, type = "cairo")
  pheatmap(act[head(tt$TF,30),], annotation_col=ann, main="DoRothEA TF activity (top 30 by disease effect)",
           col=colorRampPalette(rev(brewer.pal(11,"RdBu")))(255), fontsize_row=8); dev.off()
}
cat(sprintf("\n저장 완료 -> %s\n", RES))
