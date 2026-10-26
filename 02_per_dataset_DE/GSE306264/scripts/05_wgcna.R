#!/usr/bin/env Rscript
# =============================================================================
# GSE306264 : WGCNA 공발현 모듈 -> 모듈-표현형(유전형/서명점수/하위군) 상관
#   30개 세포주(HGPS 20 / WT 10) 기준. 샘플수가 많지 않으므로 모듈 최소 크기를
#   크게 잡고(minModuleSize=50) 과분할을 피한다.
# =============================================================================
suppressPackageStartupMessages({
  library(WGCNA); library(DESeq2); library(matrixStats); library(pheatmap)
  library(RColorBrewer); library(msigdbr); library(clusterProfiler); library(org.Hs.eg.db)
})
options(stringsAsFactors = FALSE); enableWGCNAThreads(nThreads = 8)

# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
.this <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
DEPOSIT <- Sys.getenv("PROGERIA_DEPOSIT", normalizePath(file.path(dirname(.this), "..", "..", "..")))
BASE <- file.path(DEPOSIT, "02_per_dataset_DE", "GSE306264")  # [deposit] was ~/car_project/progeria_meta/GSE306264
RES  <- file.path(BASE, "results"); FIG <- file.path(RES, "figures")
set.seed(1)

o <- readRDS(file.path(RES, "GSE306264_objects.rds"))
vm <- assay(o$vsd); geno <- colData(o$dds)$genotype

# 변동 상위 8000 유전자
sel <- order(rowVars(vm), decreasing = TRUE)[1:min(8000, nrow(vm))]
datExpr <- t(vm[sel, ])
cat(sprintf("WGCNA 입력: 샘플 %d x 유전자 %d\n", nrow(datExpr), ncol(datExpr)))
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

# soft threshold
sft <- pickSoftThreshold(datExpr, powerVector = c(1:10, seq(12, 20, 2)),
                         networkType = "signed", verbose = 0)
pw <- sft$powerEstimate; if (is.na(pw)) pw <- 12
cat(sprintf("선택된 soft power = %d (R^2 = %.3f)\n", pw,
            sft$fitIndices$SFT.R.sq[sft$fitIndices$Power == pw]))

# 다른 패키지(DESeq2/S4Vectors 등)가 cor 를 마스킹하면 blockwiseModules 가 실패한다.
# WGCNA::cor 로 임시 치환 후 복원한다.
cor_backup <- cor; cor <- WGCNA::cor
net <- blockwiseModules(datExpr, power = pw, networkType = "signed",
                        TOMType = "signed", minModuleSize = 50,
                        mergeCutHeight = 0.25, numericLabels = TRUE,
                        maxBlockSize = 9000, verbose = 0, randomSeed = 1)
cor <- cor_backup
mods <- labels2colors(net$colors)
cat(sprintf("\n모듈 %d개 (grey 제외 %d개)\n", length(unique(mods)),
            length(setdiff(unique(mods), "grey"))))
print(sort(table(mods), decreasing = TRUE))

MEs <- orderMEs(moduleEigengenes(datExpr, mods)$eigengenes)

# ------------------------------------------------- 모듈-표현형 상관 ---------
score <- read.delim(file.path(RES, "signature_scores_per_line.tsv"), check.names = FALSE)
score <- as.data.frame(t(score))[rownames(datExpr), , drop = FALSE]
hs    <- read.delim(file.path(RES, "per_line_HGPS_score.tsv"))
traits <- data.frame(
  HGPS         = as.numeric(geno == "HGPS"),
  HGPS_score   = hs$score[match(rownames(datExpr), hs$line)],
  row.names    = rownames(datExpr)
)
traits <- cbind(traits, score[, c("SMC_program","Senescence_CDKN","SASP_core",
                                  "HALLMARK_MTORC1_SIGNALING","Ribosome"), drop = FALSE])

mt  <- cor(MEs, traits, use = "p")
mtp <- corPvalueStudent(mt, nrow(datExpr))
cat("\n=== 모듈-표현형 상관 (r / p) ===\n")
for (i in seq_len(nrow(mt))) {
  ln <- sprintf("%-16s", rownames(mt)[i])
  for (j in seq_len(ncol(mt))) ln <- paste0(ln, sprintf(" %s=%.2f(%.3f)", colnames(mt)[j], mt[i,j], mtp[i,j]))
  cat(ln, "\n")
}
write.table(data.frame(module = rownames(mt), mt, check.names = FALSE),
            file.path(RES, "wgcna_module_trait_cor.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

png(file.path(FIG, "wgcna_module_trait.png"), width = 7, height = max(4, .35 * nrow(mt)), units = "in", res = 300, type = "cairo")
pheatmap(mt, display_numbers = matrix(sprintf("%.2f\n(%.3f)", mt, mtp), nrow = nrow(mt)),
         cluster_cols = FALSE, fontsize_number = 6,
         col = colorRampPalette(rev(brewer.pal(11, "RdBu")))(255),
         main = "WGCNA module-trait correlation")
dev.off()

# ---------------------------------------- HGPS 와 가장 강하게 연관된 모듈 ----
ord <- order(-abs(mt[, "HGPS"]))
cat("\n=== HGPS 유전형과 상관이 큰 모듈 상위 5 ===\n")
for (i in head(ord, 5)) {
  m <- sub("^ME", "", rownames(mt)[i])
  g <- colnames(datExpr)[mods == m]
  cat(sprintf("\n[%s]  r = %.3f (p = %.3g), 유전자 %d개\n", m, mt[i,"HGPS"], mtp[i,"HGPS"], length(g)))
  # hub 유전자 (모듈 멤버십 상위)
  kME <- cor(datExpr[, g, drop = FALSE], MEs[, rownames(mt)[i]], use = "p")[,1]
  cat("  hub 유전자:", paste(head(names(sort(kME, decreasing = TRUE)), 12), collapse = ", "), "\n")
  eg <- na.omit(mapIds(org.Hs.eg.db, keys = g, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first"))
  uv <- na.omit(mapIds(org.Hs.eg.db, keys = colnames(datExpr), column = "ENTREZID", keytype = "SYMBOL", multiVals = "first"))
  go <- tryCatch(enrichGO(eg, OrgDb = org.Hs.eg.db, universe = uv, ont = "BP",
                          qvalueCutoff = .1, readable = TRUE), error = function(e) NULL)
  if (!is.null(go) && nrow(as.data.frame(go))) {
    cat("  GO:BP 상위:\n")
    print(head(as.data.frame(go)[, c("Description","p.adjust","Count")], 6), digits = 3, row.names = FALSE)
  } else cat("  GO 유의 결과 없음\n")
}

gm <- data.frame(gene = colnames(datExpr), module = mods)
write.table(gm, file.path(RES, "wgcna_gene_modules.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(net = net, MEs = MEs, mods = mods, mt = mt), file.path(RES, "wgcna_objects.rds"))
cat(sprintf("\n저장 완료 -> %s\n", RES))
