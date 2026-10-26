#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────
# STEP 17 — WGCNA 공발현 네트워크 모듈 분석
#
# ⚠️ 해석상 주의: n=8은 WGCNA 권장 최소치(보통 ≥15~20)에 못 미칩니다.
#    모듈-형질 상관은 자유도가 낮아 과대추정되기 쉬우므로,
#    "탐색적 결과"로만 보고 DEG/GSEA 결과와 교차확인해야 합니다.
#    이 경고는 리포트에도 그대로 남깁니다.
#
# 산출:
#   soft threshold 선택 그림, 모듈 덴드로그램, 모듈-형질 상관 히트맵,
#   모듈별 eigengene 박스플롯, 모듈별 GO 농축, hub 유전자 목록
# ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
    library(WGCNA); library(data.table); library(ggplot2); library(pheatmap)
    library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi)
    library(RColorBrewer)
})
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])),
                 "_plot_setup.R"))
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = max(1, as.integer(Sys.getenv("THREADS", "8"))))

# [deposit] PROJ default = dataset dir (parent of this scripts/ dir; was ~/progerin_rnaseq); overridable by PROJ
PROJ    <- Sys.getenv("PROJ", normalizePath(file.path(dirname(sub("--file=", "",
               grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "..")))
RESULTS <- file.path(PROJ, "results")
OUT     <- file.path(RESULTS, "17_wgcna")
FIG     <- file.path(OUT, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
theme_set(theme_bw(base_size = 12))
GRP_COL <- c(Control = "#0072B2", Progerin = "#D55E00")

# STEP 09가 실제로 사용한 샘플 구성(collapse 반영)을 따릅니다.
# 없으면 원본 samples.tsv로 폴백합니다.
.cd <- file.path(RESULTS, "09_deseq2/coldata_used.tsv")
meta <- if (file.exists(.cd)) fread(.cd) else fread(file.path(PROJ, "meta/samples.tsv"))
vst  <- fread(file.path(RESULTS, "09_deseq2/vst_counts.tsv"))
m <- as.matrix(vst[, -1]); rownames(m) <- vst$gene_id; m <- m[, meta$sample]

# 배치가 분산을 지배하면 모듈이 그룹이 아니라 배치를 따라 잡힙니다.
BATCH <- if ("batch" %in% names(meta) && uniqueN(meta$batch) > 1)
    factor(meta$batch) else NULL
if (!is.null(BATCH)) {
    suppressPackageStartupMessages(library(limma))
    m <- removeBatchEffect(m, batch = BATCH,
                           design = model.matrix(~ factor(meta$group)))
    cat("배치 보정 적용 (removeBatchEffect)\n")
}

# 분산 상위 5000 유전자 (WGCNA 표준 관행)
v <- apply(m, 1, var)
sel <- names(sort(v, decreasing = TRUE))[1:min(5000, length(v))]
datExpr <- t(m[sel, ])
N <- nrow(datExpr)
cat(sprintf("WGCNA 입력: %d 샘플 x %d 유전자\n", N, ncol(datExpr)))

# 표본수 경고 — 그림/리포트에 그대로 실립니다.
NWARN <- if (N < 6) {
    sprintf("⚠️ n=%d — 상관 기반 네트워크로는 신뢰할 수 없는 수준입니다. 결과는 참고용에 불과합니다.", N)
} else if (N < 15) {
    sprintf("⚠️ n=%d — WGCNA 권장치(≥15~20) 미만. 탐색적으로만 해석하세요.", N)
} else ""
if (nzchar(NWARN)) cat(NWARN, "\n")

gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
    cat("불량 샘플/유전자 제거 후:", dim(datExpr), "\n")
}

sym <- mapIds(org.Hs.eg.db, keys = colnames(datExpr), keytype = "ENSEMBL",
              column = "SYMBOL", multiVals = "first")

# ── 1) soft-threshold power 선택 ──────────────────────────────
powers <- c(1:10, seq(12, 30, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 0)
pw <- sft$powerEstimate
if (is.na(pw)) pw <- 12       # n이 작아 추정 실패 시 signed 권장 기본값
cat("선택된 soft power:", pw, "\n")

png(file.path(FIG, "01_soft_threshold.png"), 10, 5, units = "in", res = 300, type = "cairo")
par(mfrow = c(1, 2))
plot(sft$fitIndices$Power, -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
     xlab = "Soft threshold (power)", ylab = "Scale-free topology fit (signed R^2)",
     type = "n", main = "Scale independence")
text(sft$fitIndices$Power, -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
     labels = sft$fitIndices$Power, col = "red")
abline(h = 0.8, col = "blue", lty = 2)
plot(sft$fitIndices$Power, sft$fitIndices$mean.k., xlab = "Soft threshold (power)",
     ylab = "Mean connectivity", type = "n", main = "Mean connectivity")
text(sft$fitIndices$Power, sft$fitIndices$mean.k., labels = sft$fitIndices$Power, col = "red")
dev.off()

# ── 2) 네트워크 구성 + 모듈 검출 ─────────────────────────────
# ⚠️ WGCNA는 내부에서 cor()를 추가 인자(weights.x/cosine 등)와 함께 부릅니다.
# data.table/기타 패키지가 stats::cor를 앞세우면 "사용되지 않은 인자" 에러가
# 납니다. 표준 우회책대로 WGCNA::cor를 전역에 올려두고 끝나면 되돌립니다.
cor <- WGCNA::cor

net <- blockwiseModules(datExpr, power = pw, networkType = "signed",
                        TOMType = "signed", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        maxBlockSize = 6000, verbose = 0, randomSeed = 42)
mod_col <- labels2colors(net$colors)
cat("모듈 수:", length(unique(mod_col)), "\n")
print(table(mod_col))

png(file.path(FIG, "02_dendrogram.png"), 10, 6, units = "in", res = 300, type = "cairo")
plotDendroAndColors(net$dendrograms[[1]], mod_col[net$blockGenes[[1]]],
                    "Module", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "유전자 클러스터 덴드로그램 + 모듈")
dev.off()

gene_mod <- data.table(gene_id = colnames(datExpr), symbol = sym[colnames(datExpr)],
                       module = mod_col)
fwrite(gene_mod, file.path(OUT, "gene_module_assignment.tsv"), sep = "\t")

# ── 3) 모듈-형질 상관 ─────────────────────────────────────────
MEs <- orderMEs(moduleEigengenes(datExpr, mod_col)$eigengenes)
traits <- data.frame(Progerin = as.numeric(meta$group == "Progerin"),
                     row.names = meta$sample)

modTraitCor <- cor(MEs, traits, use = "p")
modTraitP   <- corPvalueStudent(modTraitCor, nrow(datExpr))
res <- data.table(module = rownames(modTraitCor), cor = modTraitCor[, 1], p = modTraitP[, 1])
res[, padj := p.adjust(p, "BH")]
setorder(res, p)
fwrite(res, file.path(OUT, "module_trait_correlation.tsv"), sep = "\t")
print(res)

txt <- matrix(sprintf("%.2f\n(p=%.2g)", modTraitCor, modTraitP), nrow = nrow(modTraitCor))
png(file.path(FIG, "03_module_trait.png"), 5, max(4, 0.35 * nrow(modTraitCor)), units = "in", res = 300, type = "cairo")
par(mar = c(5, 9, 3, 2))
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traits), yLabels = rownames(modTraitCor),
               ySymbols = rownames(modTraitCor), colorLabels = FALSE,
               colors = blueWhiteRed(50), textMatrix = txt, setStdMargins = FALSE,
               cex.text = 0.7, zlim = c(-1, 1),
               main = sprintf("모듈–형질 상관 (n=%d)", N))
dev.off()

# ── 4) 모듈 eigengene 박스플롯 ───────────────────────────────
me_long <- melt(data.table(sample = rownames(MEs), MEs, group = meta$group),
                id.vars = c("sample", "group"), variable.name = "module", value.name = "ME")
png(file.path(FIG, "04_eigengene_boxplot.png"), 11, 8, units = "in", res = 300, type = "cairo")
print(ggplot(me_long, aes(group, ME, color = group)) +
      geom_boxplot(outlier.shape = NA) + geom_jitter(width = .15, size = 2) +
      facet_wrap(~module, scales = "free_y") +
      scale_color_manual(values = GRP_COL) +
      labs(x = NULL, y = "Module eigengene", title = "모듈별 eigengene"))
dev.off()

# ── 5) 유의 모듈의 GO 농축 + hub 유전자 ──────────────────────
kME <- as.data.table(signedKME(datExpr, MEs))
kME[, gene_id := colnames(datExpr)]
kME[, symbol := sym[colnames(datExpr)]]
fwrite(kME, file.path(OUT, "kME.tsv"), sep = "\t")

sig_mods <- res[p < 0.05, module]
cor <- stats::cor    # 여기서부터는 표준 cor 사용
cat("\n유의 모듈(p<0.05):", paste(sig_mods, collapse = ", "), "\n")

universe <- na.omit(unique(AnnotationDbi::mapIds(org.Hs.eg.db, keys = colnames(datExpr),
                                                 keytype = "ENSEMBL", column = "ENTREZID")))
hub_all <- list()
for (mm in sig_mods) {
    col <- sub("^ME", "", mm)
    ids <- gene_mod[module == col, gene_id]
    if (length(ids) < 15) next

    # hub 유전자 = 해당 모듈 kME 상위
    k <- kME[gene_id %in% ids]
    kcol <- paste0("kME", col)
    if (kcol %in% names(k)) {
        setorderv(k, kcol, -1)
        hub <- head(k[, .(gene_id, symbol, kME = get(kcol))], 25)
        hub[, module := col]
        hub_all[[col]] <- hub
    }

    ent <- na.omit(AnnotationDbi::mapIds(org.Hs.eg.db, keys = ids, keytype = "ENSEMBL",
                                         column = "ENTREZID"))
    e <- tryCatch(enrichGO(ent, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
                           universe = universe, pvalueCutoff = 0.05, readable = TRUE),
                  error = function(x) NULL)
    if (!is.null(e) && nrow(as.data.frame(e)) > 0) {
        fwrite(as.data.frame(e), file.path(OUT, sprintf("GO_module_%s.tsv", col)), sep = "\t")
        png(file.path(FIG, sprintf("05_GO_module_%s.png", col)), 9, 7, units = "in", res = 300, type = "cairo")
        print(dotplot(e, showCategory = 15) + ggtitle(sprintf("모듈 %s — GO:BP (%d 유전자)", col, length(ids))))
        dev.off()
    }
}
if (length(hub_all)) fwrite(rbindlist(hub_all), file.path(OUT, "hub_genes.tsv"), sep = "\t")

if (nzchar(NWARN)) cat("\n", NWARN, "\n")
cat("완료 →", FIG, "\n")
