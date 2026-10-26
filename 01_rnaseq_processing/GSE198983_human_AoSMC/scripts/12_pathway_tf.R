#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────
# STEP 12 — 경로/전사인자 활성 추론
#   [PROGENy]  14개 신호전달 경로 활성 (반응성 유전자 footprint 기반)
#              → "경로 구성원 발현"이 아니라 "경로가 실제 켜졌는지"를 봅니다
#   [DoRothEA + decoupleR] 전사인자 활성 (regulon A~C 신뢰도)
#              → SASP를 이끄는 TF(NF-κB/RELA, STAT, E2F 등) 식별
#
# 산출: 샘플별 활성 히트맵 + 그룹 간 차등 활성 barplot
# ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
    library(data.table); library(ggplot2); library(pheatmap); library(RColorBrewer)
    library(org.Hs.eg.db); library(AnnotationDbi); library(decoupleR)
})
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])),
                 "_plot_setup.R"))

# [deposit] PROJ default = dataset dir (parent of this scripts/ dir; was ~/progerin_rnaseq); overridable by PROJ
PROJ    <- Sys.getenv("PROJ", normalizePath(file.path(dirname(sub("--file=", "",
               grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "..")))
RESULTS <- file.path(PROJ, "results")
OUT     <- file.path(RESULTS, "12_pathway_tf")
FIG     <- file.path(OUT, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
theme_set(theme_bw(base_size = 12))
GRP_COL <- c(Control = "#0072B2", Progerin = "#D55E00")

# STEP 09가 실제로 사용한 샘플 구성(collapse 반영)을 따릅니다.
# 없으면 원본 samples.tsv로 폴백합니다.
.cd <- file.path(RESULTS, "09_deseq2/coldata_used.tsv")
meta <- if (file.exists(.cd)) fread(.cd) else fread(file.path(PROJ, "meta/samples.tsv"))
grp  <- factor(meta$group, levels = c("Control", "Progerin"))
BATCH <- if ("batch" %in% names(meta) && uniqueN(meta$batch) > 1)
    factor(meta$batch) else NULL

# ── 입력 행렬 (SYMBOL) ────────────────────────────────────────
vst <- fread(file.path(RESULTS, "09_deseq2/vst_counts.tsv"))
m   <- as.matrix(vst[, -1]); rownames(m) <- vst$gene_id; m <- m[, meta$sample]
sym <- mapIds(org.Hs.eg.db, keys = rownames(m), keytype = "ENSEMBL",
              column = "SYMBOL", multiVals = "first")
keep <- !is.na(sym); ms <- m[keep, ]; rs <- sym[keep]
ord <- order(apply(ms, 1, var), decreasing = TRUE)
ms <- ms[ord, ]; rs <- rs[ord]; ms <- ms[!duplicated(rs), ]; rownames(ms) <- rs[!duplicated(rs)]

# 샘플별 활성 추론은 배치를 제거한 발현으로 합니다. 배치가 분산의 75.6%를
# 차지해서 보정 없이는 히트맵이 그룹이 아니라 배치로 갈립니다.
# (차등 활성은 아래에서 DESeq2의 배치 보정된 통계량을 쓰므로 영향 없음)
if (!is.null(BATCH)) {
    suppressPackageStartupMessages(library(limma))
    ms <- removeBatchEffect(ms, batch = BATCH, design = model.matrix(~ grp))
    cat("배치 보정 적용 (removeBatchEffect)\n")
}

# DE 통계 (t-like 통계량으로 차등 활성 추론)
de <- fread(file.path(RESULTS, "09_deseq2/DE_Progerin_vs_Control.tsv"))
de <- de[!is.na(symbol) & !is.na(log2FoldChange) & !is.na(lfcSE) & lfcSE > 0]
de[, stat := log2FoldChange / lfcSE]
de <- de[order(-abs(stat))][!duplicated(symbol)]
demat <- matrix(de$stat, ncol = 1, dimnames = list(de$symbol, "Progerin_vs_Control"))

ann <- data.frame(group = grp, row.names = colnames(ms))
if (!is.null(BATCH)) ann$batch <- BATCH

plot_acts <- function(acts_mat, tag, title) {
    png(file.path(FIG, sprintf("%s_heatmap.png", tag)), 8, min(2 + 0.25 * nrow(acts_mat), 20), units = "in", res = 300, type = "cairo")
    pheatmap(acts_mat, annotation_col = ann, annotation_colors = list(group = GRP_COL),
             color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(255),
             scale = "row", fontsize_row = 8, main = title)
    dev.off()
}

# ── PROGENy ───────────────────────────────────────────────────
if (requireNamespace("progeny", quietly = TRUE)) {
    net <- tryCatch(decoupleR::get_progeny(organism = "human", top = 500),
                    error = function(e) NULL)
    if (is.null(net)) {
        pm <- progeny::getModel("Human", top = 500)
        net <- data.table(source = rep(colnames(pm), each = nrow(pm)),
                          target = rep(rownames(pm), times = ncol(pm)),
                          weight = as.vector(as.matrix(pm)))
        net <- net[weight != 0]
    }
    setDT(net)

    # 샘플별 활성
    acts <- decoupleR::run_mlm(mat = ms, net = net, .source = "source", .target = "target",
                               .mor = "weight", minsize = 5)
    setDT(acts)
    am <- dcast(acts[statistic == "mlm"], source ~ condition, value.var = "score")
    amm <- as.matrix(am[, -1]); rownames(amm) <- am$source
    fwrite(data.table(pathway = rownames(amm), as.data.frame(amm)),
           file.path(OUT, "progeny_sample_activity.tsv"), sep = "\t")
    plot_acts(amm, "progeny", "PROGENy 경로 활성 (샘플별, row z-score)")

    # 차등 활성
    dacts <- decoupleR::run_mlm(mat = demat, net = net, .source = "source", .target = "target",
                                .mor = "weight", minsize = 5)
    setDT(dacts); dacts <- dacts[statistic == "mlm"]
    dacts[, sig := ifelse(p_value < 0.05, "p<0.05", "ns")]
    fwrite(dacts, file.path(OUT, "progeny_diff_activity.tsv"), sep = "\t")
    png(file.path(FIG, "progeny_diff_activity.png"), 7, 6, units = "in", res = 300, type = "cairo")
    print(ggplot(dacts, aes(reorder(source, score), score, fill = score)) +
          geom_col() + coord_flip() +
          scale_fill_gradient2(low = "#0072B2", mid = "grey90", high = "#D55E00") +
          geom_text(aes(label = ifelse(p_value < 0.05, "*", "")), hjust = -0.3, size = 5) +
          labs(x = NULL, y = "활성 점수 (Progerin vs Control)",
               title = "PROGENy 차등 경로 활성", subtitle = "* p<0.05"))
    dev.off()
    cat("PROGENy 완료 — 유의 경로:", sum(dacts$p_value < 0.05), "개\n")
    print(dacts[order(p_value), .(source, score, p_value)][1:8])
} else cat("progeny 미설치 — 건너뜀\n")

# ── DoRothEA 전사인자 활성 ───────────────────────────────────
tfnet <- tryCatch(decoupleR::get_collectri(organism = "human", split_complexes = FALSE),
                  error = function(e) NULL)
src <- "CollecTRI"
if (is.null(tfnet) && requireNamespace("dorothea", quietly = TRUE)) {
    d <- as.data.table(dorothea::dorothea_hs)
    tfnet <- d[confidence %in% c("A", "B", "C"), .(source = tf, target, mor)]
    src <- "DoRothEA (A/B/C)"
}
if (!is.null(tfnet)) {
    setDT(tfnet)
    if (!"mor" %in% names(tfnet)) tfnet[, mor := 1]

    acts <- decoupleR::run_ulm(mat = ms, net = tfnet, .source = "source", .target = "target",
                               .mor = "mor", minsize = 5)
    setDT(acts); acts <- acts[statistic == "ulm"]
    am <- dcast(acts, source ~ condition, value.var = "score")
    amm <- as.matrix(am[, -1]); rownames(amm) <- am$source
    fwrite(data.table(tf = rownames(amm), as.data.frame(amm)),
           file.path(OUT, "tf_sample_activity.tsv"), sep = "\t")

    # 분산 상위 30개 TF만 히트맵
    v <- apply(amm, 1, var)
    top <- names(sort(v, decreasing = TRUE))[1:min(30, length(v))]
    plot_acts(amm[top, , drop = FALSE], "tf",
              sprintf("전사인자 활성 상위 30 (%s, row z-score)", src))

    dacts <- decoupleR::run_ulm(mat = demat, net = tfnet, .source = "source", .target = "target",
                                .mor = "mor", minsize = 5)
    setDT(dacts); dacts <- dacts[statistic == "ulm"]
    dacts[, padj := p.adjust(p_value, "BH")]
    setorder(dacts, p_value)
    fwrite(dacts, file.path(OUT, "tf_diff_activity.tsv"), sep = "\t")

    show <- head(dacts, 30)
    png(file.path(FIG, "tf_diff_activity.png"), 7, 9, units = "in", res = 300, type = "cairo")
    print(ggplot(show, aes(reorder(source, score), score, fill = score)) +
          geom_col() + coord_flip() +
          scale_fill_gradient2(low = "#0072B2", mid = "grey90", high = "#D55E00") +
          geom_text(aes(label = ifelse(padj < 0.05, "*", "")), hjust = -0.3, size = 5) +
          labs(x = NULL, y = "활성 점수 (Progerin vs Control)",
               title = sprintf("차등 전사인자 활성 상위 30 (%s)", src),
               subtitle = "* FDR<0.05"))
    dev.off()
    cat("\nTF 활성 완료 — FDR<0.05:", sum(dacts$padj < 0.05), "개\n")
    print(dacts[1:12, .(source, score, p_value, padj)])
} else cat("CollecTRI/DoRothEA 없음 — TF 분석 건너뜀\n")

cat("\n완료 →", FIG, "\n")
