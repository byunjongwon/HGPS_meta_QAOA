#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────
# STEP 12. 기능 분석 — GO over-representation + GSEA
#
# 산출 그림:
#   01_go_dotplot_{비교}.pdf     GO BP over-representation (up/down 분리)
#   02_go_barplot_{비교}.pdf     상위 GO term bar
#   03_gsea_dotplot_{비교}.pdf   GSEA (Hallmark + GO BP)
#   04_gsea_running_{비교}.pdf   상위 gene set의 running enrichment score
#   05_gsea_heatmap_NES.pdf      세 비교의 NES 비교 히트맵
#   06_go_emapplot_{비교}.pdf    term 간 유사도 네트워크
#   07_splicing_gene_GO.pdf      ★ rMATS 유의 event 유전자의 GO
#
# 산출 표: GO_{비교}.tsv, GSEA_{비교}.tsv
# ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
    library(clusterProfiler); library(org.Mm.eg.db); library(enrichplot)
    library(fgsea); library(msigdbr); library(data.table); library(ggplot2)
})

# UTF-8(Δ, 한글) 렌더링을 위한 cairo PDF 장치 설정
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)[1])), "_plot_setup.R"))

# [deposit] default PROJ = this dataset dir (scripts/..), still overridable by PROJ (was ~/lmna_rnaseq)
.this <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
PROJ <- Sys.getenv("PROJ", normalizePath(file.path(dirname(.this), "..")))
DE   <- file.path(PROJ, "results/08_deseq2")
OUT  <- file.path(PROJ, "results/12_enrichment")
FIG  <- file.path(OUT, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

COMPS <- c("LMNA_vs_WTyng", "WTold_vs_WTyng", "LMNA_vs_WTold")
FDR <- 0.05; LFC <- 1

# Hallmark gene set (마우스)
hm <- tryCatch(msigdbr(species = "Mus musculus", collection = "H"),
               error = function(e) msigdbr(species = "Mus musculus", category = "H"))
hm_list <- split(hm$ensembl_gene, hm$gs_name)

nes_all <- list()

for (cp in COMPS) {
    f <- file.path(DE, sprintf("DE_%s.tsv", cp))
    if (!file.exists(f)) { message("건너뜀: ", f); next }
    d <- fread(f)
    d <- d[!is.na(padj) & !is.na(log2FoldChange)]

    universe <- d$gene_id
    up   <- d[padj < FDR & log2FoldChange >  LFC, gene_id]
    down <- d[padj < FDR & log2FoldChange < -LFC, gene_id]
    cat(sprintf("%-16s up %d / down %d\n", cp, length(up), length(down)))

    # ── GO over-representation ────────────────────────────────
    go_res <- list()
    for (dir in c("up", "down")) {
        genes <- if (dir == "up") up else down
        if (length(genes) < 10) next
        e <- tryCatch(enrichGO(genes, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL",
                               ont = "BP", universe = universe,
                               pAdjustMethod = "BH", qvalueCutoff = 0.1,
                               readable = TRUE),
                      error = function(x) NULL)
        if (!is.null(e) && nrow(as.data.frame(e)) > 0) {
            go_res[[dir]] <- e
            dt <- as.data.table(as.data.frame(e)); dt[, direction := dir]
            fwrite(dt, file.path(OUT, sprintf("GO_%s_%s.tsv", cp, dir)), sep = "\t")
        }
    }
    if (length(go_res) > 0) {
        pdf(file.path(FIG, sprintf("01_go_dotplot_%s.pdf", cp)), 9, 4 + 3 * length(go_res))
        for (dir in names(go_res)) {
            print(dotplot(go_res[[dir]], showCategory = 15) +
                  ggtitle(sprintf("GO BP — %s (%s regulated)", cp, dir)))
        }
        dev.off()

        pdf(file.path(FIG, sprintf("02_go_barplot_%s.pdf", cp)), 9, 4 + 3 * length(go_res))
        for (dir in names(go_res)) {
            print(barplot(go_res[[dir]], showCategory = 15) +
                  ggtitle(sprintf("GO BP — %s (%s)", cp, dir)))
        }
        dev.off()

        for (dir in names(go_res)) {
            tryCatch({
                e2 <- pairwise_termsim(go_res[[dir]])
                pdf(file.path(FIG, sprintf("06_go_emapplot_%s_%s.pdf", cp, dir)), 9, 8)
                print(emapplot(e2, showCategory = 30) +
                      ggtitle(sprintf("GO term network — %s (%s)", cp, dir)))
                dev.off()
            }, error = function(x) { try(dev.off(), silent = TRUE) })
        }
    }

    # ── GSEA ──────────────────────────────────────────────────
    # lfcShrink(type="ashr")는 Wald 통계량(stat) 컬럼을 반환하지 않으므로
    # 방향(log2FC 부호) x 유의성(-log10 p)으로 랭킹 지표를 만듭니다.
    rk_dt <- d[!is.na(pvalue) & !is.na(log2FoldChange)]
    rk_dt[, rank_metric := sign(log2FoldChange) * -log10(pmax(pvalue, 1e-300))]
    rk <- setNames(rk_dt$rank_metric, rk_dt$gene_id)
    rk <- sort(rk[!duplicated(names(rk))], decreasing = TRUE)

    gse_go <- tryCatch(gseGO(rk, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL",
                             ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                             verbose = FALSE, eps = 0),
                       error = function(x) NULL)
    fg <- tryCatch(fgsea(hm_list, rk, minSize = 10, maxSize = 500, nPermSimple = 10000),
                   error = function(x) NULL)

    if (!is.null(fg) && nrow(fg) > 0) {
        fg <- as.data.table(fg)[order(padj)]
        fwrite(fg[, .(pathway, pval, padj, NES, size)],
               file.path(OUT, sprintf("GSEA_hallmark_%s.tsv", cp)), sep = "\t")
        sig <- fg[padj < 0.25][order(-abs(NES))][seq_len(min(20, .N))]
        if (nrow(sig) > 0) {
            sig[, pathway := gsub("^HALLMARK_", "", pathway)]
            pdf(file.path(FIG, sprintf("03_gsea_dotplot_%s.pdf", cp)), 9, 7)
            print(ggplot(sig, aes(NES, reorder(pathway, NES),
                                  size = size, color = padj)) +
                  geom_point() + geom_vline(xintercept = 0, lty = 2) +
                  scale_color_gradient(low = "#D55E00", high = "grey70") +
                  theme_bw(base_size = 11) +
                  labs(y = NULL, title = sprintf("GSEA Hallmark — %s", cp)))
            dev.off()
            nes_all[[cp]] <- fg[, .(pathway = gsub("^HALLMARK_", "", pathway), NES, padj)]
        }
    }

    if (!is.null(gse_go) && nrow(as.data.frame(gse_go)) > 0) {
        fwrite(as.data.table(as.data.frame(gse_go)),
               file.path(OUT, sprintf("GSEA_GO_%s.tsv", cp)), sep = "\t")
        pdf(file.path(FIG, sprintf("04_gsea_running_%s.pdf", cp)), 8, 6)
        for (i in seq_len(min(5, nrow(as.data.frame(gse_go))))) {
            tryCatch(print(gseaplot2(gse_go, geneSetID = i,
                                     title = as.data.frame(gse_go)$Description[i])),
                     error = function(x) NULL)
        }
        dev.off()
    }
}

# ── 05 NES 히트맵 ─────────────────────────────────────────────
if (length(nes_all) >= 2) {
    m <- Reduce(function(a, b) merge(a, b, by = "pathway", all = TRUE),
                lapply(names(nes_all), function(n) {
                    x <- nes_all[[n]][, .(pathway, NES)]; setnames(x, "NES", n); x
                }))
    mm <- as.matrix(m[, -1]); rownames(mm) <- m$pathway
    keep <- rowSums(!is.na(mm)) == ncol(mm) & apply(abs(mm), 1, max, na.rm = TRUE) > 1
    if (sum(keep) >= 3) {
        pdf(file.path(FIG, "05_gsea_heatmap_NES.pdf"), 7, 0.25 * sum(keep) + 3)
        pheatmap::pheatmap(mm[keep, , drop = FALSE], cluster_cols = FALSE,
                           color = colorRampPalette(c("#0072B2", "white", "#D55E00"))(100),
                           main = "GSEA Hallmark NES", fontsize_row = 8)
        dev.off()
    }
}

# ── 07 ★ splicing 유전자의 GO ─────────────────────────────────
sig_f <- file.path(PROJ, "results/09_rmats_figures/significant_events_LMNA_vs_WTyng.tsv")
if (file.exists(sig_f)) {
    sp <- fread(sig_f)
    sp_genes <- unique(sub("\\..*$", "", gsub('"', '', sp$GeneID)))
    bg <- unique(sub("\\..*$", "", fread(file.path(DE, "DE_LMNA_vs_WTyng.tsv"))$gene_id))
    e <- tryCatch(enrichGO(sp_genes, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL",
                           ont = "BP", universe = bg, qvalueCutoff = 0.1, readable = TRUE),
                  error = function(x) NULL)
    if (!is.null(e) && nrow(as.data.frame(e)) > 0) {
        fwrite(as.data.table(as.data.frame(e)),
               file.path(OUT, "GO_splicing_genes_LMNA_vs_WTyng.tsv"), sep = "\t")
        pdf(file.path(FIG, "07_splicing_gene_GO.pdf"), 9, 7)
        print(dotplot(e, showCategory = 20) +
              ggtitle("차등 스플라이싱 유전자의 GO BP (LMNA vs WTyng)"))
        dev.off()
        cat("★ splicing 유전자 GO 완료\n")
    }
}

cat("완료 → ", FIG, "\n")
