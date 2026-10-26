#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────
# STEP 10 — 기능 분석
#   [ORA]  GO(BP/MF/CC), KEGG, Reactome — up/down 각각 + 전체
#   [GSEA] Hallmark(H), C2:CP(KEGG/Reactome/WikiPathways), C5:BP, C6, C8
#          → 순위 기반이라 임계값에 의존하지 않음 (권장 해석 근거)
#
# 산출 그림: dotplot, barplot, cnetplot, emapplot, ridgeplot, GSEA running score
# 산출 표:   각 분석의 전체 결과 tsv
# ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
    library(clusterProfiler); library(org.Hs.eg.db); library(ReactomePA)
    library(enrichplot); library(msigdbr); library(data.table); library(ggplot2)
    library(DOSE); library(AnnotationDbi)
})
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])),
                 "_plot_setup.R"))

# [deposit] PROJ default = dataset dir (parent of this scripts/ dir; was ~/progerin_rnaseq); overridable by PROJ
PROJ    <- Sys.getenv("PROJ", normalizePath(file.path(dirname(sub("--file=", "",
               grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "..")))
RESULTS <- file.path(PROJ, "results")
OUT     <- file.path(RESULTS, "10_enrichment")
FIG     <- file.path(OUT, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
theme_set(theme_bw(base_size = 11))

FDR <- 0.05; LFC <- 1
de <- fread(file.path(RESULTS, "09_deseq2/DE_Progerin_vs_Control.tsv"))
de <- de[!is.na(padj)]

# ENSEMBL → ENTREZ
e2 <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = de$gene_id,
                                             keytype = "ENSEMBL", columns = "ENTREZID"))
setDT(e2); e2 <- e2[!is.na(ENTREZID)][!duplicated(ENSEMBL)]
de <- merge(de, e2, by.x = "gene_id", by.y = "ENSEMBL", all.x = TRUE)

universe <- unique(de[!is.na(ENTREZID), ENTREZID])
sig      <- de[padj < FDR & abs(log2FoldChange) > LFC & !is.na(ENTREZID)]
up       <- unique(sig[log2FoldChange > 0, ENTREZID])
down     <- unique(sig[log2FoldChange < 0, ENTREZID])
allsig   <- unique(sig$ENTREZID)
cat(sprintf("ORA 입력 — up %d / down %d / 전체 %d (universe %d)\n",
            length(up), length(down), length(allsig), length(universe)))

save_plot <- function(p, file, w = 9, h = 7) {
    png(file.path(FIG, file), w, h, units = "in", res = 300, type = "cairo"); print(p); dev.off()
}

# ── ORA ───────────────────────────────────────────────────────
run_ora <- function(genes, tag) {
    if (length(genes) < 10) { cat("건너뜀(유전자 부족):", tag, "\n"); return(invisible(NULL)) }

    for (ont in c("BP", "MF", "CC")) {
        e <- tryCatch(enrichGO(genes, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                               ont = ont, universe = universe,
                               pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                               readable = TRUE), error = function(e) NULL)
        if (is.null(e) || nrow(as.data.frame(e)) == 0) next
        fwrite(as.data.frame(e), file.path(OUT, sprintf("ORA_GO%s_%s.tsv", ont, tag)), sep = "\t")
        es <- tryCatch(simplify(e, cutoff = 0.7, by = "p.adjust"), error = function(x) e)
        save_plot(dotplot(es, showCategory = 20) + ggtitle(sprintf("GO:%s — %s", ont, tag)),
                  sprintf("ORA_GO%s_%s_dot.pdf", ont, tag))
        if (ont == "BP") {
            save_plot(barplot(es, showCategory = 20) + ggtitle(sprintf("GO:BP — %s", tag)),
                      sprintf("ORA_GOBP_%s_bar.pdf", tag))
            ex <- tryCatch(pairwise_termsim(es), error = function(x) NULL)
            if (!is.null(ex) && nrow(as.data.frame(ex)) > 2)
                save_plot(emapplot(ex, showCategory = 30) + ggtitle(sprintf("GO:BP network — %s", tag)),
                          sprintf("ORA_GOBP_%s_emap.pdf", tag), 10, 9)
            # enrichplot 1.30에서 cex_label_gene 인자가 제거되어 사용하지 않습니다.
            fc <- setNames(sig$log2FoldChange, sig$ENTREZID)
            save_plot(cnetplot(es, showCategory = 8, foldChange = fc) +
                          ggtitle(sprintf("GO:BP gene-concept — %s", tag)),
                      sprintf("ORA_GOBP_%s_cnet.pdf", tag), 11, 9)
        }
    }

    k <- tryCatch(enrichKEGG(genes, organism = "hsa", universe = universe,
                             pvalueCutoff = 0.05), error = function(e) NULL)
    if (!is.null(k) && nrow(as.data.frame(k)) > 0) {
        k <- setReadable(k, org.Hs.eg.db, "ENTREZID")
        fwrite(as.data.frame(k), file.path(OUT, sprintf("ORA_KEGG_%s.tsv", tag)), sep = "\t")
        save_plot(dotplot(k, showCategory = 20) + ggtitle(sprintf("KEGG — %s", tag)),
                  sprintf("ORA_KEGG_%s_dot.pdf", tag))
    }

    r <- tryCatch(enrichPathway(genes, organism = "human", universe = universe,
                                pvalueCutoff = 0.05, readable = TRUE), error = function(e) NULL)
    if (!is.null(r) && nrow(as.data.frame(r)) > 0) {
        fwrite(as.data.frame(r), file.path(OUT, sprintf("ORA_Reactome_%s.tsv", tag)), sep = "\t")
        save_plot(dotplot(r, showCategory = 20) + ggtitle(sprintf("Reactome — %s", tag)),
                  sprintf("ORA_Reactome_%s_dot.pdf", tag))
    }
}

run_ora(up,     "up")
run_ora(down,   "down")
run_ora(allsig, "all")

# ── GSEA (순위 기반, 임계값 비의존) ──────────────────────────
rk <- de[!is.na(ENTREZID) & !is.na(log2FoldChange)]
# 동일 ENTREZ 중복 시 |log2FC| 최대만 사용
rk <- rk[order(-abs(log2FoldChange))][!duplicated(ENTREZID)]
gl <- setNames(rk$log2FoldChange, rk$ENTREZID)
gl <- sort(gl, decreasing = TRUE)
cat(sprintf("GSEA 입력: %d 유전자\n", length(gl)))

set.seed(42)

gse_go <- tryCatch(gseGO(gl, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
                         pvalueCutoff = 0.05, eps = 0, seed = TRUE), error = function(e) NULL)
if (!is.null(gse_go) && nrow(as.data.frame(gse_go)) > 0) {
    gse_go <- setReadable(gse_go, org.Hs.eg.db, "ENTREZID")
    fwrite(as.data.frame(gse_go), file.path(OUT, "GSEA_GOBP.tsv"), sep = "\t")
    save_plot(dotplot(gse_go, showCategory = 20, split = ".sign") + facet_grid(. ~ .sign) +
                  ggtitle("GSEA GO:BP"), "GSEA_GOBP_dot.pdf", 11, 8)
    save_plot(ridgeplot(gse_go, showCategory = 20) + ggtitle("GSEA GO:BP — 분포"),
              "GSEA_GOBP_ridge.pdf", 10, 9)
}

gse_kegg <- tryCatch(gseKEGG(gl, organism = "hsa", pvalueCutoff = 0.05, eps = 0, seed = TRUE),
                     error = function(e) NULL)
if (!is.null(gse_kegg) && nrow(as.data.frame(gse_kegg)) > 0) {
    gse_kegg <- setReadable(gse_kegg, org.Hs.eg.db, "ENTREZID")
    fwrite(as.data.frame(gse_kegg), file.path(OUT, "GSEA_KEGG.tsv"), sep = "\t")
    save_plot(dotplot(gse_kegg, showCategory = 20, split = ".sign") + facet_grid(. ~ .sign) +
                  ggtitle("GSEA KEGG"), "GSEA_KEGG_dot.pdf", 11, 8)
}

gse_re <- tryCatch(gsePathway(gl, organism = "human", pvalueCutoff = 0.05, eps = 0, seed = TRUE),
                   error = function(e) NULL)
if (!is.null(gse_re) && nrow(as.data.frame(gse_re)) > 0) {
    gse_re <- setReadable(gse_re, org.Hs.eg.db, "ENTREZID")
    fwrite(as.data.frame(gse_re), file.path(OUT, "GSEA_Reactome.tsv"), sep = "\t")
    save_plot(dotplot(gse_re, showCategory = 20, split = ".sign") + facet_grid(. ~ .sign) +
                  ggtitle("GSEA Reactome"), "GSEA_Reactome_dot.pdf", 11, 8)
}

# ── MSigDB 컬렉션 GSEA ───────────────────────────────────────
# msigdbr 버전에 따라 인자명이 다릅니다 (category/collection).
msig <- function(...) {
    a <- list(...)
    r <- tryCatch(do.call(msigdbr, c(list(species = "Homo sapiens"), a)), error = function(e) NULL)
    if (is.null(r)) {
        names(a)[names(a) == "category"]     <- "collection"
        names(a)[names(a) == "subcategory"]  <- "subcollection"
        r <- tryCatch(do.call(msigdbr, c(list(species = "Homo sapiens"), a)), error = function(e) NULL)
    }
    r
}

COLLECTIONS <- list(
    Hallmark      = list(category = "H"),
    C2_CP_KEGG    = list(category = "C2", subcategory = "CP:KEGG_LEGACY"),
    C2_CP_REACTOME= list(category = "C2", subcategory = "CP:REACTOME"),
    C2_CP_WP      = list(category = "C2", subcategory = "CP:WIKIPATHWAYS"),
    C5_GOBP       = list(category = "C5", subcategory = "GO:BP"),
    C6_oncogenic  = list(category = "C6")
)

for (cn in names(COLLECTIONS)) {
    m <- do.call(msig, COLLECTIONS[[cn]])
    if (is.null(m) || !nrow(m)) { cat("컬렉션 없음:", cn, "\n"); next }
    idcol <- intersect(c("ncbi_gene", "entrez_gene"), names(m))[1]
    if (is.na(idcol)) { cat("ENTREZ 컬럼 없음:", cn, "\n"); next }
    t2g <- unique(as.data.frame(m)[, c("gs_name", idcol)])
    names(t2g) <- c("term", "gene")
    t2g$gene <- as.character(t2g$gene)

    g <- tryCatch(GSEA(gl, TERM2GENE = t2g, pvalueCutoff = 0.05, eps = 0, seed = TRUE),
                  error = function(e) NULL)
    if (is.null(g) || nrow(as.data.frame(g)) == 0) { cat("유의 결과 없음:", cn, "\n"); next }
    fwrite(as.data.frame(g), file.path(OUT, sprintf("GSEA_%s.tsv", cn)), sep = "\t")
    cat(sprintf("%-16s 유의 gene set %d개\n", cn, nrow(as.data.frame(g))))

    save_plot(dotplot(g, showCategory = 20, split = ".sign") + facet_grid(. ~ .sign) +
                  ggtitle(sprintf("GSEA — %s", cn)), sprintf("GSEA_%s_dot.pdf", cn), 12, 8)

    # 상위 gene set의 running enrichment score
    d <- as.data.frame(g)
    topn <- head(d[order(d$p.adjust), "ID"], 6)
    if (length(topn)) {
        png(file.path(FIG, sprintf("GSEA_%s_runningscore.png", cn)), 8, 6, units = "in", res = 300, type = "cairo")
        for (id in topn) print(gseaplot2(g, geneSetID = id, title = id))
        dev.off()
    }
}

cat("\n완료 →", FIG, "\n")
