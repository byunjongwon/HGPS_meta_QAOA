# ─────────────────────────────────────────────────────────────
# 28. Sensitivity analysis — repeated measurements from the same donor in GSE231761
#
# Why
#   GSE231761 has four donor lines (HGPS 0031C/AG1B, control 0901C/BJ1C), each profiled
#   at passages 7 and 14. The primary model (02_per_dataset_DE/GSE231761, 01_deseq2)
#   collapses technical duplicates within donor x passage (16 -> 8 libraries) and fits
#   ~ passage + condition. Donor is nested within condition, so it cannot be a fixed
#   effect next to condition, and the two passages of a donor enter as if independent.
#
# What
#   The GSE231761 HGPS effect is re-estimated with two donor-aware models and the k=3
#   random-effects meta-analysis (05_meta_k3.R) is repeated on the same 12,290 genes:
#     S1 donor-level aggregation : counts summed over passages within donor -> 2 vs 2,
#                                  DESeq2 ~ condition (unshrunken MLE log2FC + lfcSE)
#     S2 limma-voom, donor block : ~ passage + condition, donor as blocking factor
#                                  (duplicateCorrelation consensus), moderated SE
#     S2b limma-voom, no block   : same as S2 without the block (separates the effect of
#                                  the block from the change of method DESeq2 -> limma)
#   GSE306264 and GSE198983 inputs are unchanged. Gene sets are defined as in 05_meta_k3.R;
#   Hallmark and GO:BP GSEA on z = meta_lfc / meta_se repeat 22_meta_GSEA.R exactly.
#
# Output (03_meta_analysis/results/)
#   sensitivity_GSE231761_summary.tsv     one row per model
#   sensitivity_GSE231761_hallmark.tsv    Hallmark NES / padj per model
#   sensitivity_GSE231761_genes.tsv       per-gene meta results per model
# ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
    library(DESeq2); library(metafor); library(limma); library(edgeR)
    library(data.table); library(msigdbr); library(fgsea)
    library(clusterProfiler); library(org.Hs.eg.db)
})

# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
.this <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
DEPOSIT <- Sys.getenv("PROGERIA_DEPOSIT", normalizePath(file.path(dirname(.this), "..", "..")))
ROOT <- file.path(DEPOSIT, "02_per_dataset_DE")
RES  <- file.path(DEPOSIT, "03_meta_analysis", "results")
set.seed(1)

# ------------------------------------------------ inputs (as in 05_meta_k3.R) ----
o1 <- readRDS(file.path(ROOT, "GSE306264/results/GSE306264_objects.rds"))
d1 <- as.data.frame(results(o1$dds, contrast = c("genotype","HGPS","WT"), alpha = .05))
d1$gene <- rownames(d1)
d1 <- d1[!is.na(d1$padj) & is.finite(d1$lfcSE) & d1$lfcSE > 0, ]
d1 <- d1[!duplicated(d1$gene), ]

o2 <- readRDS(file.path(ROOT, "GSE231761/results/GSE231761_objects.rds"))
d2 <- as.data.frame(results(o2$dds, contrast = c("condition","HGPS","Control"), alpha = .05))
d2$ensembl <- rownames(d2)
d2 <- merge(d2, unique(o2$res[, c("ensembl","gene")]), by = "ensembl")
d2 <- d2[!is.na(d2$padj) & is.finite(d2$lfcSE) & d2$lfcSE > 0 & !is.na(d2$gene) & d2$gene != "", ]
d2 <- d2[order(d2$padj), ]; d2 <- d2[!duplicated(d2$gene), ]
cnv <- read.delim(file.path(ROOT, "GSE231761/results/CNV_polymorphism_candidates.tsv"), check.names = FALSE)
d2 <- d2[!d2$gene %in% cnv$gene[cnv$splits_by_condition %in% c(TRUE, "TRUE")], ]

o3 <- readRDS(file.path(DEPOSIT, "01_rnaseq_processing/GSE198983_human_AoSMC/results/09_deseq2/dds.rds"))
d3 <- as.data.frame(results(o3, contrast = c("group","Progerin","Control"), alpha = .05))
d3$gene_id <- rownames(d3)
ann3 <- read.delim(file.path(DEPOSIT, "01_rnaseq_processing/GSE198983_human_AoSMC/results/09_deseq2/DE_Progerin_vs_Control_blocked.tsv"))
d3 <- merge(d3, ann3[, c("gene_id","symbol")], by = "gene_id")
names(d3)[names(d3) == "symbol"] <- "gene"
d3 <- d3[!is.na(d3$padj) & is.finite(d3$lfcSE) & d3$lfcSE > 0 & !is.na(d3$gene) & d3$gene != "", ]
d3 <- d3[order(d3$padj), ]; d3 <- d3[!duplicated(d3$gene), ]

common <- Reduce(intersect, list(d1$gene, d2$gene, d3$gene))
a  <- d1[match(common, d1$gene), ]
c3 <- d3[match(common, d3$gene), ]
ens <- d2$ensembl[match(common, d2$gene)]
cat(sprintf("genes evaluated in all three human datasets: %d\n", length(common)))

# ------------------------------------------- GSE231761 alternative estimates ----
dds <- o2$dds
cd  <- as.data.frame(colData(dds))
cat("\nGSE231761 libraries after technical-duplicate collapse:\n")
print(cd[, c("line","passage","condition")])

est <- list()
est$primary <- data.frame(lfc = d2$log2FoldChange[match(common, d2$gene)],
                          se  = d2$lfcSE[match(common, d2$gene)])

# S1: donor-level aggregation (each donor contributes one sample)
ddl <- collapseReplicates(dds, groupby = factor(cd$line))
colData(ddl)$condition <- droplevels(factor(colData(ddl)$condition, levels = c("Control","HGPS")))
design(ddl) <- ~ condition
ddl <- DESeq(ddl, quiet = TRUE)
r1  <- results(ddl, contrast = c("condition","HGPS","Control"))
est$donor_aggregated_DESeq2 <- data.frame(lfc = r1[ens, "log2FoldChange"], se = r1[ens, "lfcSE"])

# S2: limma-voom with donor as blocking factor
y  <- calcNormFactors(DGEList(counts(dds)))
X  <- model.matrix(~ passage + condition, cd)
v  <- voom(y, X)
dc <- duplicateCorrelation(v, X, block = cd$line)
v  <- voom(y, X, block = cd$line, correlation = dc$consensus)
dc <- duplicateCorrelation(v, X, block = cd$line)
rho_donor <- dc$consensus
cat(sprintf("\nwithin-donor correlation (duplicateCorrelation consensus) = %.3f\n", rho_donor))
fit <- eBayes(lmFit(v, X, block = cd$line, correlation = rho_donor))
cf  <- "conditionHGPS"
est$limma_donor_block <- data.frame(lfc = fit$coefficients[ens, cf],
                                    se  = (sqrt(fit$s2.post) * fit$stdev.unscaled[, cf])[ens])

# S2b: identical limma model without the donor block
fitn <- eBayes(lmFit(voom(y, X), X))
est$limma_no_block <- data.frame(lfc = fitn$coefficients[ens, cf],
                                 se  = (sqrt(fitn$s2.post) * fitn$stdev.unscaled[, cf])[ens])

# ------------------------------------------------------- k=3 meta-analysis ----
meta <- function(s) {
    mt <- t(vapply(seq_along(common), function(i) {
        yv <- c(a$log2FoldChange[i], s$lfc[i], c3$log2FoldChange[i])
        vv <- c(a$lfcSE[i], s$se[i], c3$lfcSE[i])^2
        if (any(!is.finite(c(yv, vv)))) return(rep(NA_real_, 4))
        re <- tryCatch(rma.uni(yi = yv, vi = vv, method = "DL"), error = function(e) NULL)
        if (is.null(re)) return(rep(NA_real_, 4))
        c(re$beta[1], re$se[1], re$pval, re$I2)
    }, numeric(4)))
    M <- data.frame(gene = common, meta_lfc = mt[, 1], meta_se = mt[, 2], meta_p = mt[, 3], I2 = mt[, 4])
    M$meta_padj <- p.adjust(M$meta_p, "BH")
    M$all_same_dir <- abs(rowSums(sign(cbind(a$log2FoldChange, s$lfc, c3$log2FoldChange)))) == 3
    M$supported  <- !is.na(M$meta_padj) & M$meta_padj < .05
    M$concordant <- M$supported & M$all_same_dir & M$I2 < 75
    M$effect05   <- M$concordant & abs(M$meta_lfc) >= .5
    M$effect1    <- M$concordant & abs(M$meta_lfc) >= 1
    M
}
R <- lapply(est, meta)
P <- R$primary

# the primary model must reproduce 05_meta_k3.R exactly
ref <- fread(file.path(RES, "meta_k3_final.tsv"))
ref <- ref[match(common, ref$gene)]
stopifnot(isTRUE(all.equal(ref$meta_lfc, P$meta_lfc, tolerance = 1e-8)),
          sum(ref$meta_padj < .05, na.rm = TRUE) == sum(P$supported))
cat(sprintf("primary model reproduces 05_meta_k3.R: %d supported / %d concordant / %d |LFC|>=0.5 / %d |LFC|>=1\n",
            sum(P$supported), sum(P$concordant), sum(P$effect05), sum(P$effect1)))

# ------------------------------------------------- Hallmark GSEA (as in 22) ----
hs <- msigdbr(species = "Homo sapiens", category = "H")
pw <- split(hs$gene_symbol, hs$gs_name)
# Same call sequence as 22_meta_GSEA.R (seed, then gseGO, then fgsea), so that the
# primary row reproduces GSEA_meta_Hallmark.tsv / GSEA_meta_GOBP.tsv exactly.
GO <- list()
gsea <- function(M) {
    M <- M[is.finite(M$meta_lfc) & is.finite(M$meta_se) & M$meta_se > 0, ]
    rk <- setNames(M$meta_lfc / M$meta_se, M$gene)
    rk <- sort(rk[!duplicated(names(rk))], decreasing = TRUE)
    set.seed(11)
    g <- gseGO(geneList = rk, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
               minGSSize = 15, maxGSSize = 500, pAdjustMethod = "BH",
               pvalueCutoff = 0.05, eps = 0, seed = TRUE)
    go <- as.data.frame(clusterProfiler::simplify(g, cutoff = 0.7, by = "p.adjust"))
    fg <- as.data.table(fgsea(pathways = pw, stats = rk, minSize = 10, maxSize = 500, nPermSimple = 10000))
    list(hallmark = fg, n_go = nrow(go))
}
GS <- lapply(R, gsea)
G  <- lapply(GS, `[[`, "hallmark")

off <- fread(file.path(RES, "GSEA_meta_Hallmark.tsv"))
chk <- merge(off[, .(pathway, padj)], G$primary[, .(pathway, padj)], by = "pathway")
stopifnot(isTRUE(all.equal(chk$padj.x, chk$padj.y, tolerance = 1e-8)))
cat(sprintf("primary Hallmark GSEA reproduces 22_meta_GSEA.R: %d sets FDR<0.05\n", sum(G$primary$padj < .05)))

# ----------------------------------------------------------------- outputs ----
sets <- c(supported = "supported", concordant = "concordant", effect05 = "effect05", effect1 = "effect1")
seven <- c("OGN","GCNT4","THRB","CAMK2B","G0S2","PPL","PTGES")
summ <- rbindlist(lapply(names(R), function(k) {
    M <- R[[k]]
    row <- list(model = k)
    for (s in names(sets)) {
        row[[paste0("n_", s)]] <- sum(M[[s]])
        row[[paste0("overlap_", s)]] <- sum(M[[s]] & P[[s]])
    }
    row$spearman_meta_lfc        <- cor(M$meta_lfc, P$meta_lfc, method = "spearman", use = "complete")
    row$spearman_GSE231761_lfc   <- cor(est[[k]]$lfc, est$primary$lfc, method = "spearman", use = "complete")
    row$median_GSE231761_se_ratio <- median(est[[k]]$se / est$primary$se, na.rm = TRUE)
    row$median_I2_all_genes      <- median(M$I2, na.rm = TRUE)
    row$cross_species7_supported <- sum(M$supported[match(seven, M$gene)])
    row$cross_species7_effect1   <- sum(M$effect1[match(seven, M$gene)])
    row$n_hallmark_fdr05         <- sum(G[[k]]$padj < .05)
    row$n_GOBP_nonredundant      <- GS[[k]]$n_go
    row$within_donor_correlation <- if (k == "limma_donor_block") rho_donor else NA_real_
    row
}))
fwrite(summ, file.path(RES, "sensitivity_GSE231761_summary.tsv"), sep = "\t")

hm <- rbindlist(lapply(names(G), function(k) G[[k]][, .(model = k, pathway, NES, pval, padj, size)]))
fwrite(hm[order(model, padj)], file.path(RES, "sensitivity_GSE231761_hallmark.tsv"), sep = "\t")

gt <- data.table(gene = common)
for (k in names(R)) {
    M <- R[[k]]
    gt[, (paste0(k, c("_GSE231761_lfc","_GSE231761_se","_meta_lfc","_meta_padj","_I2","_concordant"))) :=
           list(est[[k]]$lfc, est[[k]]$se, M$meta_lfc, M$meta_padj, M$I2, M$concordant)]
}
gt <- gt[R$primary$supported | Reduce(`|`, lapply(R, `[[`, "supported"))]
fwrite(gt, file.path(RES, "sensitivity_GSE231761_genes.tsv"), sep = "\t")

cat("\n================ summary ================\n")
print(summ[, .(model, n_supported, overlap_supported, n_concordant, overlap_concordant,
               n_effect1, overlap_effect1, rho = round(spearman_meta_lfc, 3),
               se_ratio = round(median_GSE231761_se_ratio, 2), cross_species7_supported,
               n_hallmark_fdr05, n_GOBP_nonredundant)])
key <- c("HALLMARK_MTORC1_SIGNALING","HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
         "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY","HALLMARK_PI3K_AKT_MTOR_SIGNALING",
         "HALLMARK_UV_RESPONSE_UP","HALLMARK_INTERFERON_GAMMA_RESPONSE")
cat("\nrecurrent pathways (NES / padj):\n")
print(dcast(hm[pathway %in% key, .(model, pathway = sub("HALLMARK_", "", pathway),
                                    v = sprintf("%.2f / %.1e", NES, padj))], pathway ~ model, value.var = "v"))
for (k in setdiff(names(R), "primary"))
    cat(sprintf("\n[%s] larger-effect (|LFC|>=1) genes not retained: %s", k,
                paste(setdiff(P$gene[P$effect1], R[[k]]$gene[R[[k]]$effect1]), collapse = ", ")))
cat("\n\ndone ->", RES, "\n")
