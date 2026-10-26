#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────
# 22. 메타분석 순위 기반 GSEA (GO:BP + Hallmark)
#
# 왜 ORA 가 아니라 GSEA 인가
#   21번에서 287개 코어를 배경 12,290개에 대해 ORA 했더니 유의 term 이 4개뿐이었다.
#   원인은 검정력이다. 임계값으로 잘라낸 287개만 쓰면 나머지 12,003개가 담고 있는
#   효과 순서 정보를 통째로 버리게 된다.
#   메타분석은 유전자마다 통합 효과와 그 표준오차를 이미 갖고 있으므로,
#   z = meta_lfc / meta_se 로 12,290개 전체를 부호 있는 순위로 세울 수 있다.
#   이 순위를 그대로 GSEA 에 넣는 것이 메타분석 결과의 표준 후속 분석이다.
#
# 산출
#   GSEA_meta_GOBP.tsv, GSEA_meta_Hallmark.tsv
# ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
    library(data.table); library(clusterProfiler); library(org.Hs.eg.db)
    library(msigdbr); library(fgsea)
})

# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
.this <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
DEPOSIT <- Sys.getenv("PROGERIA_DEPOSIT", normalizePath(file.path(dirname(.this), "..", "..")))
P <- file.path(DEPOSIT, "03_meta_analysis", "results")  # [deposit] was .../progeria_gse_원본/meta/results
M <- fread(file.path(P, "meta_k3_final.tsv"))
M <- M[is.finite(meta_lfc) & is.finite(meta_se) & meta_se > 0]

# 부호 있는 통합 통계량으로 순위를 만든다 (효과 크기 x 정밀도)
M[, z := meta_lfc / meta_se]
rk <- setNames(M$z, M$gene)
rk <- sort(rk[!duplicated(names(rk))], decreasing = TRUE)
cat(sprintf("순위 유전자 %d개  (z 범위 %.2f ~ %.2f)\n", length(rk), min(rk), max(rk)))

set.seed(11)
g <- gseGO(geneList = rk, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
           minGSSize = 15, maxGSSize = 500, pAdjustMethod = "BH",
           pvalueCutoff = 0.05, eps = 0, seed = TRUE)
gd <- as.data.frame(clusterProfiler::simplify(g, cutoff = 0.7, by = "p.adjust"))
fwrite(gd, file.path(P, "GSEA_meta_GOBP.tsv"), sep = "\t")
cat(sprintf("GO:BP  유의 %d개 (중복 제거 후)\n", nrow(gd)))

hs <- msigdbr(species = "Homo sapiens", category = "H")
pw <- split(hs$gene_symbol, hs$gs_name)
fg <- fgsea(pathways = pw, stats = rk, minSize = 10, maxSize = 500, nPermSimple = 10000)
fg <- as.data.table(fg)[order(padj)]
fg[, leadingEdge := sapply(leadingEdge, paste, collapse = ";")]
fwrite(fg, file.path(P, "GSEA_meta_Hallmark.tsv"), sep = "\t")
cat(sprintf("Hallmark 유의(padj<0.05) %d개 / 전체 %d개\n", sum(fg$padj < 0.05), nrow(fg)))
print(head(fg[padj < 0.05, .(pathway, NES, padj)], 12))
cat("완료 ->", P, "\n")
