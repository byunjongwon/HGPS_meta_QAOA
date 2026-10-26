#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────
# 21. 일치 코어(287개) 의 GO:BP 과대표현 분석
#
# 왜 새로 도는가
#   기존 ORA 는 core16(2 terms) / panel34(4 terms) / panel97(63 terms) 만 있었다.
#   앞의 두 개는 유전자 수가 적어 유의 term 이 거의 없고, 그림으로 쓸 수 없다.
#   논문에서 "공통 유전자의 기능" 을 보여줄 대상은
#   '3개 데이터 방향 일치 + I2<75%' 를 만족한 287개 코어다.
#
# 배경집합
#   전사체 전체가 아니라 '3개 데이터에서 모두 추정 가능했던 12,290개' 를 쓴다.
#   메타분석 대상이 될 자격 자체가 없던 유전자를 배경에 넣으면
#   검정 가능 여부가 만들어낸 편향이 기능 신호로 둔갑한다.
# ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
    library(data.table); library(clusterProfiler); library(org.Hs.eg.db)
})

# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
.this <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
DEPOSIT <- Sys.getenv("PROGERIA_DEPOSIT", normalizePath(file.path(dirname(.this), "..", "..")))
P <- file.path(DEPOSIT, "03_meta_analysis", "results")  # [deposit] was .../progeria_gse_원본/meta/results
M   <- fread(file.path(P, "meta_k3_final.tsv"))

universe <- M$gene                                   # 12,290 evaluable genes
core <- M[meta_padj < 0.05 & all_same_dir == TRUE & I2 < 75]
up   <- core[meta_lfc > 0]; dn <- core[meta_lfc < 0]

cat(sprintf("배경 %d개 / 코어 %d개 (up %d, down %d)\n",
            length(universe), nrow(core), nrow(up), nrow(dn)))

run <- function(g, tag) {
    e <- tryCatch(enrichGO(gene = g, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                           universe = universe, ont = "BP",
                           pAdjustMethod = "BH", pvalueCutoff = 0.05,
                           qvalueCutoff = 0.20, minGSSize = 10, maxGSSize = 500),
                  error = function(x) NULL)
    if (is.null(e) || nrow(as.data.frame(e)) == 0) {
        cat(sprintf("  %-14s 유의 term 없음\n", tag)); return(invisible(NULL))
    }
    d <- as.data.frame(clusterProfiler::simplify(e, cutoff = 0.7, by = "p.adjust"))
    fwrite(d, file.path(P, sprintf("ORA_GO_core287_%s.tsv", tag)), sep = "\t")
    cat(sprintf("  %-14s %d terms (중복 제거 후)\n", tag, nrow(d)))
}

run(core$gene, "all")
run(up$gene,   "up")
run(dn$gene,   "down")
cat("완료 ->", P, "\n")
