#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────
# 27. QUBO 유전자 주석 추출
#
# 표에 들어갈 '유전자 이름'과 '기능'을 기억에 의존해 쓰지 않고
# org.Hs.eg.db 에서 공식 명칭과 GO:BP term 을 직접 뽑는다.
# 표의 서술은 여기서 나온 term 을 근거로만 작성한다.
# ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
    library(org.Hs.eg.db); library(AnnotationDbi); library(data.table)
})

# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
.this <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
DEPOSIT <- Sys.getenv("PROGERIA_DEPOSIT", normalizePath(file.path(dirname(.this), "..", "..")))
OUT <- file.path(DEPOSIT, "03_meta_analysis", "results")  # [deposit] was .../progeria_gse_원본/meta/results

GENES <- unique(c(
    "ACVRL1","AGTR1","CAMK2B","FGFR3","GAS6","IL1A","LYN","MAFB","OGN","PDE4D",
    "PECAM1","PLCB1","PTGES","PTN","SFRP1","SNAP25","TENT5C","WNT2","WNT5A",
    "GCNT4","TSHZ2","CSGALNACT1","NDRG2","RAB39B","PPL","ST6GALNAC3","SELENOP",
    "G0S2","THRB"))

ann <- AnnotationDbi::select(org.Hs.eg.db, keys = GENES, keytype = "SYMBOL",
                             columns = c("GENENAME", "ENTREZID"))
ann <- unique(as.data.table(ann))
cat(sprintf("유전자 %d개 중 주석 %d개\n", length(GENES), nrow(ann)))

# GO:BP term 을 유전자마다 최대 6개까지
go <- AnnotationDbi::select(org.Hs.eg.db, keys = GENES, keytype = "SYMBOL",
                            columns = c("GO", "ONTOLOGY"))
go <- as.data.table(go)[ONTOLOGY == "BP" & !is.na(GO)]
suppressPackageStartupMessages(library(GO.db))
go[, term := AnnotationDbi::mapIds(GO.db, keys = GO, keytype = "GOID",
                                   column = "TERM")]
go <- unique(go[, .(SYMBOL, term)])
agg <- go[, .(go_bp = paste(head(term, 6), collapse = " | ")), by = SYMBOL]

res <- merge(ann[, .(SYMBOL, GENENAME)], agg, by = "SYMBOL", all.x = TRUE)
fwrite(res, file.path(OUT, "qubo_gene_annotation.tsv"), sep = "\t")
cat("저장 ->", file.path(OUT, "qubo_gene_annotation.tsv"), "\n")
print(res[1:5])
