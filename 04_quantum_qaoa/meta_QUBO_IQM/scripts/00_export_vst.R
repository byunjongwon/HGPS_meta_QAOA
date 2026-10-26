#!/usr/bin/env Rscript
# GSE306264 VST 행렬을 양자 파이프라인이 쓸 수 있게 gene symbol 기준으로 내보낸다.
suppressPackageStartupMessages(library(DESeq2))
# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
.this <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
DEPOSIT <- Sys.getenv("PROGERIA_DEPOSIT", normalizePath(file.path(dirname(.this), "..", "..", "..")))
GSE_RES <- file.path(DEPOSIT, "02_per_dataset_DE", "GSE306264", "results")  # [deposit] was ~/car_project/progeria_meta/GSE306264/results
o <- readRDS(file.path(GSE_RES, "GSE306264_objects.rds"))
vm <- assay(o$vsd)
write.table(vm, file.path(GSE_RES, "vst_for_quantum.tsv"),
            sep = "\t", quote = FALSE)
cat(sprintf("VST 내보내기: %d 유전자 x %d 세포주\n", nrow(vm), ncol(vm)))
