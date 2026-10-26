#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────
# GTF → BED12 변환 (RSeQC 입력용)
#
# RSeQC의 geneBody_coverage.py, read_distribution.py,
# junction_annotation.py, junction_saturation.py, inner_distance.py 는
# 모두 BED12 형식의 유전자 모델을 요구합니다.
#
# 추가로 geneBody_coverage 전용 축소판(무작위 2000 transcript)도 만듭니다.
# 전체 14만 transcript로 돌리면 샘플당 수 시간이 걸리기 때문입니다.
# ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
    library(GenomicFeatures)
    library(rtracklayer)
    # Bioconductor 3.20부터 makeTxDbFromGFF()가 GenomicFeatures에서
    # txdbmaker로 이전되었고, GenomicFeatures >= 1.61.1에서는 defunct입니다.
    library(txdbmaker)
})

args <- commandArgs(trailingOnly = TRUE)
gtf_file <- args[1]
out_bed  <- args[2]
out_sub  <- args[3]   # 축소판 경로

stopifnot(file.exists(gtf_file))

message("TxDb 생성 중 (수 분 소요)...")
txdb <- txdbmaker::makeTxDbFromGFF(gtf_file, format = "gtf")

message("transcript별 exon 추출...")
ex <- exonsBy(txdb, by = "tx", use.names = TRUE)

message("BED12 변환...")
bed <- asBED(ex)
# asBED는 name 컬럼을 transcript ID로 채워줍니다
export(bed, out_bed, format = "BED")
message("저장: ", out_bed, "  (", length(bed), " transcripts)")

# ── geneBody_coverage 전용 축소판 ─────────────────────────────
# 짧은 transcript는 커버리지 곡선을 왜곡하므로 1000 bp 이상만 사용
w <- sum(width(ex))
keep <- names(w)[w >= 1000]
set.seed(42)
sel <- sample(keep, min(2000, length(keep)))
bed_sub <- bed[bed$name %in% sel]
export(bed_sub, out_sub, format = "BED")
message("저장: ", out_sub, "  (", length(bed_sub), " transcripts, geneBody_coverage용)")

message("완료")
