#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────
# STEP 07 — featureCounts 유전자 수준 counts 행렬
#
# 산출:
#   counts_raw.txt        featureCounts 원본
#   counts_matrix.tsv     gene x 8샘플 (DESeq2 입력)
#   counts_summary.tsv    assigned/unassigned 통계
#   star_vs_fc.tsv        STAR ReadsPerGene와의 교차검증
#
# ⚠️ STRANDNESS는 STEP 06 판정 결과를 자동으로 읽습니다.
# ─────────────────────────────────────────────────────────────
source "$(dirname "$0")/00_config.sh"

BAM_DIR="$RESULT_DIR/05_star"
OUT="$RESULT_DIR/07_counts"
mkdir -p "$OUT"

# ── strandedness 자동 판정 ────────────────────────────────────
STRAND_FILE="$RESULT_DIR/06_qc_align/strandedness_summary.tsv"
if [[ -z "${STRANDNESS:-}" && -s "$STRAND_FILE" ]]; then
    read -r AVG_FWD AVG_REV < <(awk -F'\t' 'NR>1 && $3!="NA" {f+=$3; r+=$4; n++} END{printf "%.3f %.3f\n", f/n, r/n}' "$STRAND_FILE")
    if   (( $(echo "$AVG_REV > 0.75" | bc -l) )); then STRANDNESS=2
    elif (( $(echo "$AVG_FWD > 0.75" | bc -l) )); then STRANDNESS=1
    else STRANDNESS=0; fi
    log "strandedness 자동 판정: fwd=$AVG_FWD rev=$AVG_REV → featureCounts -s $STRANDNESS"
fi
STRANDNESS="${STRANDNESS:-0}"
echo "$STRANDNESS" > "$OUT/strandness_used.txt"

BAMS=("$BAM_DIR"/*.Aligned.sortedByCoord.out.bam)
log "featureCounts — BAM ${#BAMS[@]}개, -s $STRANDNESS, paired-end"

"$FEATURECOUNTS" \
    -T "$THREADS" \
    -p --countReadPairs \
    -s "$STRANDNESS" \
    -t exon -g gene_id \
    -a "$GTF" \
    -o "$OUT/counts_raw.txt" \
    "${BAMS[@]}" 2>&1 | tee -a "$LOG_DIR/07_featurecounts.log"

log "counts 행렬 정리"
awk -F'\t' 'NR==2{
        printf "gene_id"
        for(i=7;i<=NF;i++){
            n=split($i,p,"/"); s=p[n]; sub(/\.Aligned\.sortedByCoord\.out\.bam$/,"",s)
            printf "\t%s", s
        }
        printf "\n"; next
    }
    NR>2{
        printf "%s", $1
        for(i=7;i<=NF;i++) printf "\t%s", $i
        printf "\n"
    }' "$OUT/counts_raw.txt" > "$OUT/counts_matrix.tsv"

# 유전자 길이 (TPM 계산용)
awk -F'\t' 'NR>2{print $1"\t"$6}' "$OUT/counts_raw.txt" > "$OUT/gene_lengths.tsv"

cp "$OUT/counts_raw.txt.summary" "$OUT/counts_summary.tsv"

log "행렬: $(( $(wc -l < "$OUT/counts_matrix.tsv") - 1 )) 유전자 x $(( $(head -1 "$OUT/counts_matrix.tsv" | awk -F'\t' '{print NF}') - 1 )) 샘플"

# ── STAR GeneCounts와 교차검증 ────────────────────────────────
log "STAR ReadsPerGene와 상관 비교 (파이프라인 검증)"
"$RSCRIPT" - <<'RS' 2>&1 | tee -a "$LOG_DIR/07_featurecounts.log"
suppressPackageStartupMessages(library(data.table))
res_dir <- Sys.getenv("RESULT_DIR")
meta <- fread(Sys.getenv("SAMPLES"))
strand <- as.integer(readLines(file.path(res_dir, "07_counts/strandness_used.txt"))[1])
# STAR ReadsPerGene 컬럼: 2=unstranded, 3=forward(s1), 4=reverse(s2)
col <- c("2" = 3L, "1" = 3L, "0" = 2L)[as.character(strand)]
if (strand == 2L) col <- 4L

fc <- fread(file.path(res_dir, "07_counts/counts_matrix.tsv"))
out <- data.table(sample = character(), spearman = numeric(), pearson_log = numeric())
for (s in meta$sample) {
    f <- file.path(res_dir, "05_star", paste0(s, ".ReadsPerGene.out.tab"))
    if (!file.exists(f)) next
    st <- fread(f, header = FALSE)
    st <- st[!V1 %like% "^N_"]
    m <- merge(fc[, .(gene_id, fc = get(s))], st[, .(gene_id = V1, star = get(paste0("V", col)))],
               by = "gene_id")
    out <- rbind(out, data.table(sample = s,
                                 spearman = cor(m$fc, m$star, method = "spearman"),
                                 pearson_log = cor(log1p(m$fc), log1p(m$star))))
}
fwrite(out, file.path(res_dir, "07_counts/star_vs_fc.tsv"), sep = "\t")
print(out)
cat(sprintf("\n평균 Spearman: %.4f  → 0.98 이상이면 정상\n", mean(out$spearman, na.rm = TRUE)))
RS

log "완료 → $OUT"
