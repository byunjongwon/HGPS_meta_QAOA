#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────
# STEP 5. featureCounts — 유전자 수준 counts 행렬
#
# 산출:
#   counts_raw.txt          featureCounts 원본
#   counts_matrix.tsv       유전자 x 18샘플 정리본 (DESeq2 입력)
#   counts_summary.tsv      assigned/unassigned 통계 → QC bar plot 재료
#   author_vs_ours.tsv      저자 counts와의 비교 (파이프라인 검증)
#
# ⚠️ STRANDNESS: 04 스크립트의 판정 결과를 반드시 넣어야 합니다.
#    잘못 넣으면 counts가 절반 이하로 떨어집니다.
#    환경변수로 넘기거나 04가 만든 파일에서 자동으로 읽습니다.
#
# 예상 시간: 10~20분
# ─────────────────────────────────────────────────────────────
source "$(dirname "$0")/00_config.sh"

BAM_DIR="$RESULT_DIR/03_star"
OUT="$RESULT_DIR/05_counts"
mkdir -p "$OUT"

# ── strandedness 자동 판정 ────────────────────────────────────
STRAND_FILE="$RESULT_DIR/04_qc_align/strandedness_summary.tsv"
if [[ -z "${STRANDNESS:-}" && -s "$STRAND_FILE" ]]; then
    # 모든 샘플의 fwd/rev 평균으로 판정
    # printf에 개행을 반드시 넣어야 합니다. 개행이 없으면 read가 EOF를 만나
    # 값을 정상적으로 읽고도 exit 1을 반환해 set -e가 스크립트를 중단시킵니다.
    read -r AVG_FWD AVG_REV < <(awk -F'\t' 'NR>1 && $3!="NA" {f+=$3; r+=$4; n++} END{printf "%.3f %.3f\n", f/n, r/n}' "$STRAND_FILE")
    if   (( $(echo "$AVG_REV > 0.75" | bc -l) )); then STRANDNESS=2
    elif (( $(echo "$AVG_FWD > 0.75" | bc -l) )); then STRANDNESS=1
    else STRANDNESS=0; fi
    log "strandedness 자동 판정: fwd=$AVG_FWD rev=$AVG_REV → featureCounts -s $STRANDNESS"
fi
STRANDNESS="${STRANDNESS:-0}"
echo "$STRANDNESS" > "$OUT/strandness_used.txt"

BAMS=("$BAM_DIR"/*.Aligned.sortedByCoord.out.bam)
log "featureCounts 시작 — BAM ${#BAMS[@]}개, -s $STRANDNESS, paired-end"

"$FEATURECOUNTS" \
    -T "$THREADS" \
    -p --countReadPairs \
    -s "$STRANDNESS" \
    -t exon -g gene_id \
    -a "$GTF" \
    -o "$OUT/counts_raw.txt" \
    "${BAMS[@]}" 2>&1 | tee -a "$LOG_DIR/05_featurecounts.log"

# ── 정리본 생성 (샘플명만 컬럼으로) ───────────────────────────
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

# featureCounts summary → QC 그림 재료
cp "$OUT/counts_raw.txt.summary" "$OUT/counts_summary.tsv"

log "행렬 크기: $(( $(wc -l < "$OUT/counts_matrix.tsv") - 1 )) 유전자 x $(( $(head -1 "$OUT/counts_matrix.tsv" | awk -F'\t' '{print NF}') - 1 )) 샘플"

# ── 저자 counts와 비교 (파이프라인 검증) ──────────────────────
log "저자 counts와 상관 비교"
"$RSCRIPT" - <<'RS' 2>&1 | tee -a "$LOG_DIR/05_featurecounts.log"
suppressPackageStartupMessages(library(data.table))
proj <- Sys.getenv("RESULT_DIR")
geo  <- Sys.getenv("GEO_DIR")
meta <- fread(Sys.getenv("SAMPLES"))

ours <- fread(file.path(proj, "05_counts/counts_matrix.tsv"))
ours[, gene_id := sub("\\..*$", "", gene_id)]          # GENCODE 버전 접미사 제거
ours <- ours[, lapply(.SD, sum), by = gene_id]          # PAR_Y 중복 합치기

auth <- fread(file.path(geo, "GSE165409_LMNA_rawCounts_18samp.txt.gz"))
setnames(auth, 1, "gene_id")

# 샘플명 매칭: 우리 컬럼(sample) → 저자 컬럼(counts_col)
map <- setNames(meta$counts_col, meta$sample)
common_genes <- intersect(ours$gene_id, auth$gene_id)
cat(sprintf("공통 유전자: %d개 (ours %d / author %d)\n",
            length(common_genes), nrow(ours), nrow(auth)))

o <- ours[match(common_genes, gene_id)]
a <- auth[match(common_genes, gene_id)]

res <- data.table(sample = names(map), author_col = unname(map))
res[, spearman := mapply(function(s, ac) {
        if (!(s %in% names(o)) || !(ac %in% names(a))) return(NA_real_)
        suppressWarnings(cor(o[[s]], a[[ac]], method = "spearman"))
    }, sample, author_col)]
res[, pearson_log := mapply(function(s, ac) {
        if (!(s %in% names(o)) || !(ac %in% names(a))) return(NA_real_)
        suppressWarnings(cor(log1p(o[[s]]), log1p(a[[ac]])))
    }, sample, author_col)]

fwrite(res, file.path(proj, "05_counts/author_vs_ours.tsv"), sep = "\t")
print(res)
cat(sprintf("\n평균 Spearman: %.4f\n", mean(res$spearman, na.rm = TRUE)))
cat("→ 0.95 이상이면 파이프라인 설정(GTF, strandedness)이 저자와 일치합니다.\n")
RS

log "완료 → $OUT"
