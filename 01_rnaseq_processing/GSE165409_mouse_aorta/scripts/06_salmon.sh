#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────
# STEP 6. Salmon — transcript 수준 정량 (isoform/DTU 분석용)
#
# selective alignment 모드 사용 (decoy = 전체 게놈).
# 이 결과가 IsoformSwitchAnalyzeR / tximport 입력이 됩니다.
#
# 산출:
#   salmon_index/           인덱스 (약 15 G, 빌드 30~40분)
#   quant/{sample}/quant.sf transcript별 TPM/NumReads
#   → isoform switch plot, DTU volcano, isoform usage bar 재료
#
# 예상 시간: 인덱스 40분 + 정량 18샘플 x 5분 = 2시간
# ─────────────────────────────────────────────────────────────
source "$(dirname "$0")/00_config.sh"

OUT="$RESULT_DIR/06_salmon"
IDX="$REF_DIR/salmon_index"
TX_FA="$REF_DIR/gencode.vM25.transcripts.fa"
mkdir -p "$OUT/quant"

INPUT_DIR="${INPUT_DIR:-$RESULT_DIR/02_trimmed}"
SUFFIX="${SUFFIX:-.trim.fastq.gz}"

# ── transcriptome fasta 준비 ──────────────────────────────────
if [[ ! -s "$TX_FA" ]]; then
    log "GENCODE vM25 transcript fasta 다운로드"
    /usr/bin/curl -L --max-time 900 -o "${TX_FA}.gz" \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.transcripts.fa.gz"
    gunzip -f "${TX_FA}.gz"
fi
require_file "$TX_FA"

# ── decoy-aware 인덱스 ────────────────────────────────────────
if [[ ! -s "$IDX/info.json" ]]; then
    log "salmon 인덱스 빌드 (decoy = 전체 게놈, 30~40분)"
    grep "^>" "$GENOME_FA" | cut -d " " -f1 | sed 's/>//g' > "$REF_DIR/decoys.txt"
    cat "$TX_FA" "$GENOME_FA" > "$REF_DIR/gentrome.fa"
    "$SALMON" index -t "$REF_DIR/gentrome.fa" -d "$REF_DIR/decoys.txt" \
        -p "$THREADS" -i "$IDX" --gencode 2>&1 | tee -a "$LOG_DIR/06_salmon.log"
    rm -f "$REF_DIR/gentrome.fa"     # 3 G 절약
fi

# ── 정량 ──────────────────────────────────────────────────────
log "salmon 정량 시작 (libType=A 자동 판정)"
sample_rows | while IFS=$'\t' read -r sample gsm group genotype age_group counts_col run1 run2; do
    if [[ -s "$OUT/quant/$sample/quant.sf" ]]; then
        log "  $sample — 이미 완료, 건너뜀"
        continue
    fi
    log "  $sample 정량 중"
    "$SALMON" quant -i "$IDX" -l A \
        -1 "$INPUT_DIR/${run1}_1${SUFFIX}" "$INPUT_DIR/${run2}_1${SUFFIX}" \
        -2 "$INPUT_DIR/${run1}_2${SUFFIX}" "$INPUT_DIR/${run2}_2${SUFFIX}" \
        -p "$THREADS" \
        --validateMappings --seqBias --gcBias --posBias \
        --numGibbsSamples 20 \
        -o "$OUT/quant/$sample" 2>&1 | tail -3 | tee -a "$LOG_DIR/06_salmon.log"
done

log "매핑률 요약"
{
    echo -e "sample\tmapping_rate"
    for d in "$OUT"/quant/*/; do
        s=$(basename "$d")
        r=$(grep -oP '"percent_mapped":\s*\K[0-9.]+' "$d/aux_info/meta_info.json" 2>/dev/null || echo NA)
        echo -e "$s\t$r"
    done
} | tee "$OUT/mapping_rates.tsv" | column -t

log "MultiQC"
"$MULTIQC" "$OUT" -o "$OUT" -n multiqc_salmon --force 2>&1 | tee -a "$LOG_DIR/06_salmon.log"

log "완료 → $OUT"
