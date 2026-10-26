#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────
# STEP 02 — GRCh38 / GENCODE v47 레퍼런스 준비
#   1) genome fasta / GTF / transcript fasta 다운로드
#   2) STAR 인덱스 (sjdbOverhang 150 — 2x151bp read)
#   3) Salmon decoy-aware 인덱스
#   각 산출물이 이미 있으면 건너뜁니다.
# ─────────────────────────────────────────────────────────────
source "$(dirname "$0")/00_config.sh"

BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_REL}"
cd "$REF_DIR"

fetch() {   # $1=url  $2=출력파일(gz)
    if [[ -s "$2" ]]; then log "SKIP  $2 (이미 존재)"; return; fi
    log "GET   $2"
    curl -fsS --retry 5 --retry-delay 10 -o "$2.part" "$1"
    mv "$2.part" "$2"
}

# ── 1) 다운로드 ──────────────────────────────────────────────
fetch "$BASE/GRCh38.primary_assembly.genome.fa.gz"                    "GRCh38.primary_assembly.genome.fa.gz"
fetch "$BASE/gencode.v${GENCODE_REL}.primary_assembly.annotation.gtf.gz" "gencode.v${GENCODE_REL}.primary_assembly.annotation.gtf.gz"
fetch "$BASE/gencode.v${GENCODE_REL}.transcripts.fa.gz"               "gencode.v${GENCODE_REL}.transcripts.fa.gz"

for f in GRCh38.primary_assembly.genome gencode.v${GENCODE_REL}.primary_assembly.annotation gencode.v${GENCODE_REL}.transcripts; do
    ext="fa"; [[ "$f" == *annotation ]] && ext="gtf"
    if [[ ! -s "$f.$ext" ]]; then
        log "UNZIP $f.$ext"
        "$PIGZ_BIN" -dc -p "$THREADS" "$f.$ext.gz" > "$f.$ext"
    fi
done

require_file "$GENOME_FA"; require_file "$GTF"; require_file "$TXFA"
log "genome: $(grep -c '^>' "$GENOME_FA") 개 시퀀스 / GTF: $(grep -vc '^#' "$GTF") 줄"

# ── 2) STAR 인덱스 ───────────────────────────────────────────
if [[ -s "$STAR_INDEX/SA" ]]; then
    log "SKIP  STAR 인덱스 (이미 존재)"
else
    log "BUILD STAR 인덱스 (수십 분 소요)"
    mkdir -p "$STAR_INDEX"
    "$STAR" --runMode genomeGenerate \
        --runThreadN "$THREADS" \
        --genomeDir "$STAR_INDEX" \
        --genomeFastaFiles "$GENOME_FA" \
        --sjdbGTFfile "$GTF" \
        --sjdbOverhang 150 \
        --outTmpDir "$TMP_DIR/star_genome_tmp" \
        --outFileNamePrefix "$REF_DIR/star_index_" \
        ${STAR_GENOME_EXTRA:-}   # [deposit] optional extra args (test: --genomeSAindexNbases 8); default none
    log "STAR 인덱스 완료"
fi

# ── 3) Salmon decoy-aware 인덱스 ─────────────────────────────
if [[ -s "$SALMON_INDEX/info.json" ]]; then
    log "SKIP  Salmon 인덱스 (이미 존재)"
else
    log "BUILD Salmon decoy-aware 인덱스"
    grep '^>' "$GENOME_FA" | cut -d' ' -f1 | sed 's/^>//' > "$REF_DIR/decoys.txt"
    cat "$TXFA" "$GENOME_FA" > "$TMP_DIR/gentrome.fa"
    "$SALMON" index \
        -t "$TMP_DIR/gentrome.fa" \
        -d "$REF_DIR/decoys.txt" \
        -i "$SALMON_INDEX" \
        -k 31 -p "$THREADS" --gencode
    rm -f "$TMP_DIR/gentrome.fa"
    log "Salmon 인덱스 완료"
fi

log "════════ STEP 02 완료 ════════"
du -sh "$STAR_INDEX" "$SALMON_INDEX"
