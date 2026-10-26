#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────
# [deposit] created for the code deposit (the original project has no
# reference-building script; the full STAR index was built once by hand, see
# the command below). Builds the reference files for the TEST SUBSET into
# $PROJ/ref under the file names that ../scripts/00_config.sh expects:
#
#   ref/GRCm38.primary_assembly.genome.fa           <- test/ref/mini_genome.fa.gz
#   ref/gencode.vM25.primary_assembly.annotation.gtf <- test/ref/mini_annotation.gtf.gz
#   ref/gencode.vM25.transcripts.fa                  <- test/ref/mini_transcripts.fa.gz
#   ref/star_index/                                  STAR genomeGenerate (small --genomeSAindexNbases)
#   ref/gencode.vM25.bed12, ref/gencode.vM25.genebody2000.bed12  via scripts/03b_make_bed12.R
# (the salmon index is built by scripts/06_salmon.sh itself)
#
# Full-size equivalent used for the real analysis (STAR 2.7.11b, from logs/star_index.log):
#   STAR --runMode genomeGenerate --runThreadN 18 --genomeDir ref/star_index \
#        --genomeFastaFiles ref/GRCm38.primary_assembly.genome.fa \
#        --sjdbGTFfile ref/gencode.vM25.primary_assembly.annotation.gtf --sjdbOverhang 149 \
#        --limitGenomeGenerateRAM 100000000000
#
# Usage: PROJ=/path/to/workdir bash build_test_reference.sh      (normally called by run_test.sh)
# ─────────────────────────────────────────────────────────────
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$(cd "$HERE/../scripts" && pwd)"
: "${PROJ:?set PROJ to the test work directory}"
source "$SCRIPTS/00_config.sh"
mkdir -p "$REF_DIR"

log "test reference -> $REF_DIR"
zcat "$HERE/ref/mini_genome.fa.gz"      > "$GENOME_FA"
zcat "$HERE/ref/mini_annotation.gtf.gz" > "$GTF"
zcat "$HERE/ref/mini_transcripts.fa.gz" > "$REF_DIR/gencode.vM25.transcripts.fa"

# genomeSAindexNbases = min(14, log2(genome length)/2 - 1)  (STAR manual, small genomes)
GLEN=$(grep -v '^>' "$GENOME_FA" | tr -d '\n' | wc -c)
SAIDX=$(awk -v g="$GLEN" 'BEGIN { n = int(log(g) / log(2) / 2 - 1); if (n > 14) n = 14; print n }')
log "STAR genomeGenerate (genome ${GLEN} bp, --genomeSAindexNbases $SAIDX)"
mkdir -p "$STAR_INDEX"
"$STAR" --runMode genomeGenerate --runThreadN "$THREADS" \
    --genomeDir "$STAR_INDEX" --genomeFastaFiles "$GENOME_FA" \
    --sjdbGTFfile "$GTF" --sjdbOverhang 149 \
    --genomeSAindexNbases "$SAIDX" \
    --outFileNamePrefix "$LOG_DIR/starIdx_" > "$LOG_DIR/star_index.log" 2>&1
require_file "$STAR_INDEX/SA"

log "BED12 (03b_make_bed12.R)"
"$RSCRIPT" "$SCRIPTS/03b_make_bed12.R" "$GTF" "$BED12" "$REF_DIR/gencode.vM25.genebody2000.bed12" \
    > "$LOG_DIR/03b_make_bed12.log" 2>&1
require_file "$BED12"
log "done"
