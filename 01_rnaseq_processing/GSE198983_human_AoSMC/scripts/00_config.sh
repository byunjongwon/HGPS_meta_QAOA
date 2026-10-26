#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────
# 공통 설정 — 모든 스크립트가 이 파일을 source 합니다.
# 프로젝트: GSE198983 / PRJNA817614
#   Human primary aortic SMC, lentiviral Flag-progerin vs Flag-only control
#   8 샘플 (4 vs 4), paired-end 2x151bp, Illumina NovaSeq 6000
#   논문: Coll-Bonfill et al. (PMID 36482259), Gonzalo lab
# ─────────────────────────────────────────────────────────────
set -euo pipefail

# ── 경로 ──────────────────────────────────────────────────────
# [deposit] PROJ default = dataset dir (parent of this scripts/ dir); still overridable by PROJ
export PROJ="${PROJ:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export FASTQ_DIR="$PROJ/raw_data/fastq"
export GEO_DIR="$PROJ/raw_data/geo"
export REF_DIR="$PROJ/ref"
export META_DIR="$PROJ/meta"
export LOG_DIR="$PROJ/logs"
export RESULT_DIR="$PROJ/results"
export TMP_DIR="$PROJ/tmp"

# ── 레퍼런스 (GRCh38 / GENCODE v47) ───────────────────────────
export GENCODE_REL="47"
export GENOME_FA="$REF_DIR/GRCh38.primary_assembly.genome.fa"
export GTF="$REF_DIR/gencode.v${GENCODE_REL}.primary_assembly.annotation.gtf"
export TXFA="$REF_DIR/gencode.v${GENCODE_REL}.transcripts.fa"
# [deposit] the deposit ships only the .gz annotation/transcripts; 16_isoformswitch.R falls back to
#           "$GTF.gz"/"$TXFA.gz" when the uncompressed files are absent (02_reference.sh creates them).
export STAR_INDEX="$REF_DIR/star_index"
export SALMON_INDEX="$REF_DIR/salmon_index"
# [deposit] extra STAR genomeGenerate args (e.g. "--genomeSAindexNbases 8" for the tiny test genome);
#           default empty = original behaviour
export STAR_GENOME_EXTRA="${STAR_GENOME_EXTRA:-}"
export BED12="$REF_DIR/gencode.v${GENCODE_REL}.bed12"           # RSeQC용
export BED12_HK="$REF_DIR/gencode.v${GENCODE_REL}.genebody2000.bed12"  # geneBody_coverage용 서브셋

export SAMPLES="$META_DIR/samples.tsv"

# ── 리소스 ────────────────────────────────────────────────────
export THREADS="${THREADS:-18}"          # 전체 20코어 중 18 사용
export SORT_MEM="${SORT_MEM:-4G}"        # samtools sort 스레드당 메모리

# ── 실행 파일 (env가 여러 개로 나뉘어 있어 절대경로로 고정) ────
CONDA_ENVS="${CONDA_ENVS:-$HOME/miniforge3/envs}"   # [deposit] overridable
export STAR="$CONDA_ENVS/lmna_rmats/bin/STAR"
export RMATS="$CONDA_ENVS/lmna_rmats/bin/rmats.py"
export SAMTOOLS="$CONDA_ENVS/lmna_rmats/bin/samtools"

export FASTQC="$CONDA_ENVS/rnaseq/bin/fastqc"
export MULTIQC="$CONDA_ENVS/rnaseq/bin/multiqc"
export FEATURECOUNTS="$CONDA_ENVS/rnaseq/bin/featureCounts"
export SALMON="$CONDA_ENVS/rnaseq/bin/salmon"
export PICARD="$CONDA_ENVS/rnaseq/bin/picard"
export PIGZ_BIN="$CONDA_ENVS/rnaseq/bin/pigz"
export RSEQC_BIN="$CONDA_ENVS/rnaseq/bin"          # infer_experiment.py 등

export FASTP="$CONDA_ENVS/lmna_tools/bin/fastp"
export QUALIMAP="$CONDA_ENVS/lmna_tools/bin/qualimap"
export BEDTOOLS="$CONDA_ENVS/lmna_tools/bin/bedtools"
export STRINGTIE="$CONDA_ENVS/lmna_tools/bin/stringtie"

export RSCRIPT="$CONDA_ENVS/lmna_r/bin/Rscript"

# BiocManager가 이 네트워크의 TLS 프록시를 통과하도록
export CURL_CA_BUNDLE="${CURL_CA_BUNDLE:-/etc/ssl/certs/ca-certificates.crt}"   # [deposit] overridable

# ── 헬퍼 ──────────────────────────────────────────────────────
mkdir -p "$LOG_DIR" "$RESULT_DIR" "$TMP_DIR"

log() { echo "[$(date '+%F %T')] $*"; }

# samples.tsv에서 헤더를 뺀 본문만 출력
sample_rows() { tail -n +2 "$SAMPLES"; }

# 샘플명 목록
sample_names() { sample_rows | cut -f1; }

# 필수 파일 존재 확인
require_file() {
    [[ -s "$1" ]] || { echo "ERROR: 파일이 없습니다: $1" >&2; exit 1; }
}
