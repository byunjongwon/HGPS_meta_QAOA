#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────
# 공통 설정 — 모든 스크립트가 이 파일을 source 합니다.
# 프로젝트: GSE165409 / SRP303066
#   LMNA G609G (progeria) vs WT young vs WT old, 마우스 하행대동맥
#   18 샘플 x 2 run, paired-end 2x150bp, HiSeq 4000
# ─────────────────────────────────────────────────────────────
set -euo pipefail

# ── 경로 ──────────────────────────────────────────────────────
# [deposit] default PROJ = this dataset dir (scripts/..), still overridable by PROJ
export PROJ="${PROJ:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export FASTQ_DIR="$PROJ/raw_data/fastq"
export GEO_DIR="$PROJ/raw_data/geo"
export REF_DIR="$PROJ/ref"
export META_DIR="$PROJ/meta"
export LOG_DIR="$PROJ/logs"
export RESULT_DIR="$PROJ/results"

export GENOME_FA="$REF_DIR/GRCm38.primary_assembly.genome.fa"
export GTF="$REF_DIR/gencode.vM25.primary_assembly.annotation.gtf"
export STAR_INDEX="$REF_DIR/star_index"
export BED12="$REF_DIR/gencode.vM25.bed12"          # RSeQC용 (03b에서 생성)

export SAMPLES="$META_DIR/samples.tsv"

# ── 리소스 ────────────────────────────────────────────────────
export THREADS="${THREADS:-18}"          # 전체 20코어 중 18 사용
export SORT_MEM="${SORT_MEM:-4G}"        # samtools sort 스레드당 메모리

# ── 실행 파일 (env가 여러 개로 나뉘어 있어 절대경로로 고정) ────
# [deposit] conda env root overridable by CONDA_ENVS
export CONDA_ENVS="${CONDA_ENVS:-$HOME/miniforge3/envs}"
export STAR="$CONDA_ENVS/lmna_rmats/bin/STAR"
export RMATS="$CONDA_ENVS/lmna_rmats/bin/rmats.py"
export SAMTOOLS="$CONDA_ENVS/lmna_rmats/bin/samtools"

export FASTQC="$CONDA_ENVS/rnaseq/bin/fastqc"
export MULTIQC="$CONDA_ENVS/rnaseq/bin/multiqc"
export FEATURECOUNTS="$CONDA_ENVS/rnaseq/bin/featureCounts"
export SALMON="$CONDA_ENVS/rnaseq/bin/salmon"
export PICARD="$CONDA_ENVS/rnaseq/bin/picard"
export PIGZ="$CONDA_ENVS/rnaseq/bin/pigz"
export RSEQC_BIN="$CONDA_ENVS/rnaseq/bin"          # infer_experiment.py 등

export FASTP="$CONDA_ENVS/lmna_tools/bin/fastp"
export QUALIMAP="$CONDA_ENVS/lmna_tools/bin/qualimap"
export BEDTOOLS="$CONDA_ENVS/lmna_tools/bin/bedtools"
export STRINGTIE="$CONDA_ENVS/lmna_tools/bin/stringtie"

export RSCRIPT="$CONDA_ENVS/lmna_r/bin/Rscript"

# BiocManager가 이 네트워크의 TLS 프록시를 통과하도록
export CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt

# ── 헬퍼 ──────────────────────────────────────────────────────
mkdir -p "$LOG_DIR" "$RESULT_DIR"

log() { echo "[$(date '+%F %T')] $*"; }

# samples.tsv에서 헤더를 뺀 본문만 출력
sample_rows() { tail -n +2 "$SAMPLES"; }

# 샘플 하나의 R1/R2 파일 목록을 콤마로 이어 출력 (STAR에 그대로 전달)
#   $1 = run1, $2 = run2, $3 = mate(1|2), $4 = 디렉토리
mate_list() {
    local r1="$1" r2="$2" mate="$3" dir="$4"
    echo "$dir/${r1}_${mate}.fastq.gz,$dir/${r2}_${mate}.fastq.gz"
}

# 필수 파일 존재 확인
require_file() {
    [[ -s "$1" ]] || { echo "ERROR: 파일이 없습니다: $1" >&2; exit 1; }
}
