#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────
# [deposit] created for the code deposit — documents and reproduces how the
# tiny test subset in this folder (test/) was made. It is NOT part of the
# original analysis.
#
# Source (read-only): the original STAR BAMs of the full run
#   $SRC_PROJ/results/03_star/<sample>.Aligned.sortedByCoord.out.bam   (18 samples)
#   $SRC_PROJ/ref/GRCm38.primary_assembly.genome.fa                     (GENCODE M25 / GRCm38)
#   GENCODE vM25 GTF + transcripts (.gz copies shipped in ../ref/)
#
# What it makes (in OUT_DIR, default = the folder of this script):
#   raw_data/fastq/<SRR>_{1,2}.fastq.gz   2 runs x 18 samples, a few thousand read pairs each
#   ref/mini_genome.fa.gz                 3 windows of GRCm38, contigs keep the names chr3/chr5/chr6
#   ref/mini_annotation.gtf.gz            GENCODE vM25 genes fully inside the windows, coordinates
#                                         shifted to the window (pos - window_start + 1)
#   ref/mini_transcripts.fa.gz            GENCODE vM25 transcripts of those genes (headers unchanged)
#   meta/samples.tsv                      identical to ../meta/samples.tsv
#
# Windows (GRCm38, 1-based inclusive):
#   chr3:87,800,001-89,000,000   1.2 Mb around Lmna (chr3:88,480,147-88,509,956, minus strand);
#                                contains the progerin A5SS of Lmna exon 11 (chr3:88,482,345-88,482,598)
#   chr5:142,880,001-142,930,000 Actb (housekeeping)
#   chr6:125,140,001-125,190,000 Gapdh (housekeeping)
#
# Read selection per sample (properly paired, primary alignments only: -f 2 -F 0x904):
#   pairs are sub-sampled per region by read-name hash (samtools --subsample, seed 42) to about
#     Lmna exon 10-12 core (chr3:88,481,500-88,483,700) : 600 pairs (so rMATS/sashimi see the
#                                                         progerin junctions)
#     whole chr3 window                                 : 2500 pairs
#     Actb window / Gapdh window                        : 300 / 200 pairs
#   then both mates are written back to FASTQ (samtools collate | samtools fastq; reverse-strand
#   reads are reverse-complemented back to the sequencing orientation), and split into the two SRA
#   runs of the sample by the run accession that STAR kept in the read name (SRRxxxx.N).
#   NOTE: the reads in the BAMs are already fastp-trimmed, so the test FASTQs are trimmed reads.
#
# Usage:  bash make_test_subset.sh [OUT_DIR]
#   env: SRC_PROJ (default: ~/lmna_rnaseq, the original project dir), CONDA_ENVS, SEED (42)
# ─────────────────────────────────────────────────────────────
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATASET="$(cd "$HERE/.." && pwd)"
OUT_DIR="${1:-$HERE}"
SRC_PROJ="${SRC_PROJ:?set SRC_PROJ to the full-size original run (STAR BAMs + GRCm38 reference)}"   # [deposit] originally ~/lmna_rnaseq
CONDA_ENVS="${CONDA_ENVS:-$HOME/miniforge3/envs}"
SAMTOOLS="$CONDA_ENVS/lmna_rmats/bin/samtools"
PY="${PYTHON:-python3}"
SEED="${SEED:-42}"

BAM_DIR="$SRC_PROJ/results/03_star"
GENOME_FA="$SRC_PROJ/ref/GRCm38.primary_assembly.genome.fa"
GTF_GZ="$DATASET/ref/gencode.vM25.primary_assembly.annotation.gtf.gz"
TX_GZ="$DATASET/ref/gencode.vM25.transcripts.fa.gz"
SAMPLES="$DATASET/meta/samples.tsv"

# contig  start(1-based)  end
WINDOWS="chr3 87800001 89000000
chr5 142880001 142930000
chr6 125140001 125190000"
# region  target_pairs
TARGETS="chr3:88481500-88483700 600
chr3:87800001-89000000 2500
chr5:142880001-142930000 300
chr6:125140001-125190000 200"
FULL_REGIONS="chr3:87800001-89000000 chr5:142880001-142930000 chr6:125140001-125190000"

for f in "$GENOME_FA" "$GTF_GZ" "$TX_GZ" "$SAMPLES" "$SAMTOOLS"; do
    [[ -s "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }
done

mkdir -p "$OUT_DIR"/{raw_data/fastq,ref,meta}
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
cp "$SAMPLES" "$OUT_DIR/meta/samples.tsv"

# ── 1. mini genome (window sequences; contig names kept) ─────
echo "[subset] mini genome"
"$PY" - "$GENOME_FA" "$TMP/mini_genome.fa" "$WINDOWS" <<'PY'
import sys
fa, out, wins = sys.argv[1], sys.argv[2], sys.argv[3]
win = {}
for line in wins.strip().splitlines():
    c, s, e = line.split(); win[c] = (int(s), int(e))
seqs, cur, buf = {}, None, []
with open(fa) as fh:
    for line in fh:
        if line.startswith(">"):
            if cur in win: seqs[cur] = "".join(buf)
            cur = line[1:].split()[0]; buf = []
        elif cur in win:
            buf.append(line.strip())
    if cur in win: seqs[cur] = "".join(buf)
with open(out, "w") as o:
    for c, (s, e) in win.items():
        sub = seqs[c][s - 1:e]
        o.write(f">{c}\n")
        for i in range(0, len(sub), 60):
            o.write(sub[i:i + 60] + "\n")
PY

# ── 2. mini GTF (genes completely inside a window; coordinates shifted) ──
echo "[subset] mini GTF"
zcat "$GTF_GZ" | awk -F'\t' -v OFS='\t' -v W="$WINDOWS" '
    BEGIN { n = split(W, L, "\n"); for (i = 1; i <= n; i++) { split(L[i], f, " "); S[f[1]] = f[2]; E[f[1]] = f[3] } }
    /^#/ { print; next }
    ($1 in S) && $3 == "gene" && $4 >= S[$1] && $5 <= E[$1] {
        match($9, /gene_id "[^"]+"/); keep[substr($9, RSTART, RLENGTH)] = 1 }
    ($1 in S) {
        match($9, /gene_id "[^"]+"/); g = substr($9, RSTART, RLENGTH)
        if (g in keep) { $4 = $4 - S[$1] + 1; $5 = $5 - S[$1] + 1; print } }
' > "$TMP/mini_annotation.gtf"

# ── 3. mini transcripts FASTA (transcripts of those genes) ──
echo "[subset] mini transcripts"
grep -oP 'transcript_id "\K[^"]+' "$TMP/mini_annotation.gtf" | sort -u > "$TMP/tx_ids.txt"
zcat "$TX_GZ" | awk 'NR == FNR { k[$1] = 1; next }
    /^>/ { split(substr($0, 2), h, "|"); p = (h[1] in k) }
    p' "$TMP/tx_ids.txt" - > "$TMP/mini_transcripts.fa"

gzip -n -9 -c "$TMP/mini_genome.fa"       > "$OUT_DIR/ref/mini_genome.fa.gz"
gzip -n -9 -c "$TMP/mini_annotation.gtf"  > "$OUT_DIR/ref/mini_annotation.gtf.gz"
gzip -n -9 -c "$TMP/mini_transcripts.fa"  > "$OUT_DIR/ref/mini_transcripts.fa.gz"

# ── 4. reads ─────────────────────────────────────────────────
tail -n +2 "$SAMPLES" | while IFS=$'\t' read -r sample gsm group genotype age_group counts_col run1 run2; do
    bam="$BAM_DIR/${sample}.Aligned.sortedByCoord.out.bam"
    [[ -s "$bam" ]] || { echo "ERROR: missing $bam" >&2; exit 1; }
    : > "$TMP/names.txt"
    while read -r reg target; do
        n=$("$SAMTOOLS" view -c -f 2 -F 0x904 "$bam" "$reg")
        frac=$(awk -v t="$target" -v n="$n" 'BEGIN { f = (n > 0) ? 2 * t / n : 1; if (f > 1) f = 1; printf "%.6f", f }')
        "$SAMTOOLS" view -f 2 -F 0x904 --subsample "$frac" --subsample-seed "$SEED" "$bam" "$reg" \
            | cut -f1 >> "$TMP/names.txt"
    done <<< "$TARGETS"
    sort -u "$TMP/names.txt" > "$TMP/names.u.txt"

    # shellcheck disable=SC2086
    "$SAMTOOLS" view -u -f 2 -F 0x904 -N "$TMP/names.u.txt" "$bam" $FULL_REGIONS \
        | "$SAMTOOLS" collate -u -O - "$TMP/collate" \
        | "$SAMTOOLS" fastq -n -1 "$TMP/r1.fq" -2 "$TMP/r2.fq" -0 /dev/null -s /dev/null - 2>/dev/null

    for run in "$run1" "$run2"; do
        for m in 1 2; do
            awk -v r="@$run." 'NR % 4 == 1 { p = (index($0, r) == 1) } p' "$TMP/r$m.fq" \
                | gzip -n -6 > "$OUT_DIR/raw_data/fastq/${run}_${m}.fastq.gz"
        done
    done
    echo "[subset] $sample: $(( $(wc -l < "$TMP/r1.fq") / 4 )) pairs  ($run1 + $run2)"
done

echo "[subset] done -> $OUT_DIR  ($(du -sh "$OUT_DIR" | cut -f1))"
