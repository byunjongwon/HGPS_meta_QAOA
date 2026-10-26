#!/usr/bin/env bash
# [deposit] created for the code deposit: builds the tiny test subset in test/ (a PROJ-like dir) that lets
# the raw-data steps 02-15 run end-to-end in minutes instead of on ~27 GB FASTQ + ~43 GB reference.
#
# What it makes (from the ORIGINAL full-size run; read-only access):
#   ref/GRCh38.primary_assembly.genome.fa.gz              1 Mb window chr1:155,600,001-156,600,000 (GRCh38),
#                                                         contig still named "chr1", coordinates shifted by
#                                                         -155,600,000 (window start becomes position 1)
#   ref/gencode.v47.primary_assembly.annotation.gtf.gz    GENCODE v47 genes lying entirely inside the window,
#                                                         coordinates shifted the same way
#   ref/gencode.v47.transcripts.fa.gz                     GENCODE v47 transcripts of those genes
#   raw_data/fastq/<SRR>_{1,2}.fastq.gz                   per sample: properly paired, primary alignments in the
#                                                         window, subsampled (seeded) to ~N_LMNA pairs touching
#                                                         LMNA (chr1:156,082,572-156,140,081, incl. exon 11) +
#                                                         ~N_REST pairs elsewhere in the window
#   meta/samples.tsv                                      copy of the real sample sheet (8 samples, 4 vs 4)
#   meta_download/ena_filereport.tsv                      ENA report of one genuinely tiny public run
#                                                         (SRR19598819, 2 x ~1 kB) to exercise 01_download_fastq.sh
#
# NOTE: reads come from the STAR BAMs (i.e. already fastp-trimmed reads, converted back to FASTQ with
# samtools fastq, which restores the original read orientation). They are re-trimmed by 04_trim.sh in the test.
#
# Usage:  bash make_test_subset.sh                       (writes next to this script)
#         OUT=/some/dir SRC_PROJ=/path/to/full/run bash make_test_subset.sh
set -euo pipefail

SRC_PROJ="${SRC_PROJ:?set SRC_PROJ to the full-size original run (STAR BAMs + GRCh38 reference)}"   # [deposit] no hard-coded user path
OUT="${OUT:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
CONDA_ENVS="${CONDA_ENVS:-$HOME/miniforge3/envs}"
SAMTOOLS="${SAMTOOLS:-$CONDA_ENVS/lmna_rmats/bin/samtools}"
THREADS="${THREADS:-4}"

CHR=chr1
WSTART=155600000            # 0-based offset of the window (window = WSTART+1 .. WEND, 1-based inclusive)
WEND=156600000
LMNA_REG="chr1:156082572-156140081"   # LMNA gene (GENCODE v47)
N_LMNA="${N_LMNA:-3000}"    # target read pairs touching LMNA per sample
N_REST="${N_REST:-3000}"    # target read pairs elsewhere in the window per sample
SEED=42

GENOME="$SRC_PROJ/ref/GRCh38.primary_assembly.genome.fa"
GTFGZ="$SRC_PROJ/ref/gencode.v47.primary_assembly.annotation.gtf.gz"
TXGZ="$SRC_PROJ/ref/gencode.v47.transcripts.fa.gz"
BAMDIR="$SRC_PROJ/results/05_star"
SAMPLES="$SRC_PROJ/meta/samples.tsv"
for f in "$GENOME" "$GTFGZ" "$TXGZ" "$SAMPLES"; do [[ -s "$f" ]] || { echo "missing $f" >&2; exit 1; }; done

W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
mkdir -p "$OUT/ref" "$OUT/raw_data/fastq" "$OUT/meta" "$OUT/meta_download"
cp "$SAMPLES" "$OUT/meta/samples.tsv"

# ── mini genome ───────────────────────────────────────────────
echo "[genome] extracting $CHR:$((WSTART+1))-$WEND"
awk -v c=">$CHR" '/^>/{ if (p) exit; p = ($1 == c) } p' "$GENOME" > "$W/chr.fa"   # chr1 is the first record
"$SAMTOOLS" faidx "$W/chr.fa"
"$SAMTOOLS" faidx "$W/chr.fa" "$CHR:$((WSTART+1))-$WEND" | sed "1s/.*/>$CHR/" | gzip -n -9 \
    > "$OUT/ref/GRCh38.primary_assembly.genome.fa.gz"

# ── mini GTF (genes fully inside the window; shifted coordinates) ──
echo "[gtf] genes fully inside the window"
zcat "$GTFGZ" | awk -F'\t' -v OFS='\t' -v c="$CHR" -v s="$WSTART" -v e="$WEND" '
    /^#/ { print; next }
    $1 != c { next }
    $3 == "gene" { match($9, /gene_id "[^"]+"/); g = substr($9, RSTART, RLENGTH); keep[g] = ($4 > s && $5 <= e) }
    { match($9, /gene_id "[^"]+"/); g = substr($9, RSTART, RLENGTH)
      if (keep[g]) { $4 -= s; $5 -= s; print } }' | gzip -n -9 \
    > "$OUT/ref/gencode.v47.primary_assembly.annotation.gtf.gz"
zcat "$OUT/ref/gencode.v47.primary_assembly.annotation.gtf.gz" | awk -F'\t' '$3=="transcript"{
    match($9, /transcript_id "[^"]+"/); print substr($9, RSTART+15, RLENGTH-16) }' | sort -u > "$W/tx.ids"
echo "        $(zcat "$OUT/ref/gencode.v47.primary_assembly.annotation.gtf.gz" | awk -F'\t' '$3=="gene"' | wc -l) genes, $(wc -l < "$W/tx.ids") transcripts"

# ── mini transcript FASTA ─────────────────────────────────────
zcat "$TXGZ" | awk 'NR==FNR { k[$1]; next }
    /^>/ { split(substr($0, 2), h, "|"); p = (h[1] in k) } p' "$W/tx.ids" - | gzip -n -9 \
    > "$OUT/ref/gencode.v47.transcripts.fa.gz"

# ── reads ────────────────────────────────────────────────────
tail -n +2 "$SAMPLES" | while IFS=$'\t' read -r sample gsm srr group replicate; do
    bam="$BAMDIR/${sample}.Aligned.sortedByCoord.out.bam"
    [[ -s "$bam" ]] || { echo "missing $bam" >&2; exit 1; }
    # primary, properly paired alignments in the window
    "$SAMTOOLS" view -@ "$THREADS" -b -F 0x904 -f 0x2 -o "$W/win.bam" "$bam" "$CHR:$((WSTART+1))-$WEND"
    "$SAMTOOLS" index "$W/win.bam"
    "$SAMTOOLS" view "$W/win.bam" "$LMNA_REG" | cut -f1 | sort -u > "$W/lmna.names"
    "$SAMTOOLS" view "$W/win.bam" | cut -f1 | sort -u > "$W/all.names"
    comm -23 "$W/all.names" "$W/lmna.names" > "$W/rest.names"
    nl=$(wc -l < "$W/lmna.names"); nr=$(wc -l < "$W/rest.names")
    fl=$(awk -v n="$N_LMNA" -v t="$nl" 'BEGIN{ f = (t > 0) ? n / t : 1; if (f > 1) f = 1; printf "%.6f", f }')
    fr=$(awk -v n="$N_REST" -v t="$nr" 'BEGIN{ f = (t > 0) ? n / t : 1; if (f > 1) f = 1; printf "%.6f", f }')
    # name-hash based subsampling keeps both mates of a pair together and is deterministic for a seed
    "$SAMTOOLS" view -b -N "$W/lmna.names" --subsample "$fl" --subsample-seed "$SEED" -o "$W/l.bam" "$W/win.bam"
    "$SAMTOOLS" view -b -N "$W/rest.names" --subsample "$fr" --subsample-seed "$SEED" -o "$W/r.bam" "$W/win.bam"
    "$SAMTOOLS" cat -o "$W/sel.bam" "$W/l.bam" "$W/r.bam"
    "$SAMTOOLS" sort -n -@ "$THREADS" -o "$W/sel.n.bam" "$W/sel.bam"
    # -n: keep identical read names for both mates; mates outside the window (singletons) are dropped
    "$SAMTOOLS" fastq -n -1 "$W/r1.fq" -2 "$W/r2.fq" -0 /dev/null -s /dev/null "$W/sel.n.bam" 2>/dev/null
    gzip -n -c "$W/r1.fq" > "$OUT/raw_data/fastq/${srr}_1.fastq.gz"
    gzip -n -c "$W/r2.fq" > "$OUT/raw_data/fastq/${srr}_2.fastq.gz"
    echo "[reads] $sample ($srr): LMNA pairs $nl -> frac $fl, other pairs $nr -> frac $fr; written $(( $(wc -l < "$W/r1.fq") / 4 )) pairs"
    rm -f "$W"/*.bam "$W"/*.bai "$W"/*.fq "$W"/*.names
done

# ── download test: one genuinely tiny public ENA run (same 13 columns as meta/ena_filereport.tsv) ──
cat > "$OUT/meta_download/ena_filereport.tsv" <<'EOF'
run_accession	experiment_accession	sample_accession	sample_title	sample_alias	library_layout	library_strategy	instrument_model	read_count	base_count	fastq_ftp	fastq_md5	fastq_bytes
SRR19598819	SRX15650195	SAMN28938782	sc-Patch-seq,Red only,Regenerating CSMN98	GSM6226012	PAIRED	RNA-Seq	Illumina NovaSeq 6000	15	4500	ftp.sra.ebi.ac.uk/vol1/fastq/SRR195/019/SRR19598819/SRR19598819_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/SRR195/019/SRR19598819/SRR19598819_2.fastq.gz	8e745ba9f4b5a2666cbb195ffea69c5d;498fcc8c03c7c0e89d1f002784a7c499	983;1075
EOF

echo "done -> $OUT ($(du -sh "$OUT" | cut -f1))"
