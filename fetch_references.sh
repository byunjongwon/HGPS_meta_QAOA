#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# Download the GENCODE reference files that are not shipped in this repository
# (public, unmodified files; ~245 MB). They are needed only by the isoform-switch
# scripts of stage 01 (GSE165409 11_isoformswitch.R, GSE198983 16_isoformswitch.R)
# and by the full raw-read reprocessing. No other result depends on them.
#
#   bash fetch_references.sh
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
M="$ROOT/01_rnaseq_processing/GSE165409_mouse_aorta/ref"
H="$ROOT/01_rnaseq_processing/GSE198983_human_AoSMC/ref"
GM=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25
GH=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47

# md5 = GENCODE MD5SUMS (identical content to the files used in the analysis)
fetch() {   # url dest md5
    local url=$1 dest=$2 md5=$3
    mkdir -p "$(dirname "$dest")"
    if [[ -s "$dest" ]] && echo "$md5  $dest" | md5sum -c --quiet 2>/dev/null; then
        echo "ok (cached) $dest"; return
    fi
    echo "download $url"
    curl -fL --retry 3 -o "$dest" "$url"
    echo "$md5  $dest" | md5sum -c
}

fetch $GM/gencode.vM25.primary_assembly.annotation.gtf.gz $M/gencode.vM25.primary_assembly.annotation.gtf.gz c5125258a0a2c5250ddb4c192abbf4e8
fetch $GM/gencode.vM25.transcripts.fa.gz                  $M/gencode.vM25.transcripts.fa.gz                  a821c0dde39c48b9d2c4b48d36b0180c
fetch $GH/gencode.v47.primary_assembly.annotation.gtf.gz  $H/gencode.v47.primary_assembly.annotation.gtf.gz  81cf2337bb31aeaddba7059a3b0f14f5
fetch $GH/gencode.v47.transcripts.fa.gz                   $H/gencode.v47.transcripts.fa.gz                   9c3ce49c6b6f13f3e4f4deb53ca80bf9
echo "all GENCODE references present"
