#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────
# STEP 7. rMATS — 차등 스플라이싱 (이 데이터셋의 핵심 분석)
#
# LMNA G609G는 그 자체가 스플라이싱 변이입니다. exon 11의 잠재
# 5' splice donor를 활성화시켜 150 nt가 잘린 progerin 전사체를 만듭니다.
# counts 행렬에서는 Lmna 유전자 발현 하나로 뭉뚱그려져 progerin이
# 전혀 보이지 않으므로, 이 단계가 main figure의 근거가 됩니다.
#
# 3개 비교:
#   LMNA_vs_WTyng   progeria 효과
#   WTold_vs_WTyng  정상 노화 효과 (intron retention 증가 예상)
#   LMNA_vs_WTold   progeria vs 노화
#
# 산출 (비교마다):
#   SE/MXE/A5SS/A3SS/RI .MATS.JC.txt, .MATS.JCEC.txt
#   → ΔPSI volcano, AS event type bar, PSI heatmap, sashimi plot 재료
#
# ⚠️ --variable-read-length 필수: 02에서 트리밍해 길이가 가변이 됨
#
# 예상 시간: 비교당 30~60분
# ─────────────────────────────────────────────────────────────
source "$(dirname "$0")/00_config.sh"

BAM_DIR="$RESULT_DIR/03_star"
OUT="$RESULT_DIR/07_rmats"
mkdir -p "$OUT"

# ── libType 결정 (04 판정 결과 반영) ──────────────────────────
STRAND_FILE="$RESULT_DIR/04_qc_align/strandedness_summary.tsv"
if [[ -z "${RMATS_LIBTYPE:-}" && -s "$STRAND_FILE" ]]; then
    # printf 끝의 개행 필수 (없으면 read가 EOF로 exit 1 → set -e가 중단시킴)
    read -r AVG_FWD AVG_REV < <(awk -F'\t' 'NR>1 && $3!="NA" {f+=$3; r+=$4; n++} END{printf "%.3f %.3f\n", f/n, r/n}' "$STRAND_FILE")
    if   (( $(echo "$AVG_REV > 0.75" | bc -l) )); then RMATS_LIBTYPE="fr-firststrand"
    elif (( $(echo "$AVG_FWD > 0.75" | bc -l) )); then RMATS_LIBTYPE="fr-secondstrand"
    else RMATS_LIBTYPE="fr-unstranded"; fi
    log "rMATS libType 자동 판정: $RMATS_LIBTYPE"
fi
RMATS_LIBTYPE="${RMATS_LIBTYPE:-fr-unstranded}"

# ── 그룹별 BAM 목록 생성 ──────────────────────────────────────
make_list() {
    local grp="$1" outfile="$2"
    sample_rows | awk -F'\t' -v g="$grp" '$3==g {print $1}' | \
        while read -r s; do echo -n "$BAM_DIR/${s}.Aligned.sortedByCoord.out.bam,"; done | \
        sed 's/,$//' > "$outfile"
    log "  $grp: $(tr ',' '\n' < "$outfile" | wc -l) 샘플"
}

log "그룹별 BAM 목록 생성"
make_list LMNA  "$OUT/bam_LMNA.txt"
make_list WTyng "$OUT/bam_WTyng.txt"
make_list WTold "$OUT/bam_WTold.txt"

# ── read 길이 (트리밍 후 최대값) ──────────────────────────────
READ_LEN="${READ_LEN:-150}"

run_comparison() {
    local name="$1" b1="$2" b2="$3"
    local odir="$OUT/$name"
    if [[ -s "$odir/SE.MATS.JC.txt" ]]; then
        log "  $name — 이미 완료, 건너뜀"
        return
    fi
    mkdir -p "$odir" "$OUT/tmp_$name"
    log "  $name 실행 중..."
    "$RMATS" \
        --b1 "$b1" --b2 "$b2" \
        --gtf "$GTF" \
        -t paired \
        --readLength "$READ_LEN" \
        --variable-read-length \
        --libType "$RMATS_LIBTYPE" \
        --nthread "$THREADS" \
        --od "$odir" \
        --tmp "$OUT/tmp_$name" \
        --novelSS \
        2>&1 | tail -20 | tee -a "$LOG_DIR/07_rmats.log"
    log "  $name 완료"
}

log "rMATS 시작 (libType=$RMATS_LIBTYPE, readLength=$READ_LEN, variable-read-length)"
run_comparison "LMNA_vs_WTyng"  "$OUT/bam_LMNA.txt"  "$OUT/bam_WTyng.txt"
run_comparison "WTold_vs_WTyng" "$OUT/bam_WTold.txt" "$OUT/bam_WTyng.txt"
run_comparison "LMNA_vs_WTold"  "$OUT/bam_LMNA.txt"  "$OUT/bam_WTold.txt"

# ── 요약 ──────────────────────────────────────────────────────
log "────────── 유의 event 개수 (FDR<0.05, |ΔPSI|>0.1) ──────────"
{
    echo -e "comparison\tevent\ttotal\tsignificant"
    for c in LMNA_vs_WTyng WTold_vs_WTyng LMNA_vs_WTold; do
        for e in SE MXE A5SS A3SS RI; do
            f="$OUT/$c/${e}.MATS.JC.txt"
            [[ -s "$f" ]] || continue
            tot=$(( $(wc -l < "$f") - 1 ))
            sig=$(awk -F'\t' 'NR>1 && $(NF-2)<0.05 && ($(NF)>0.1 || $(NF)< -0.1)' "$f" | wc -l)
            echo -e "$c\t$e\t$tot\t$sig"
        done
    done
} | tee "$OUT/event_summary.tsv" | column -t

log "완료 → $OUT"
log "※ --novelSS 를 켜 두었습니다. progerin의 cryptic splice site는 annotation에"
log "  없으므로 이 옵션 없이는 검출되지 않습니다."
