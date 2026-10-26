#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────
# STEP 13 — rMATS 차등 스플라이싱 (Progerin vs Control)
#
# progerin은 그 자체가 스플라이싱 산물입니다. LMNA exon 11의 잠재
# 5' splice donor가 활성화되어 150 nt가 잘린 전사체가 만들어집니다.
# 다만 이 데이터셋은 렌티바이러스로 progerin cDNA를 넣은 것이라
# 내인성 LMNA 스플라이싱 변화 자체가 주 관심은 아니고,
# "progerin 발현이 세포 전반의 스플라이싱을 어떻게 바꾸는가"를 봅니다.
#
# 산출: SE/MXE/A5SS/A3SS/RI .MATS.JC(EC).txt
#       → ΔPSI volcano, event type bar, PSI heatmap, sashimi 재료
#
# ⚠️ --variable-read-length 필수 (STEP 04에서 트리밍해 길이 가변)
# ⚠️ --novelSS : annotation에 없는 cryptic splice site 검출용
# ─────────────────────────────────────────────────────────────
source "$(dirname "$0")/00_config.sh"

BAM_DIR="$RESULT_DIR/05_star"
OUT="$RESULT_DIR/13_rmats"
mkdir -p "$OUT"

# ── libType 결정 ──────────────────────────────────────────────
STRAND_FILE="$RESULT_DIR/06_qc_align/strandedness_summary.tsv"
if [[ -z "${RMATS_LIBTYPE:-}" && -s "$STRAND_FILE" ]]; then
    read -r AVG_FWD AVG_REV < <(awk -F'\t' 'NR>1 && $3!="NA" {f+=$3; r+=$4; n++} END{printf "%.3f %.3f\n", f/n, r/n}' "$STRAND_FILE")
    if   (( $(echo "$AVG_REV > 0.75" | bc -l) )); then RMATS_LIBTYPE="fr-firststrand"
    elif (( $(echo "$AVG_FWD > 0.75" | bc -l) )); then RMATS_LIBTYPE="fr-secondstrand"
    else RMATS_LIBTYPE="fr-unstranded"; fi
    log "rMATS libType 자동 판정: $RMATS_LIBTYPE (fwd=$AVG_FWD rev=$AVG_REV)"
fi
RMATS_LIBTYPE="${RMATS_LIBTYPE:-fr-unstranded}"

# ── 그룹별 BAM 목록 ───────────────────────────────────────────
make_list() {
    local grp="$1" outfile="$2"
    sample_rows | awk -F'\t' -v g="$grp" '$4==g {print $1}' | \
        while read -r s; do echo -n "$BAM_DIR/${s}.Aligned.sortedByCoord.out.bam,"; done | \
        sed 's/,$//' > "$outfile"
    log "  $grp: $(tr ',' '\n' < "$outfile" | wc -l) 샘플"
}
log "그룹별 BAM 목록 생성"
make_list Progerin "$OUT/bam_Progerin.txt"
make_list Control  "$OUT/bam_Control.txt"

READ_LEN="${READ_LEN:-151}"

NAME="Progerin_vs_Control"
ODIR="$OUT/$NAME"
if [[ -s "$ODIR/SE.MATS.JC.txt" ]]; then
    log "  $NAME — 이미 완료, 건너뜀"
else
    # ⚠️ tmp에 이전(중단된) 실행의 .rmats가 남아 있으면 BAM 개수와 어긋나
    #    통계 단계에서 죽습니다. 항상 비우고 시작합니다.
    rm -rf "$OUT/tmp_$NAME"
    mkdir -p "$ODIR" "$OUT/tmp_$NAME"

    log "rMATS 실행 (libType=$RMATS_LIBTYPE, readLength=$READ_LEN, novelSS)"
    # 전체 출력을 로그 파일에 남깁니다. tail로 자르면 에러 메시지가 사라져
    # "왜 죽었는지 모르는" 상황이 됩니다.
    if "$RMATS" \
        --b1 "$OUT/bam_Progerin.txt" --b2 "$OUT/bam_Control.txt" \
        --gtf "$GTF" \
        -t paired \
        --readLength "$READ_LEN" \
        --variable-read-length \
        --libType "$RMATS_LIBTYPE" \
        --nthread "$THREADS" \
        --od "$ODIR" \
        --tmp "$OUT/tmp_$NAME" \
        --novelSS \
        > "$LOG_DIR/13_rmats_full.log" 2>&1; then
        log "  $NAME 완료"
    else
        rc=$?
        log "  ❌ rMATS 실패 (exit $rc) — 마지막 30줄:"
        tail -30 "$LOG_DIR/13_rmats_full.log" | tee -a "$LOG_DIR/13_rmats.log"
        exit "$rc"
    fi
fi

# ── 요약 ──────────────────────────────────────────────────────
log "────────── 유의 event (FDR<0.05, |ΔPSI|>0.1) ──────────"
{
    echo -e "comparison\tevent\ttotal\tsignificant"
    for e in SE MXE A5SS A3SS RI; do
        f="$ODIR/${e}.MATS.JC.txt"
        [[ -s "$f" ]] || continue
        tot=$(( $(wc -l < "$f") - 1 ))
        sig=$(awk -F'\t' 'NR>1 && $(NF-2)<0.05 && ($(NF)>0.1 || $(NF)< -0.1)' "$f" | wc -l)
        echo -e "$NAME\t$e\t$tot\t$sig"
    done
} | tee "$OUT/event_summary.tsv" | column -t

log "완료 → $OUT"
log "※ SE=exon skipping, MXE=mutually exclusive, A5SS/A3SS=대체 splice site, RI=intron retention"
