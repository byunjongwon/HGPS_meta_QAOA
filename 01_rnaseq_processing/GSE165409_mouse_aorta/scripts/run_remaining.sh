#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────
# STEP 4(마무리) → 5 → 6 → 7 → 8 자동 연속 실행
#
# 에러가 나면 그 자리에서 멈추고 FAILED 마커를 남깁니다.
# 각 스크립트는 완료된 산출물을 건너뛰므로 재실행에 안전합니다.
#
# 상태 확인:
#   cat $PROJ/logs/PIPELINE_STATUS
#   tail -f $PROJ/logs/run_remaining.log   # [deposit] was ~/lmna_rnaseq
# ─────────────────────────────────────────────────────────────
# [deposit] SCRIPTS = directory of this runner; default PROJ = dataset dir (scripts/..), overridable by PROJ
SCRIPTS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PROJ="${PROJ:-$(cd "$SCRIPTS/.." && pwd)}"
export CONDA_ENVS="${CONDA_ENVS:-$HOME/miniforge3/envs}"
LOGS="$PROJ/logs"
STATUS="$LOGS/PIPELINE_STATUS"
mkdir -p "$LOGS"

# Qualimap은 건너뜀 (RSeQC 결과와 중복, 샘플당 11분)
export SKIP_QUALIMAP=1

stamp() { date '+%F %T'; }

set_status() { echo "$(stamp)  $*" | tee -a "$STATUS"; }

run_step() {
    local num="$1" name="$2" cmd="$3"
    set_status "▶ START  STEP $num — $name"
    if eval "$cmd" >> "$LOGS/${num}_stdout.log" 2>&1; then
        set_status "✅ DONE   STEP $num — $name"
    else
        local rc=$?
        set_status "❌ FAILED STEP $num — $name (exit $rc)"
        set_status "   로그 확인: $LOGS/${num}_stdout.log"
        echo "FAILED at STEP $num" > "$LOGS/PIPELINE_FAILED"
        exit "$rc"
    fi
}

rm -f "$LOGS/PIPELINE_FAILED"
set_status "════════ 파이프라인 자동 실행 시작 ════════"

run_step 04 "정렬 후 QC 마무리 (strandedness + MultiQC)" \
    "bash '$SCRIPTS/04_post_align_qc.sh'"

run_step 05 "featureCounts — counts 행렬 + 저자 데이터 검증" \
    "bash '$SCRIPTS/05_featurecounts.sh'"

run_step 06 "Salmon — transcript 정량 (DTU용)" \
    "bash '$SCRIPTS/06_salmon.sh'"

run_step 07 "rMATS — 차등 스플라이싱 3개 비교" \
    "bash '$SCRIPTS/07_rmats.sh'"

run_step 08 "DESeq2 — DEG 분석 + 그림 13종" \
    "'$CONDA_ENVS/lmna_r/bin/Rscript' '$SCRIPTS/08_deseq2.R'"

set_status "════════ 전체 완료 ════════"
echo "ALL DONE" > "$LOGS/PIPELINE_DONE"
