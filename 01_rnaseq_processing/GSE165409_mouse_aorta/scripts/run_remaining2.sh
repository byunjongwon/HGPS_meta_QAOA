#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────
# STEP 9 → 13 자동 연속 실행 (rMATS 그림 / sashimi / DTU / 기능분석 / dupRadar)
# 에러 시 그 자리에서 중단하고 FAILED 마커를 남깁니다.
# ─────────────────────────────────────────────────────────────
# [deposit] SCRIPTS = directory of this runner; default PROJ = dataset dir (scripts/..), overridable by PROJ
SCRIPTS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PROJ="${PROJ:-$(cd "$SCRIPTS/.." && pwd)}"
export CONDA_ENVS="${CONDA_ENVS:-$HOME/miniforge3/envs}"
LOGS="$PROJ/logs"
STATUS="$LOGS/PIPELINE2_STATUS"
RSCRIPT="$CONDA_ENVS/lmna_r/bin/Rscript"   # [deposit] was $HOME/miniforge3/envs/...
mkdir -p "$LOGS"

set_status() { echo "$(date '+%F %T')  $*" | tee -a "$STATUS"; }

run_step() {
    local num="$1" name="$2" cmd="$3"
    set_status "▶ START  STEP $num — $name"
    if eval "$cmd" >> "$LOGS/${num}_stdout.log" 2>&1; then
        set_status "✅ DONE   STEP $num — $name"
    else
        local rc=$?
        set_status "❌ FAILED STEP $num — $name (exit $rc)"
        echo "FAILED at STEP $num" > "$LOGS/PIPELINE2_FAILED"
        exit "$rc"
    fi
}

rm -f "$LOGS/PIPELINE2_FAILED" "$LOGS/PIPELINE2_DONE"
set_status "════════ 2차 파이프라인 시작 ════════"

# 08은 이미 실행됐지만 그림에 Δ·한글이 깨져 들어갔으므로 cairo 장치로 다시 그립니다.
run_step 08 "DESeq2 그림 재생성 (폰트 수정)" \
    "PROJ='$PROJ' '$RSCRIPT' '$SCRIPTS/08_deseq2.R'"

run_step 09 "rMATS 시각화 (ΔPSI volcano, PSI heatmap, Lmna exon11)" \
    "PROJ='$PROJ' '$RSCRIPT' '$SCRIPTS/09_rmats_figures.R'"

run_step 12 "기능 분석 (GO + GSEA)" \
    "PROJ='$PROJ' '$RSCRIPT' '$SCRIPTS/12_enrichment.R'"

# [deposit] STEP 10 (10_sashimi.sh, figure-only) removed from the code deposit
# run_step 10 "Sashimi plot (Lmna exon11 + 상위 event)" "bash '$SCRIPTS/10_sashimi.sh'"

run_step 11 "IsoformSwitchAnalyzeR — DTU" \
    "PROJ='$PROJ' '$RSCRIPT' '$SCRIPTS/11_isoformswitch.R'"

run_step 13 "dupRadar — PCR duplicate 판별" \
    "bash '$SCRIPTS/13_dupradar.sh'"

set_status "════════ 전체 완료 ════════"
echo "ALL DONE" > "$LOGS/PIPELINE2_DONE"
