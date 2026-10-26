#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────
# STEP 12 → 10 → 11 → 13 재개 (08, 09는 완료됨)
# 에러 시 중단하고 마커를 남깁니다.
# ─────────────────────────────────────────────────────────────
# [deposit] SCRIPTS = directory of this runner; default PROJ = dataset dir (scripts/..), overridable by PROJ
SCRIPTS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PROJ="${PROJ:-$(cd "$SCRIPTS/.." && pwd)}"
export CONDA_ENVS="${CONDA_ENVS:-$HOME/miniforge3/envs}"
LOGS="$PROJ/logs"
STATUS="$LOGS/PIPELINE2_STATUS"
RSCRIPT="$CONDA_ENVS/lmna_r/bin/Rscript"   # [deposit] was $HOME/miniforge3/envs/...

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
set_status "════════ 재개 (STEP 12부터) ════════"

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
