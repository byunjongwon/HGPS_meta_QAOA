#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────
# 전체 파이프라인 순차 실행
#
#   각 단계는 완료된 산출물을 건너뛰므로 재실행에 안전합니다.
#   실패하면 그 자리에서 멈추고 FAILED 마커를 남깁니다.
#
# 사용법:
#   bash scripts/run_pipeline.sh            전체 실행
#   START=09 bash scripts/run_pipeline.sh   09단계부터 재개
#   ONLY=13   bash scripts/run_pipeline.sh  13단계만 실행
#
# 상태 확인:
#   cat $PROJ/logs/PIPELINE_STATUS            # [deposit] was ~/progerin_rnaseq
#   tail -f $PROJ/logs/<번호>_stdout.log
# ─────────────────────────────────────────────────────────────
# [deposit] scripts = this runner's own dir; PROJ default = its parent (overridable, e.g. a test dir)
SCRIPTS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PROJ="${PROJ:-$(dirname "$SCRIPTS")}"
LOGS="$PROJ/logs"
STATUS="$LOGS/PIPELINE_STATUS"
mkdir -p "$LOGS"

# 00_config.sh의 변수를 R 스크립트에서도 쓰도록 export
set -a; source "$SCRIPTS/00_config.sh"; set +a

START="${START:-03}"
ONLY="${ONLY:-}"

stamp() { date '+%F %T'; }
set_status() { echo "$(stamp)  $*" | tee -a "$STATUS"; }

run_step() {
    local num="$1" name="$2" cmd="$3"

    if [[ -n "$ONLY" ]]; then
        [[ "$num" == "$ONLY" ]] || return 0
    else
        [[ "$num" < "$START" ]] && { echo "$(stamp)  ⏭  STEP $num 건너뜀 (START=$START)"; return 0; }
    fi

    set_status "▶ START  STEP $num — $name"
    if eval "$cmd" >> "$LOGS/${num}_stdout.log" 2>&1; then
        set_status "✅ DONE   STEP $num — $name"
    else
        local rc=$?
        set_status "❌ FAILED STEP $num — $name (exit $rc)"
        set_status "   로그: tail -50 $LOGS/${num}_stdout.log"
        echo "FAILED at STEP $num" > "$LOGS/PIPELINE_FAILED"
        exit "$rc"
    fi
}

rm -f "$LOGS/PIPELINE_FAILED" "$LOGS/PIPELINE_DONE"
set_status "════════ 파이프라인 시작 (START=$START ONLY=${ONLY:-none}) ════════"

run_step 03 "FastQC (raw) + MultiQC"                    "bash '$SCRIPTS/03_fastqc_raw.sh'"
run_step 04 "fastp 트리밍 + 트림 후 QC"                  "bash '$SCRIPTS/04_trim.sh'"
run_step 05 "STAR 정렬"                                  "bash '$SCRIPTS/05_star_align.sh'"
run_step 06 "정렬 후 QC (RSeQC/Picard/dupRadar)"         "bash '$SCRIPTS/06_post_align_qc.sh'"
run_step 07 "featureCounts counts 행렬"                  "bash '$SCRIPTS/07_featurecounts.sh'"
run_step 08 "Salmon transcript 정량"                     "bash '$SCRIPTS/08_salmon.sh'"
run_step 09 "DESeq2 DEG + 그림"                          "'$RSCRIPT' '$SCRIPTS/09_deseq2.R'"
run_step 10 "기능 분석 (GO/KEGG/Reactome/GSEA)"          "'$RSCRIPT' '$SCRIPTS/10_enrichment.R'"
run_step 11 "시그니처 스코어링"                          "'$RSCRIPT' '$SCRIPTS/11_signatures.R'"
run_step 12 "경로·전사인자 활성"                         "'$RSCRIPT' '$SCRIPTS/12_pathway_tf.R'"
run_step 13 "rMATS 차등 스플라이싱"                      "bash '$SCRIPTS/13_rmats.sh'"
# [deposit] steps 14 (rMATS figures), 15 (sashimi), 18 (report) removed: figure/report-only code not in the deposit
run_step 16 "IsoformSwitchAnalyzeR DTU"                  "'$RSCRIPT' '$SCRIPTS/16_isoformswitch.R'"
run_step 17 "WGCNA 공발현 모듈"                          "'$RSCRIPT' '$SCRIPTS/17_wgcna.R'"

set_status "════════ 전체 완료 ════════"
echo "ALL DONE" > "$LOGS/PIPELINE_DONE"
