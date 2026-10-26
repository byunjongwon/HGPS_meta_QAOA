#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────
# STEP 01 — ENA에서 FASTQ 다운로드 (병렬) + md5 검증
#   ENA는 연결당 대역폭을 제한하므로 PAR개 연결을 동시에 사용합니다.
#   이미 받은 파일은 md5 확인 후 건너뜁니다 (재실행 안전).
#   md5 불일치 시 해당 파일만 삭제 후 재시도 (최대 3회).
# ─────────────────────────────────────────────────────────────
source "$(dirname "$0")/00_config.sh"

PAR="${PAR:-8}"                     # 동시 연결 수
REPORT="$META_DIR/ena_filereport.tsv"
require_file "$REPORT"
mkdir -p "$FASTQ_DIR"

URLS="$META_DIR/fastq_urls.tsv"     # url<TAB>md5<TAB>filename
awk -F'\t' 'NR>1{
    n=split($11,u,";"); split($12,m,";");
    for(i=1;i<=n;i++){
        split(u[i],p,"/"); f=p[length(p)];
        print "https://"u[i]"\t"m[i]"\t"f
    }
}' "$REPORT" | sort -k3,3 > "$URLS"

log "다운로드 대상 $(wc -l < "$URLS") 개 파일 / 동시 연결 $PAR"

# 한 파일을 받는 워커 (xargs가 호출)
grab() {
    local url="$1" md5="$2" fname="$3"
    local dest="$FASTQ_DIR/$fname"

    if [[ -s "$dest" ]] && [[ "$(md5sum "$dest" | cut -d' ' -f1)" == "$md5" ]]; then
        echo "[$(date '+%T')] SKIP  $fname (md5 OK)"
        return 0
    fi

    for attempt in 1 2 3; do
        echo "[$(date '+%T')] GET   $fname (시도 $attempt)"
        if curl -fsS -C - --retry 5 --retry-delay 10 --connect-timeout 30 \
                -o "$dest" "$url"; then
            if [[ "$(md5sum "$dest" | cut -d' ' -f1)" == "$md5" ]]; then
                echo "[$(date '+%T')] OK    $fname ($(du -h "$dest" | cut -f1))"
                return 0
            fi
            echo "[$(date '+%T')] WARN  $fname md5 불일치 → 삭제 후 재시도"
            rm -f "$dest"
        else
            echo "[$(date '+%T')] WARN  $fname 전송 실패 → 재시도"
        fi
        sleep 5
    done
    echo "[$(date '+%T')] ERROR $fname 3회 실패"
    return 1
}
export -f grab
export FASTQ_DIR

# xargs로 병렬 실행. 하나라도 실패하면 xargs가 비0 종료.
if xargs -a "$URLS" -d'\n' -P "$PAR" -I{} bash -c 'IFS=$'"'"'\t'"'"' read -r u m f <<< "{}"; grab "$u" "$m" "$f"'; then
    STATUS=0
else
    STATUS=1
fi

log "────────────────────────────────────────"
log "총 용량: $(du -sh "$FASTQ_DIR" | cut -f1)  파일 수: $(ls -1 "$FASTQ_DIR"/*.fastq.gz 2>/dev/null | wc -l)/16"
[[ $STATUS -eq 0 ]] && log "✅ 전부 성공" || log "❌ 실패한 파일 있음"
exit $STATUS
