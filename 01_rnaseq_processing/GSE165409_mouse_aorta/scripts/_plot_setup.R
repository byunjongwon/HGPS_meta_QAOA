# ─────────────────────────────────────────────────────────────
# 그림 공통 설정 — 모든 R 스크립트가 맨 위에서 source 합니다.
#
# 기본 pdf() 장치는 8비트 인코딩만 처리해서 그리스 문자(Δ)나 한글이
# "conversion failure in 'mbcsToSbcs'" 경고와 함께 깨져 들어갑니다.
# cairo 기반 장치로 바꾸면 UTF-8과 시스템 폰트를 그대로 씁니다.
# ─────────────────────────────────────────────────────────────

# 여러 페이지를 한 파일에 담으므로 onefile=TRUE 가 필수입니다.
pdf <- function(file, width = 7, height = 7, ...) {
    grDevices::cairo_pdf(filename = file, width = width, height = height,
                         onefile = TRUE, ...)
}

# 한글이 있는 시스템에서는 CJK 폰트를 기본으로 잡아 둡니다.
.pick_font <- function() {
    fams <- tryCatch(system2("fc-list", c(":lang=ko", "family"), stdout = TRUE),
                     error = function(e) character(0))
    for (want in c("Noto Sans CJK KR", "NanumGothic", "나눔고딕", "Noto Sans CJK SC")) {
        if (any(grepl(want, fams, fixed = TRUE))) return(want)
    }
    ""
}
CJK_FONT <- .pick_font()
if (nzchar(CJK_FONT) && requireNamespace("ggplot2", quietly = TRUE)) {
    # 스크립트가 나중에 theme_bw()를 다시 호출해도 폰트가 풀리지 않도록
    # 전역 환경에 base_family가 채워진 theme_bw를 정의해 가려둡니다.
    theme_bw <- function(base_size = 11, base_family = CJK_FONT, ...)
        ggplot2::theme_bw(base_size = base_size, base_family = base_family, ...)
    ggplot2::theme_set(theme_bw(base_size = 12))
}
# 주의: par(family = CJK_FONT) 를 전역으로 걸면 안 됩니다.
# 패키지 내부(satuRn 등)가 base graphics로 hist()를 그릴 때 그 장치가
# 해당 폰트를 모르면 "유효하지 않은 폰트타입입니다" 에러로 죽습니다.
# 내가 그리는 base 그래픽에는 한글이 없으므로 ggplot 설정만으로 충분합니다.
