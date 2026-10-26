#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────
# STEP 9. rMATS 결과 시각화
#
# 산출 그림:
#   01_event_type_bar.pdf      비교별 AS event 종류별 유의 개수
#   02_dpsi_volcano_{비교}.pdf ΔPSI volcano (x=ΔPSI, y=-log10 FDR)
#   03_psi_heatmap_{비교}.pdf  유의 event PSI 히트맵 (샘플 x event)
#   04_dpsi_distribution.pdf   ΔPSI 분포 (비교별, event type별)
#   05_lmna_exon11_PSI.pdf     ★ progerin cryptic splice site
#   06_upset_AS_events.pdf     세 비교 간 유의 event 중첩
#   07_ri_comparison.pdf       intron retention 비교 (노화 지표)
#
# 산출 표:
#   significant_events_{비교}.tsv
#   lmna_events.tsv
# ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
    library(data.table); library(ggplot2); library(ggrepel)
    library(pheatmap); library(RColorBrewer); library(UpSetR)
})

# UTF-8(Δ, 한글) 렌더링을 위한 cairo PDF 장치 설정
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value=TRUE)[1])), "_plot_setup.R"))

# [deposit] default PROJ = this dataset dir (scripts/..), still overridable by PROJ (was ~/lmna_rnaseq)
.this <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
PROJ    <- Sys.getenv("PROJ", normalizePath(file.path(dirname(.this), "..")))
RMATS   <- file.path(PROJ, "results/07_rmats")
OUT     <- file.path(PROJ, "results/09_rmats_figures")
FIG     <- file.path(OUT, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_bw(base_size = 12))
COMPS  <- c("LMNA_vs_WTyng", "WTold_vs_WTyng", "LMNA_vs_WTold")
EVENTS <- c("SE", "MXE", "A5SS", "A3SS", "RI")
FDR_CUT <- 0.05; DPSI_CUT <- 0.1
meta <- fread(file.path(PROJ, "meta/samples.tsv"))

# 비교별 그룹1/그룹2 샘플 순서 (rMATS의 IncLevel1/2 순서와 동일)
grp_samples <- function(g) meta[group == g, sample]
COMP_GRP <- list(
    LMNA_vs_WTyng  = c("LMNA", "WTyng"),
    WTold_vs_WTyng = c("WTold", "WTyng"),
    LMNA_vs_WTold  = c("LMNA", "WTold")
)

# ── 전체 읽기 ─────────────────────────────────────────────────
read_mats <- function(comp, ev) {
    f <- file.path(RMATS, comp, sprintf("%s.MATS.JC.txt", ev))
    if (!file.exists(f)) return(NULL)
    d <- fread(f)
    d[, `:=`(comparison = comp, event = ev)]
    d[, geneSymbol := gsub('"', '', geneSymbol)]
    d
}

all_list <- list()
for (cp in COMPS) for (ev in EVENTS) {
    d <- read_mats(cp, ev)
    if (!is.null(d)) all_list[[paste(cp, ev)]] <- d[, .(
        comparison, event, ID, GeneID, geneSymbol, chr, strand,
        PValue, FDR, IncLevel1, IncLevel2, IncLevelDifference)]
}
all <- rbindlist(all_list, fill = TRUE)
all[, sig := !is.na(FDR) & FDR < FDR_CUT & abs(IncLevelDifference) > DPSI_CUT]
cat(sprintf("전체 event %d개, 유의 %d개\n", nrow(all), sum(all$sig)))

for (cp in COMPS) {
    fwrite(all[comparison == cp & sig == TRUE][order(FDR)],
           file.path(OUT, sprintf("significant_events_%s.tsv", cp)), sep = "\t")
}

# ── 01 event type bar ─────────────────────────────────────────
cnt <- all[sig == TRUE, .N, by = .(comparison, event)]
cnt[, comparison := factor(comparison, levels = COMPS)]
cnt[, event := factor(event, levels = EVENTS)]
pdf(file.path(FIG, "01_event_type_bar.pdf"), 8, 5)
print(ggplot(cnt, aes(event, N, fill = comparison)) +
      geom_col(position = position_dodge(0.8), width = 0.7) +
      geom_text(aes(label = N), position = position_dodge(0.8),
                vjust = -0.3, size = 3) +
      scale_fill_brewer(palette = "Set2") +
      labs(x = "Alternative splicing event type", y = "Significant events",
           title = sprintf("rMATS significant events (FDR<%g, |ΔPSI|>%g)",
                           FDR_CUT, DPSI_CUT)))
dev.off()

# ── 02 ΔPSI volcano ───────────────────────────────────────────
for (cp in COMPS) {
    d <- all[comparison == cp & !is.na(FDR)]
    d[, negLog := -log10(pmax(FDR, 1e-300))]
    d[, lab := ifelse(sig & rank(-negLog) <= 15, geneSymbol, NA_character_)]
    pdf(file.path(FIG, sprintf("02_dpsi_volcano_%s.pdf", cp)), 7.5, 6)
    print(ggplot(d, aes(IncLevelDifference, negLog, color = sig)) +
          geom_point(size = 0.8, alpha = 0.6) +
          geom_vline(xintercept = c(-DPSI_CUT, DPSI_CUT), lty = 2, col = "grey40") +
          geom_hline(yintercept = -log10(FDR_CUT), lty = 2, col = "grey40") +
          geom_text_repel(aes(label = lab), size = 2.8, max.overlaps = 30,
                          show.legend = FALSE) +
          scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "#D55E00"),
                             name = "significant") +
          labs(x = "ΔPSI (IncLevelDifference)", y = "-log10(FDR)", title = cp))
    dev.off()
}

# ── 03 PSI 히트맵 ─────────────────────────────────────────────
psi_matrix <- function(cp, top = 60) {
    d <- all[comparison == cp & sig == TRUE][order(FDR)][seq_len(min(top, .N))]
    if (nrow(d) < 3) return(NULL)
    g1 <- grp_samples(COMP_GRP[[cp]][1]); g2 <- grp_samples(COMP_GRP[[cp]][2])
    m <- t(apply(d, 1, function(r) {
        as.numeric(c(strsplit(r[["IncLevel1"]], ",")[[1]],
                     strsplit(r[["IncLevel2"]], ",")[[1]]))
    }))
    colnames(m) <- c(g1, g2)
    rownames(m) <- make.unique(paste0(d$geneSymbol, "_", d$event))
    m[, colSums(is.na(m)) < nrow(m), drop = FALSE]
}
for (cp in COMPS) {
    m <- psi_matrix(cp)
    if (is.null(m)) next
    keep <- rowSums(is.na(m)) == 0
    m <- m[keep, , drop = FALSE]
    if (nrow(m) < 3) next
    ann <- data.frame(group = meta[match(colnames(m), sample), group],
                      row.names = colnames(m))
    pdf(file.path(FIG, sprintf("03_psi_heatmap_%s.pdf", cp)), 8, 11)
    pheatmap(m, annotation_col = ann, cluster_cols = FALSE,
             color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),
             fontsize_row = 6, main = sprintf("PSI — top significant events (%s)", cp))
    dev.off()
}

# ── 04 ΔPSI 분포 ──────────────────────────────────────────────
pdf(file.path(FIG, "04_dpsi_distribution.pdf"), 10, 5)
print(ggplot(all[sig == TRUE], aes(IncLevelDifference, fill = comparison)) +
      geom_histogram(bins = 60, alpha = 0.7, position = "identity") +
      facet_wrap(~event, scales = "free_y", nrow = 1) +
      scale_fill_brewer(palette = "Set2") +
      labs(x = "ΔPSI", y = "events", title = "ΔPSI distribution of significant events"))
dev.off()

# ── 05 ★ Lmna exon 11 (progerin) ──────────────────────────────
lmna <- all[geneSymbol == "Lmna"]
fwrite(lmna, file.path(OUT, "lmna_events.tsv"), sep = "\t")

lm_top <- lmna[event == "A5SS" & sig == TRUE][order(FDR)]
if (nrow(lm_top) > 0) {
    pl <- rbindlist(lapply(seq_len(nrow(lm_top)), function(i) {
        r  <- lm_top[i]
        cp <- r$comparison
        g1 <- grp_samples(COMP_GRP[[cp]][1]); g2 <- grp_samples(COMP_GRP[[cp]][2])
        data.table(
            comparison = cp,
            group = c(rep(COMP_GRP[[cp]][1], length(g1)), rep(COMP_GRP[[cp]][2], length(g2))),
            sample = c(g1, g2),
            PSI = as.numeric(c(strsplit(r$IncLevel1, ",")[[1]],
                               strsplit(r$IncLevel2, ",")[[1]])),
            label = sprintf("%s\nΔPSI=%.3f  FDR=%.2g", cp, r$IncLevelDifference, r$FDR))
    }))
    pdf(file.path(FIG, "05_lmna_exon11_PSI.pdf"), 4 + 2.2 * nrow(lm_top), 5)
    print(ggplot(pl[!is.na(PSI)], aes(group, PSI, color = group)) +
          geom_boxplot(outlier.shape = NA, width = 0.55) +
          geom_jitter(width = 0.12, size = 3) +
          facet_wrap(~label, scales = "free_x") +
          scale_color_manual(values = c(LMNA = "#D55E00", WTyng = "#0072B2",
                                        WTold = "#009E73")) +
          ylim(0, 1.05) +
          labs(x = NULL, y = "PSI (full-length exon inclusion)",
               title = "Lmna exon 11 alternative 5' splice site",
               subtitle = paste("PSI 1.0 = 정상 스플라이싱 100%.",
                                "PSI 저하 = cryptic donor 사용 (progerin, 150 nt/50 aa 결손)")))
    dev.off()
    cat("★ Lmna A5SS 유의 event", nrow(lm_top), "개 → 05_lmna_exon11_PSI.pdf\n")
}

# ── 06 UpSet ──────────────────────────────────────────────────
key <- function(d) paste(d$event, d$GeneID, d$ID, sep = "|")
sets <- lapply(COMPS, function(cp) key(all[comparison == cp & sig == TRUE]))
names(sets) <- COMPS
if (all(lengths(sets) > 0)) {
    pdf(file.path(FIG, "06_upset_AS_events.pdf"), 8, 6)
    print(upset(fromList(sets), order.by = "freq", text.scale = 1.3))
    dev.off()
}

# ── 07 intron retention (노화 지표) ───────────────────────────
ri <- all[event == "RI" & sig == TRUE]
if (nrow(ri) > 0) {
    ri[, direction := ifelse(IncLevelDifference > 0, "retention ↑", "retention ↓")]
    s <- ri[, .N, by = .(comparison, direction)]
    s[, comparison := factor(comparison, levels = COMPS)]
    pdf(file.path(FIG, "07_ri_comparison.pdf"), 7, 5)
    print(ggplot(s, aes(comparison, N, fill = direction)) +
          geom_col(position = "dodge", width = 0.7) +
          geom_text(aes(label = N), position = position_dodge(0.7), vjust = -0.3, size = 3.2) +
          scale_fill_manual(values = c("retention ↑" = "#CC79A7", "retention ↓" = "#56B4E9")) +
          labs(x = NULL, y = "significant RI events",
               title = "Intron retention — 노화·조로증 비교",
               subtitle = "노화에서 intron retention 증가가 문헌에 반복 보고됨"))
    dev.off()
}

cat("완료 → ", FIG, "\n")
