#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────
# STEP 11 — 시그니처 스코어링 (샘플 단위 활성 점수)
#
# 논문의 핵심 주장을 샘플 단위 점수로 정량합니다:
#   · VSMC contractile → synthetic phenotype switch
#   · Senescence / SASP
#   · Replication stress / DNA damage response
#   · Cell cycle, EMT, 염증, 섬유화
#
# 방법: singscore (rank 기반 단일샘플 스코어링) 우선 → GSVA/ssGSEA → z-score 평균
#       ※ GSVA는 이 aarch64 환경에서 Rhdf5lib 빌드가 실패해 설치되지 않습니다.
#         singscore는 동일하게 샘플별 점수를 주는 표준 방법이라 대체로 씁니다.
#       + limma로 그룹 간 시그니처 차이 검정
# ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
    library(data.table); library(ggplot2); library(pheatmap); library(limma)
    library(msigdbr); library(org.Hs.eg.db); library(AnnotationDbi); library(RColorBrewer)
})
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])),
                 "_plot_setup.R"))

# [deposit] PROJ default = dataset dir (parent of this scripts/ dir; was ~/progerin_rnaseq); overridable by PROJ
PROJ    <- Sys.getenv("PROJ", normalizePath(file.path(dirname(sub("--file=", "",
               grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "..")))
RESULTS <- file.path(PROJ, "results")
OUT     <- file.path(RESULTS, "11_signatures")
FIG     <- file.path(OUT, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
theme_set(theme_bw(base_size = 12))
GRP_COL <- c(Control = "#0072B2", Progerin = "#D55E00")

# STEP 09가 실제로 사용한 샘플 구성(collapse 반영)을 따릅니다.
# 없으면 원본 samples.tsv로 폴백합니다.
.cd <- file.path(RESULTS, "09_deseq2/coldata_used.tsv")
meta <- if (file.exists(.cd)) fread(.cd) else fread(file.path(PROJ, "meta/samples.tsv"))
vst  <- fread(file.path(RESULTS, "09_deseq2/vst_counts.tsv"))
m    <- as.matrix(vst[, -1]); rownames(m) <- vst$gene_id
m    <- m[, meta$sample]

sym <- mapIds(org.Hs.eg.db, keys = rownames(m), keytype = "ENSEMBL",
              column = "SYMBOL", multiVals = "first")
# SYMBOL 행렬 (중복 심볼은 분산이 큰 행 유지)
keep <- !is.na(sym)
ms   <- m[keep, ]; rs <- sym[keep]
ord  <- order(apply(ms, 1, var), decreasing = TRUE)
ms   <- ms[ord, ]; rs <- rs[ord]
ms   <- ms[!duplicated(rs), ]; rownames(ms) <- rs[!duplicated(rs)]
cat(sprintf("SYMBOL 행렬: %d 유전자 x %d 샘플\n", nrow(ms), ncol(ms)))

grp <- factor(meta$group, levels = c("Control", "Progerin"))
BATCH <- if ("batch" %in% names(meta) && uniqueN(meta$batch) > 1)
    factor(meta$batch) else NULL
if (!is.null(BATCH)) cat("배치 blocking 적용:", paste(levels(BATCH), collapse = "/"), "\n")

# ── 시그니처 정의 ─────────────────────────────────────────────
# (1) 문헌 기반 큐레이션 세트 — 이 논문 맥락에 직접 대응
CURATED <- list(
    VSMC_contractile = c("MYH11","ACTA2","CNN1","TAGLN","LMOD1","MYOCD","SRF","SMTN","DES","TPM2"),
    VSMC_synthetic   = c("KLF4","SPP1","VCAM1","COL1A1","COL1A2","COL3A1","COL8A1","FN1","MMP2","MMP9","MMP3","ELN","VIM"),
    Senescence_core  = c("CDKN1A","CDKN2A","CDKN2B","TP53","GLB1","TP53BP1","SERPINE1","IGFBP7"),
    SASP             = c("IL6","IL1A","IL1B","CXCL8","CXCL1","CXCL2","CCL2","CCL20","TNF","IGFBP3",
                         "SERPINE1","MMP1","MMP3","MMP10","VEGFA","CSF2","ICAM1"),
    Proliferation    = c("MKI67","PCNA","CCNA2","CCNB1","CCNE1","CDK1","TOP2A","E2F1","BUB1","AURKA","PLK1"),
    Replication_stress = c("ATR","ATRIP","CHEK1","CLSPN","TIMELESS","TIPIN","RPA1","RPA2","RPA3",
                           "RRM2B","CDC25A","WEE1","FANCD2","BLM","RAD51"),
    DNA_damage_response = c("ATM","CHEK2","TP53","MDM2","GADD45A","BRCA1","BRCA2","RAD50","MRE11","NBN","H2AX"),
    Telomere         = c("TERF1","TERF2","TERF2IP","TINF2","POT1","ACD","TERT","RTEL1"),
    Interferon       = c("ISG15","IFI6","IFI27","IFIT1","IFIT3","MX1","OAS1","STAT1","IRF7","CGAS","STING1"),
    Fibrosis_TGFb    = c("TGFB1","TGFB2","CTGF","CCN2","POSTN","LOX","LOXL1","TIMP1","THBS1","SMAD3"),
    NFkB_target      = c("NFKB1","RELA","NFKBIA","TNFAIP3","BIRC3","TRAF1","CXCL8","IL6","ICAM1"),
    Nuclear_envelope = c("LMNA","LMNB1","LMNB2","LBR","EMD","SUN1","SUN2","SYNE1","SYNE2","BANF1","TMPO")
)

# (2) MSigDB Hallmark — 표준 참조
msig <- function(...) {
    a <- list(...)
    r <- tryCatch(do.call(msigdbr, c(list(species = "Homo sapiens"), a)), error = function(e) NULL)
    if (is.null(r)) {
        names(a)[names(a) == "category"] <- "collection"
        r <- tryCatch(do.call(msigdbr, c(list(species = "Homo sapiens"), a)), error = function(e) NULL)
    }
    r
}
hm <- msig(category = "H")
HALLMARK <- if (!is.null(hm) && nrow(hm))
    split(as.character(hm$gene_symbol), hm$gs_name) else list()
cat(sprintf("Hallmark gene set: %d개\n", length(HALLMARK)))

score_and_test <- function(sets, tag, use_gsva = TRUE) {
    sets <- lapply(sets, function(g) intersect(g, rownames(ms)))
    sets <- sets[lengths(sets) >= 3]
    if (!length(sets)) { cat("유효 세트 없음:", tag, "\n"); return(invisible(NULL)) }

    sc <- NULL; method <- NA_character_

    # (1) singscore — rank 기반, 샘플별 독립 점수 (배치/샘플수에 강건)
    if (is.null(sc) && requireNamespace("singscore", quietly = TRUE)) {
        sc <- tryCatch({
            rk <- singscore::rankGenes(ms)
            t(sapply(sets, function(g)
                singscore::simpleScore(rk, upSet = g, centerScore = TRUE)$TotalScore))
        }, error = function(e) { message("  singscore 실패: ", conditionMessage(e)); NULL })
        if (!is.null(sc)) { colnames(sc) <- colnames(ms); method <- "singscore" }
    }

    # (2) GSVA/ssGSEA — 설치되어 있으면 사용
    if (is.null(sc) && use_gsva && requireNamespace("GSVA", quietly = TRUE)) {
        sc <- tryCatch({
            gp <- GSVA::ssgseaParam(ms, sets)
            GSVA::gsva(gp, verbose = FALSE)
        }, error = function(e) {
            tryCatch(GSVA::gsva(ms, sets, method = "ssgsea", verbose = FALSE),
                     error = function(e2) NULL)
        })
        if (!is.null(sc)) method <- "GSVA/ssGSEA"
    }

    if (!is.null(sc)) cat("  방법:", method, "—", tag, "\n")
    if (is.null(sc)) {
        # 대체: 유전자 z-score 평균 (샘플 간 표준화 후 세트 평균)
        z <- t(scale(t(ms)))
        sc <- t(sapply(sets, function(g) colMeans(z[g, , drop = FALSE], na.rm = TRUE)))
        cat("  방법: z-score 평균 (GSVA 미사용) —", tag, "\n")
    }

    fwrite(data.table(signature = rownames(sc), as.data.frame(sc)),
           file.path(OUT, sprintf("scores_%s.tsv", tag)), sep = "\t")

    # limma 검정
    # 배치가 있으면 blocking factor로 넣습니다. STEP 09에서 확인했듯 이
    # 데이터는 배치가 분산의 75.6%를 차지해서, 빼놓으면 그룹 효과가 묻힙니다.
    des <- if (!is.null(BATCH)) model.matrix(~ BATCH + grp) else model.matrix(~ grp)
    coef_i <- ncol(des)          # grpProgerin (마지막 열)
    fit <- eBayes(lmFit(sc, des))
    tt  <- topTable(fit, coef = coef_i, number = Inf, sort.by = "P")
    tt  <- data.table(signature = rownames(tt), as.data.frame(tt))
    fwrite(tt, file.path(OUT, sprintf("stats_%s.tsv", tag)), sep = "\t")
    cat(sprintf("  유의 시그니처(FDR<0.05): %d/%d\n", sum(tt$adj.P.Val < 0.05), nrow(tt)))
    print(head(tt[, .(signature, logFC, P.Value, adj.P.Val)], 12))

    # 시각화는 배치를 제거한 점수로 그립니다(검정은 위에서 blocking으로 처리).
    sc_vis <- if (!is.null(BATCH)) removeBatchEffect(sc, batch = BATCH,
                                                     design = model.matrix(~ grp)) else sc

    # 히트맵
    ann <- data.frame(group = grp, row.names = colnames(sc))
    if (!is.null(BATCH)) ann$batch <- BATCH
    h <- min(2 + 0.22 * nrow(sc), 20)
    png(file.path(FIG, sprintf("heatmap_%s.png", tag)), 8, h, units = "in", res = 300, type = "cairo")
    pheatmap(t(scale(t(sc_vis))), annotation_col = ann, annotation_colors = list(group = GRP_COL),
             fontsize_row = 7, main = sprintf("시그니처 점수 (%s, z-score)", tag))
    dev.off()

    # 박스플롯 (유의한 것 위주)
    show <- head(tt[order(P.Value), signature], 16)
    long <- rbindlist(lapply(show, function(s)
        data.table(signature = s, score = sc_vis[s, ], group = grp)))
    long[, signature := factor(signature, levels = show)]
    lab <- setNames(sprintf("%s\n(FDR=%.3g)", tt$signature, tt$adj.P.Val), tt$signature)
    png(file.path(FIG, sprintf("boxplot_%s.png", tag)), 12, 9, units = "in", res = 300, type = "cairo")
    print(ggplot(long, aes(group, score, color = group)) +
          geom_boxplot(outlier.shape = NA) + geom_jitter(width = .15, size = 2) +
          facet_wrap(~signature, scales = "free_y", ncol = 4,
                     labeller = labeller(signature = lab)) +
          scale_color_manual(values = GRP_COL) +
          theme(strip.text = element_text(size = 7)) +
          labs(x = NULL, y = "Signature score", title = sprintf("시그니처 점수 — %s", tag)))
    dev.off()

    invisible(tt)
}

cat("\n── 큐레이션 시그니처 ──\n")
tt_cur <- score_and_test(CURATED, "curated")

if (length(HALLMARK)) {
    cat("\n── MSigDB Hallmark ──\n")
    tt_hm <- score_and_test(HALLMARK, "hallmark")
}

# ── VSMC phenotype switch 요약 그림 (논문 핵심 주장) ─────────
sc_cur <- fread(file.path(OUT, "scores_curated.tsv"))
scm <- as.matrix(sc_cur[, -1]); rownames(scm) <- sc_cur$signature
if (!is.null(BATCH)) scm <- removeBatchEffect(scm, batch = BATCH, design = model.matrix(~ grp))
if (all(c("VSMC_contractile", "VSMC_synthetic") %in% rownames(scm))) {
    d <- data.table(sample = colnames(scm), group = grp,
                    contractile = scm["VSMC_contractile", ],
                    synthetic   = scm["VSMC_synthetic", ])
    d[, switch_index := synthetic - contractile]
    fwrite(d, file.path(OUT, "vsmc_phenotype_switch.tsv"), sep = "\t")

    p1 <- ggplot(d, aes(contractile, synthetic, color = group, label = sample)) +
        geom_point(size = 4) + ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
        scale_color_manual(values = GRP_COL) +
        labs(title = "VSMC 표현형 전환", x = "Contractile score", y = "Synthetic score")
    tp <- t.test(switch_index ~ group, data = d)
    p2 <- ggplot(d, aes(group, switch_index, color = group)) +
        geom_boxplot(outlier.shape = NA) + geom_jitter(width = .15, size = 3) +
        scale_color_manual(values = GRP_COL) +
        labs(title = "Phenotype switch index (synthetic − contractile)", x = NULL,
             subtitle = sprintf("Welch t-test p = %.3g", tp$p.value))
    png(file.path(FIG, "vsmc_phenotype_switch.png"), 11, 5, units = "in", res = 300, type = "cairo")
    print(cowplot::plot_grid(p1, p2, nrow = 1, rel_widths = c(1.2, 1)))
    dev.off()
    cat(sprintf("\nVSMC switch index: Welch t-test p = %.4g\n", tp$p.value))
}

cat("\n완료 →", FIG, "\n")
