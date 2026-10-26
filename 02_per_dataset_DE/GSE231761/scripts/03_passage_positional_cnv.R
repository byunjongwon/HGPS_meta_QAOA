#!/usr/bin/env Rscript
# =============================================================================
# GSE231761 종합 분석 (2/2)
#   (F) condition x passage 상호작용 + '가속노화' 가설 검정
#   (G) 위치 기반 / LAD 대리지표 분석 (featureCounts 좌표 활용)
#   (H) CNV·결실다형성 위양성 스크리닝  <- GSTM1 문제의 원시데이터 없는 해법
#   (I) 세포주별 프로파일
# =============================================================================
suppressPackageStartupMessages({
  library(DESeq2); library(ggplot2); library(ggrepel); library(pheatmap)
  library(RColorBrewer); library(matrixStats)
})

# [deposit] paths relative to the deposit root (override with PROGERIA_DEPOSIT)
.this <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
DEPOSIT <- Sys.getenv("PROGERIA_DEPOSIT", normalizePath(file.path(dirname(.this), "..", "..", "..")))
BASE <- file.path(DEPOSIT, "02_per_dataset_DE", "GSE231761")  # [deposit] was ~/car_project/progeria_meta/GSE231761
RAW  <- file.path(BASE, "raw_data"); RES <- file.path(BASE, "results"); FIG <- file.path(RES, "figures")
set.seed(1)

o   <- readRDS(file.path(RES, "GSE231761_objects.rds"))
dds <- o$dds
map <- unique(o$res[!is.na(o$res$gene) & o$res$gene != "", c("ensembl","gene")])
sym <- setNames(map$gene, map$ensembl)
cd  <- as.data.frame(colData(dds)); cond <- cd$condition; pas <- cd$passage
nc  <- counts(dds, normalized = TRUE)

# =============================================================================
# (F) condition x passage 상호작용 + 가속노화 검정
# =============================================================================
cat("############ (F) condition x passage ############\n")
d2 <- dds; design(d2) <- ~ condition * passage
d2 <- DESeq(d2, quiet = TRUE)
cat("사용 가능한 계수:\n"); print(resultsNames(d2))
ix <- grep("condition.*passage|conditionHGPS.passage", resultsNames(d2), value = TRUE)
if (length(ix)) {
  ri <- results(d2, name = ix[1], alpha = .05)
  cat(sprintf("\n상호작용 (질환효과가 passage 에 따라 달라지는가) padj<0.05 : %d\n",
              sum(ri$padj < .05, na.rm = TRUE)))
  di <- as.data.frame(ri); di$gene <- sym[rownames(di)]
  di <- di[!is.na(di$padj), ][order(di$padj[!is.na(di$padj)]), ]
  print(head(di[, c("gene","baseMean","log2FoldChange","padj")], 15), digits = 3, row.names = FALSE)
  write.table(di, file.path(RES, "DE_interaction_condition_x_passage.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
} else cat("상호작용 항 없음\n")

# 층별 질환 효과: P7 에서와 P14 에서
sub7  <- dds[, pas == "P7"];  sub7$condition  <- droplevels(sub7$condition)
sub14 <- dds[, pas == "P14"]; sub14$condition <- droplevels(sub14$condition)
design(sub7) <- ~ condition; design(sub14) <- ~ condition
sub7 <- DESeq(sub7, quiet = TRUE); sub14 <- DESeq(sub14, quiet = TRUE)
r7  <- results(sub7,  contrast = c("condition","HGPS","Control"))
r14 <- results(sub14, contrast = c("condition","HGPS","Control"))
cat(sprintf("\n질환 효과 크기: P7 평균|LFC| = %.3f,  P14 평균|LFC| = %.3f\n",
            mean(abs(r7$log2FoldChange), na.rm = TRUE), mean(abs(r14$log2FoldChange), na.rm = TRUE)))
cat(sprintf("P7 과 P14 의 질환 log2FC 상관 : %.3f\n",
            cor(r7$log2FoldChange, r14$log2FoldChange, use = "complete.obs", method = "spearman")))

# 가속노화 검정: 질환 log2FC 와 계대(P14 vs P7) log2FC 의 상관
rp <- results(dds, contrast = c("passage","P14","P7"))
rc <- results(dds, contrast = c("condition","HGPS","Control"))
ok <- !is.na(rp$log2FoldChange) & !is.na(rc$log2FoldChange) &
      (!is.na(rp$padj) & rp$padj < .05 | !is.na(rc$padj) & rc$padj < .05)
cs <- cor(rc$log2FoldChange[ok], rp$log2FoldChange[ok], method = "spearman")
cat(sprintf("\n=== 가속노화 검정 ===\n질환 log2FC vs 계대 log2FC Spearman = %.3f (n=%d)\n", cs, sum(ok)))
cat(ifelse(cs > 0.2, "  -> 양의 상관: 조로증 변화가 계대노화 방향과 일치 (가속노화 가설 지지)\n",
    ifelse(cs < -0.2, "  -> 음의 상관: 조로증 변화가 계대노화와 반대 방향 (단순 가속노화로 설명 안 됨)\n",
                      "  -> 상관 약함: 두 축이 독립적\n")))

png(file.path(FIG, "disease_vs_passage_lfc.png"), width = 6, height = 5.5, units = "in", res = 300, type = "cairo")
dd <- data.frame(disease = rc$log2FoldChange[ok], passage = rp$log2FoldChange[ok])
print(ggplot(dd, aes(disease, passage)) + geom_point(alpha = .25, size = .8, color = "#4C72B0") +
  geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 0, linetype = 2) +
  geom_smooth(method = "lm", color = "#C44E52", se = TRUE) +
  labs(title = "Disease effect vs replicative passage effect",
       subtitle = sprintf("GSE231761 VSMC;  Spearman = %.3f (n = %d)", cs, sum(ok)),
       x = "log2FC  HGPS vs Control", y = "log2FC  P14 vs P7") + theme_bw(base_size = 11))
dev.off()

# =============================================================================
# (G) 위치 기반 / LAD 대리지표
# =============================================================================
cat("\n\n############ (G) 위치 기반 / LAD 대리지표 ############\n")
f1 <- list.files(RAW, pattern = "^GSM.*\\.txt\\.gz$", full.names = TRUE)[1]
fc <- read.delim(gzfile(f1), skip = 1, check.names = FALSE)[, 1:6]
pos <- data.frame(
  ensembl = fc$Geneid,
  chr   = vapply(strsplit(as.character(fc$Chr),   ";"), function(x) x[1], character(1)),
  start = vapply(strsplit(as.character(fc$Start), ";"), function(x) as.numeric(x[1]), numeric(1)),
  end   = vapply(strsplit(as.character(fc$End),   ";"), function(x) max(as.numeric(x)), numeric(1)),
  stringsAsFactors = FALSE)
main <- paste0("chr", c(1:22, "X"))
pos <- pos[pos$chr %in% main, ]; pos$mid <- (pos$start + pos$end)/2
de <- o$res
d <- merge(de[, c("ensembl","gene","log2FoldChange","padj","sig")], pos, by = "ensembl")
d <- d[!is.na(d$padj), ]; d$chr <- factor(d$chr, levels = main)
cat(sprintf("좌표 있는 검정 유전자 %d (UP %d / DOWN %d)\n", nrow(d), sum(d$sig=="UP"), sum(d$sig=="DOWN")))

# 국소 유전자 밀도 (+-2Mb)
dens <- vapply(seq_len(nrow(d)), function(i) {
  s <- d[d$chr == d$chr[i], ]; sum(abs(s$mid - d$mid[i]) <= 2e6) - 1 }, numeric(1))
d$gene_density <- dens
d$dens_decile <- cut(rank(d$gene_density, ties.method="first"),
                     breaks = quantile(rank(d$gene_density, ties.method="first"), probs=seq(0,1,.1)),
                     include.lowest = TRUE, labels = 1:10)
ds <- data.frame(decile = 1:10,
                 n = as.numeric(table(d$dens_decile)),
                 mean_lfc = as.numeric(tapply(d$log2FoldChange, d$dens_decile, mean)),
                 frac_up  = as.numeric(tapply(d$sig=="UP",   d$dens_decile, mean)),
                 frac_dn  = as.numeric(tapply(d$sig=="DOWN", d$dens_decile, mean)),
                 med_density = as.numeric(tapply(d$gene_density, d$dens_decile, median)))
cat("\n=== 유전자 밀도 십분위 (1=희소=LAD 유사) ===\n"); print(ds, digits=3, row.names=FALSE)
lo <- d$dens_decile %in% c("1","2"); hi <- d$dens_decile %in% c("9","10")
wt <- wilcox.test(d$log2FoldChange[lo], d$log2FoldChange[hi])
cat(sprintf("\n희소(1-2) 평균 log2FC = %.4f  vs  조밀(9-10) = %.4f,  p = %.3g\n",
            mean(d$log2FoldChange[lo]), mean(d$log2FoldChange[hi]), wt$p.value))
cat(sprintf("UP 비율: 희소 %.4f vs 조밀 %.4f\n", mean(d$sig[lo]=="UP"), mean(d$sig[hi]=="UP")))
cat("  (GSE306264 섬유아세포에서는 희소 3.8%% vs 조밀 1.7%%, p=1.6e-05 였다)\n")

# 염색체별
tab <- do.call(rbind, lapply(levels(d$chr), function(cc){
  ic <- d$chr==cc
  data.frame(chr=cc, n=sum(ic), n_up=sum(ic & d$sig=="UP"), n_dn=sum(ic & d$sig=="DOWN"),
             mean_lfc=mean(d$log2FoldChange[ic]),
             p=wilcox.test(d$log2FoldChange[ic], d$log2FoldChange[!ic])$p.value)}))
tab$padj <- p.adjust(tab$p,"BH"); tab <- tab[order(tab$padj), ]
cat("\n=== 염색체별 평균 log2FC ===\n"); print(head(tab,10), digits=3, row.names=FALSE)
write.table(tab, file.path(RES,"positional_chromosome.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(d,   file.path(RES,"positional_gene_table.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

png(file.path(FIG,"LAD_proxy_density_vs_lfc.png"), width=7, height=5, units = "in", res = 300, type = "cairo")
print(ggplot(d, aes(dens_decile, log2FoldChange)) + geom_boxplot(outlier.size=.3, fill="#B0C4DE") +
  geom_hline(yintercept=0, linetype=2, color="red") +
  labs(title="Expression change vs local gene density (LAD proxy)",
       subtitle=sprintf("GSE231761 VSMC; sparse(1-2) %.3f vs dense(9-10) %.3f, p=%.2g",
                        mean(d$log2FoldChange[lo]), mean(d$log2FoldChange[hi]), wt$p.value),
       x="local gene density decile (+-2 Mb)", y="log2FC (HGPS / Control)") +
  coord_cartesian(ylim=quantile(d$log2FoldChange, c(.005,.995))) + theme_bw(base_size=11))
dev.off()

# =============================================================================
# (H) CNV / 결실다형성 위양성 스크리닝   <-- GSTM1 문제
# =============================================================================
cat("\n\n############ (H) CNV·결실다형성 스크리닝 ############\n")
lines <- unique(cd$line)
lm_ <- vapply(lines, function(L) rowMeans(nc[, cd$line == L, drop = FALSE]), numeric(nrow(nc)))
colnames(lm_) <- lines
rownames(lm_) <- sym[rownames(nc)]
lm_ <- lm_[!is.na(rownames(lm_)), ]

present <- lm_ > 50; absent <- lm_ < 3
n_pres <- rowSums(present); n_abs <- rowSums(absent)
cnv_flag <- (n_pres + n_abs == 4) & n_pres >= 1 & n_abs >= 1   # 4개 라인이 깔끔히 有/無로 갈림
cat(sprintf("전 유전자 %d 중 '라인별 全無/全有' 패턴 : %d개 (%.2f%%)\n",
            nrow(lm_), sum(cnv_flag), 100*mean(cnv_flag)))

cnv <- data.frame(gene = rownames(lm_)[cnv_flag], lm_[cnv_flag, , drop = FALSE],
                  n_present = n_pres[cnv_flag], check.names = FALSE)
# 조건과 완전히 일치하는가? (Control 2개 vs HGPS 2개로 갈리면 질환효과와 구분 불가)
ctl <- unique(cd$line[cond=="Control"]); hgp <- unique(cd$line[cond=="HGPS"])
cnv$splits_by_condition <- (rowSums(cnv[, ctl] > 50) == 2 & rowSums(cnv[, hgp] < 3) == 2) |
                           (rowSums(cnv[, ctl] < 3) == 2 & rowSums(cnv[, hgp] > 50) == 2)
cat(sprintf("그중 조건(Control/HGPS)과 완전히 일치하는 것 : %d개  <- DE 로 잘못 잡힘\n",
            sum(cnv$splits_by_condition)))

de2 <- o$res[!is.na(o$res$gene), ]
cnv <- merge(cnv, de2[, c("gene","log2FoldChange","padj","sig")], by="gene", all.x=TRUE)
cnv <- cnv[order(cnv$padj), ]
cat("\n=== 조건과 일치하는 全無/全有 유전자 (= CNV 위양성 강력 후보) ===\n")
print(head(cnv[cnv$splits_by_condition, c("gene", ctl, hgp, "log2FoldChange","padj","sig")], 30),
      digits = 3, row.names = FALSE)
write.table(cnv, file.path(RES, "CNV_polymorphism_candidates.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# 알려진 결실다형성 유전자 점검
known <- c("GSTM1","GSTT1","UGT2B17","UGT2B28","LILRA3","RHD","NAT2","AMY1A","CYP2D6",
           "FCGR3B","BTNL3","BTNL8","LCE3B","LCE3C","HBA1","HBA2","SIGLEC14","GSTM3")
kk <- intersect(known, rownames(lm_))
cat(sprintf("\n=== 알려진 결실다형성 유전자 %d개의 라인별 발현 ===\n", length(kk)))
kt <- data.frame(gene=kk, round(lm_[kk, , drop=FALSE], 1), check.names=FALSE)
kt <- merge(kt, de2[, c("gene","log2FoldChange","padj","sig")], by="gene", all.x=TRUE)
print(kt[order(kt$padj), ], digits=3, row.names=FALSE)
write.table(kt, file.path(RES, "known_CNV_genes.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# DEG 중 CNV 의심 비율
sig_genes <- de2$gene[de2$sig != "NS"]
cat(sprintf("\n유의 DEG %d개 중 CNV 의심(조건일치 全無/全有) : %d개 (%.1f%%)\n",
            length(sig_genes), sum(cnv$gene[cnv$splits_by_condition] %in% sig_genes),
            100*sum(cnv$gene[cnv$splits_by_condition] %in% sig_genes)/length(sig_genes)))
clean <- setdiff(sig_genes, cnv$gene[cnv$splits_by_condition])
cat(sprintf("CNV 의심 제거 후 DEG : %d개\n", length(clean)))
write.table(data.frame(gene=clean), file.path(RES, "DEG_CNVfiltered.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# =============================================================================
# (I) 세포주별 프로파일
# =============================================================================
cat("\n\n############ (I) 세포주별 프로파일 ############\n")
vm <- assay(o$vsd); rownames(vm) <- sym[rownames(vm)]
vm <- vm[!is.na(rownames(vm)) & !duplicated(rownames(vm)), ]
z <- t(scale(t(vm)))
up <- intersect(de2$gene[de2$sig=="UP"], rownames(z)); dn <- intersect(de2$gene[de2$sig=="DOWN"], rownames(z))
sc <- colMeans(z[up,,drop=FALSE]) - colMeans(z[dn,,drop=FALSE])
sdf <- data.frame(sample=names(sc), score=sc, line=cd$line, passage=pas, condition=cond)
sdf <- sdf[order(-sdf$score), ]
cat("=== 샘플별 HGPS 서명 점수 ===\n"); print(sdf, digits=3, row.names=FALSE)
write.table(sdf, file.path(RES,"per_sample_HGPS_score.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

cm <- cor(vm[order(rowVars(vm), decreasing=TRUE)[1:3000], ], method="spearman")
cat("\n=== 세포주 간 상관 ===\n"); print(round(cm,3))
png(file.path(FIG,"sample_correlation.png"), width=6.5, height=5.5, units = "in", res = 300, type = "cairo")
pheatmap(cm, annotation_col=data.frame(condition=cond, passage=pas, row.names=colnames(vm)),
         display_numbers=TRUE, number_format="%.3f", main="Sample correlation (top 3000 variable)",
         col=colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(255)); dev.off()
cat(sprintf("\n저장 완료 -> %s\n", RES))
