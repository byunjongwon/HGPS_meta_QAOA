options(repos = c(CRAN="https://cloud.r-project.org"), Ncpus = 8, timeout = 3600)
# [deposit] DRY_RUN=1 only reports which packages are missing (installs nothing)
DRY_RUN <- Sys.getenv("DRY_RUN", "0") == "1"
if (!DRY_RUN && !requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
need <- c("org.Hs.eg.db","GSVA","decoupleR","dorothea","progeny","ReactomePA",
          "WGCNA","PCAtools","ggpubr","viridis")
for (p in need) {
  if (!requireNamespace(p, quietly=TRUE)) {
    if (DRY_RUN) { cat("MISSING (dry run, not installed):", p, "\n"); next }   # [deposit]
    cat("\n########## installing", p, "##########\n")
    try(BiocManager::install(p, ask=FALSE, update=FALSE), silent=FALSE)
  } else cat("already:", p, "\n")
}
cat("\n===== 최종 상태 =====\n")
for (p in need) cat(sprintf("%-14s %s\n", p, requireNamespace(p, quietly=TRUE)))
