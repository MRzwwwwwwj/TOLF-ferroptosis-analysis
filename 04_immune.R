# ============================================================
# TOLF Step 04: 免疫浸润 (双算法交叉验证)
# 算法1: CIBERSORT (IOBR, 22 种免疫细胞, 与 39011365 对标)
# 算法2: ssGSEA (GSVA, 28 种免疫细胞签名, Charoentong 2017)
# 输出: results/immune_cibersort_GSE113212.csv
#       results/immune_ssgsea_GSE113212.csv
#       results/immune_groupcompare.csv
#       figures/Fig_immune_boxplot.pdf/png
# ============================================================

.libPaths("D:/ARLARS/R-4.5.3/library")
suppressPackageStartupMessages({library(GSVA); library(ggplot2); library(reshape2)})

wd <- "C:/Users/10215/WorkBuddy/2026-10-03-00-10-07"
expr <- as.matrix(read.csv(file.path(wd, "results", "expr_gene_GSE113212.csv"),
                           row.names = 1, check.names = FALSE))
group <- read.csv(file.path(wd, "results", "group_GSE113212.csv"), stringsAsFactors = FALSE)
expr <- expr[, group$sample]

# ---------- 算法1: CIBERSORT (IOBR) ----------
cib <- NULL
if (requireNamespace("IOBR", quietly = TRUE)) {
  suppressPackageStartupMessages(library(IOBR))
  # LM22 签名随 IOBR 附带 (lm22)
  data("lm22", package = "IOBR", envir = environment())
  cib <- deconvo_tme(eset = expr, method = "cibersort",
                     arrays = TRUE, perm = 100, absolute = FALSE)
  write.csv(cib, file.path(wd, "results", "immune_cibersort_GSE113212.csv"),
            row.names = FALSE)
  cat("CIBERSORT 完成:", ncol(cib) - 1, "种细胞\n")
} else {
  cat("!! IOBR 未安装, 跳过 CIBERSORT\n")
}

# ---------- 算法2: ssGSEA 免疫细胞签名 ----------
# Bindea et al. 2013 (Immunity) 签名, 从 IOBR signature_collection 导出
# 剔除 4 个非免疫细胞签名(血管/黏膜/肿瘤细胞)
sig_file <- file.path(wd, "immune_signatures_bindea.csv")
sets28 <- NULL
if (file.exists(sig_file)) {
  sig <- read.csv(sig_file, stringsAsFactors = FALSE)
  non_immune <- c("Blood_vessels", "Lymph_vessels", "Normal_mucosa", "SW480_cancer_cells")
  sig <- sig[!sig$celltype %in% non_immune, ]
  sets28 <- split(sig$gene, sig$celltype)
  cat("Bindea 免疫细胞签名:", length(sets28), "种\n")
} else {
  cat("!! 免疫签名文件缺失\n")
}
if (!is.null(sets28)) {
  sp <- ssgseaParam(expr, sets28)
  imm28 <- gsva(sp, verbose = FALSE)
  write.csv(data.frame(CellType = rownames(imm28), imm28, check.names = FALSE),
            file.path(wd, "results", "immune_ssgsea_GSE113212.csv"), row.names = FALSE)
}

# ---------- 组间比较 ----------
cmp <- function(mat, tag) {
  out <- data.frame()
  for (m in rownames(mat)) {
    y <- as.numeric(mat[m, group$group == "young_control"])
    e <- as.numeric(mat[m, group$group == "elderly_hypertrophied"])
    wt <- suppressWarnings(wilcox.test(e, y))
    out <- rbind(out, data.frame(Cell = m, method = tag,
                                 mean_young = mean(y), mean_elderly = mean(e),
                                 wilcox_P = wt$p.value))
  }
  out[order(out$wilcox_P), ]
}
res <- data.frame()
if (!is.null(cib)) {
  m <- as.matrix(cib[, -1]); rownames(m) <- colnames(cib)[-1]
  colnames(m) <- cib[[1]]
  res <- rbind(res, cmp(m, "CIBERSORT"))
}
if (exists("imm28")) res <- rbind(res, cmp(imm28, "ssGSEA28"))
write.csv(res, file.path(wd, "results", "immune_groupcompare.csv"), row.names = FALSE)
cat("\n组间差异免疫细胞 (按 P 排序, 前 15):\n"); print(head(res, 15))

# ---------- 图: 组间有差异(P<0.1)的细胞 ----------
if (exists("imm28")) {
  sigcells <- res$Cell[res$method == "ssGSEA28" & res$wilcox_P < 0.1]
  df <- melt(imm28); colnames(df) <- c("Cell", "sample", "score")
  df$group <- group$group[match(df$sample, group$sample)]
  sub <- df[df$Cell %in% sigcells, ]
  if (nrow(sub) > 0) {
    p <- ggplot(sub, aes(group, score, fill = group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.12, size = 2) +
      scale_fill_manual(values = c("young_control" = "#185FA5",
                                   "elderly_hypertrophied" = "#D85A30")) +
      facet_wrap(~Cell, scales = "free_y") +
      labs(title = "Differentially infiltrated immune cells (ssGSEA, P<0.1)",
           y = "ssGSEA score", x = NULL) +
      theme_bw(base_size = 10) + theme(legend.position = "top")
    ggsave(file.path(wd, "figures", "Fig_immune_boxplot.pdf"), p,
           width = max(6, 2.5 * length(sigcells)), height = 4.5)
    ggsave(file.path(wd, "figures", "Fig_immune_boxplot.png"), p,
           width = max(6, 2.5 * length(sigcells)), height = 4.5, dpi = 300)
  }
}
cat("\n=== 完成 ===\n")
