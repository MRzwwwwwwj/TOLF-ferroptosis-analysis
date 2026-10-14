# ============================================================
# 04b CIBERSORT 免疫浸润 (第二算法, 对标文献 39011365)
# 脚本: CIBERSORT.R (singha53/amritr 镜像, 社区标准版)
# 签名矩阵: LM22 (Newman et al. Nat Methods 2015, 22 种免疫细胞)
# 输出: 22 细胞比例 + 分组比较 + 与 PCD 评分相关性
# ============================================================
wd <- "C:/Users/10215/WorkBuddy/2026-10-03-00-10-07"
setwd(wd)
dir.create("results", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)

suppressPackageStartupMessages({
  library(e1071); library(preprocessCore); library(ggplot2)
})

# ---------- 1. 准备 mixture 文件 ----------
expr <- read.csv("results/expr_gene_GSE113212.csv", row.names = 1, check.names = FALSE)
expr <- as.matrix(expr)
# 去重基因 (保留最大方差)
rv <- apply(expr, 1, var)
expr <- expr[!duplicated(rownames(expr)) | FALSE, ]
genes <- rownames(expr)
keep <- !duplicated(genes)
expr <- expr[keep, ]; rownames(expr) <- genes[keep]
# 与 LM22 基因取交集检查
lm22 <- read.delim("LM22.txt", check.names = FALSE)
ov <- intersect(rownames(expr), lm22[[1]])
cat("LM22 基因覆盖:", length(ov), "/", nrow(lm22), "\n")

mix <- data.frame(Gene = rownames(expr), expr, check.names = FALSE)
write.table(mix, "mixture_GSE113212.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# ---------- 2. 运行 CIBERSORT ----------
source("CIBERSORT.R")
set.seed(42)
cs <- CIBERSORT("LM22.txt", "mixture_GSE113212.txt", perm = 1000, QN = TRUE)
cat("CIBERSORT 完成:", nrow(cs), "样本 x", ncol(cs), "列\n")
cat("去卷积显著性(P<0.05):", sum(cs[, "P-value"] < 0.05), "/", nrow(cs), "\n")
cat("相关性 R 范围:", range(cs[, "Correlation"]), "\n")

cs_df <- as.data.frame(cs, stringsAsFactors = FALSE)
cs_out <- data.frame(Sample = rownames(cs), cs_df, check.names = FALSE)
write.csv(cs_out, "results/immune_cibersort_GSE113212.csv", row.names = FALSE)

# ---------- 3. 分组比较 (wilcoxon) ----------
grp <- read.csv("results/group_GSE113212.csv", stringsAsFactors = FALSE)
gmap <- setNames(grp$Group, grp$Sample)
cells <- colnames(cs)[1:22]
cs <- as.data.frame(cs, stringsAsFactors = FALSE)
cmp <- data.frame()
for (cc in cells) {
  g <- gmap[rownames(cs)]
  g2 <- ifelse(grepl("old|elderly", g, ignore.case = TRUE), "Elderly", "Young")
  v_e <- cs[g2 == "Elderly", cc]; v_y <- cs[g2 == "Young", cc]
  wt <- try(wilcox.test(v_e, v_y, exact = FALSE), silent = TRUE)
  p <- if (inherits(wt, "try-error")) NA else wt$p.value
  cmp <- rbind(cmp, data.frame(Cell = cc, Elderly_mean = mean(v_e), Young_mean = mean(v_y),
                               W_p = p))
}
cmp$W_p_adj <- p.adjust(cmp$W_p, method = "BH")
cmp <- cmp[order(cmp$W_p), ]
write.csv(cmp, "results/immune_cibersort_groupcompare.csv", row.names = FALSE)
cat("\n分组差异 top5 (CIBERSORT):\n")
print(head(cmp, 5))

# ---------- 4. 与 PCD 评分相关性 ----------
pcd <- read.csv("results/PCD_scores_GSVA_GSE113212.csv", row.names = 1, check.names = FALSE)
sm <- intersect(rownames(cs), rownames(pcd))
if (length(sm) >= 4) {
  cs2 <- cs[sm, cells, drop = FALSE]; pcd2 <- pcd[sm, , drop = FALSE]
  cor_m <- matrix(NA, ncol(pcd2), ncol(cs2),
                  dimnames = list(colnames(pcd2), colnames(cs2)))
  pv_m <- matrix(NA, ncol(pcd2), ncol(cs2),
                 dimnames = list(colnames(pcd2), colnames(cs2)))
  for (i in seq_len(ncol(pcd2))) for (j in seq_len(ncol(cs2))) {
    ct <- suppressWarnings(cor.test(pcd2[, i], cs2[, j], method = "pearson"))
    cor_m[i, j] <- ct$estimate; pv_m[i, j] <- ct$p.value
  }
  write.csv(data.frame(PCD = rownames(cor_m), cor_m, check.names = FALSE),
            "results/PCD_immune_cibersort_cor.csv", row.names = FALSE)
  write.csv(data.frame(PCD = rownames(pv_m), pv_m, check.names = FALSE),
            "results/PCD_immune_cibersort_cor_p.csv", row.names = FALSE)
  sig <- which(pv_m < 0.05 & abs(cor_m) > 0.7, arr.ind = TRUE)
  cat("\nPCD x CIBERSORT 强相关对 (|r|>0.7, p<0.05):", nrow(sig), "\n")
  for (k in seq_len(nrow(sig)))
    cat(sprintf("  %s x %s: r=%.3f p=%.4f\n",
                rownames(cor_m)[sig[k,1]], colnames(cor_m)[sig[k,2]],
                cor_m[sig[k,1], sig[k,2]], pv_m[sig[k,1], sig[k,2]]))
}
cat("\nCIBERSORT 分析完成\n")
