# ============================================================
# 04c CIBERSORT 后续: 分组比较 + PCD 相关性 (读取已算好的结果)
# ============================================================
wd <- "C:/Users/10215/WorkBuddy/2026-10-03-00-10-07"
setwd(wd)
suppressPackageStartupMessages(library(ggplot2))

cs <- read.csv("results/immune_cibersort_GSE113212.csv", row.names = 1, check.names = FALSE)
grp <- read.csv("results/group_GSE113212.csv", stringsAsFactors = FALSE)
g2 <- setNames(ifelse(grepl("elderly", grp$group), "Elderly_Hypertrophy", "Young_Control"), grp$sample)
cells <- colnames(cs)[1:22]

# ---------- 分组比较 ----------
cmp <- data.frame()
for (cc in cells) {
  v_e <- cs[names(g2)[g2 == "Elderly_Hypertrophy"], cc]
  v_y <- cs[names(g2)[g2 == "Young_Control"], cc]
  wt <- try(wilcox.test(v_e, v_y, exact = FALSE), silent = TRUE)
  p <- if (inherits(wt, "try-error")) NA else wt$p.value
  cmp <- rbind(cmp, data.frame(Cell = cc, Elderly_mean = mean(v_e), Young_mean = mean(v_y), W_p = p))
}
cmp$W_p_adj <- p.adjust(cmp$W_p, method = "BH")
cmp <- cmp[order(cmp$W_p), ]
write.csv(cmp, "results/immune_cibersort_groupcompare.csv", row.names = FALSE)
cat("分组差异 top8 (CIBERSORT):\n"); print(head(cmp, 8))

# ---------- PCD 相关性 (PCD csv 是 模块x样本, 需转置) ----------
pcd <- read.csv("results/PCD_scores_GSVA_GSE113212.csv", row.names = 1, check.names = FALSE)
pcd <- t(pcd)   # 样本 x 模块
colnames(pcd) <- sub("\\.txt$", "", colnames(pcd))
sm <- intersect(rownames(cs), rownames(pcd))
cat("共同样本:", length(sm), "\n")
cs2 <- cs[sm, cells, drop = FALSE]; pcd2 <- pcd[sm, , drop = FALSE]

cor_m <- matrix(NA, ncol(pcd2), ncol(cs2), dimnames = list(colnames(pcd2), colnames(cs2)))
pv_m <- matrix(NA, ncol(pcd2), ncol(cs2), dimnames = list(colnames(pcd2), colnames(cs2)))
for (i in seq_len(ncol(pcd2))) for (j in seq_len(ncol(cs2))) {
  ct <- suppressWarnings(cor.test(pcd2[, i], cs2[, j], method = "pearson"))
  cor_m[i, j] <- ct$estimate; pv_m[i, j] <- ct$p.value
}
write.csv(data.frame(PCD = rownames(cor_m), cor_m, check.names = FALSE),
          "results/PCD_immune_cibersort_cor.csv", row.names = FALSE)
write.csv(data.frame(PCD = rownames(pv_m), pv_m, check.names = FALSE),
          "results/PCD_immune_cibersort_cor_p.csv", row.names = FALSE)

sig <- which(pv_m < 0.05 & abs(cor_m) > 0.7, arr.ind = TRUE)
cat("\n强相关对 (|r|>0.7, p<0.05):", nrow(sig), "\n")
for (k in seq_len(nrow(sig)))
  cat(sprintf("  %s x %s: r=%.3f p=%.4f\n", rownames(cor_m)[sig[k,1]],
              colnames(cor_m)[sig[k,2]], cor_m[sig[k,1], sig[k,2]], pv_m[sig[k,1], sig[k,2]]))

# ---------- 相关性热图 (CIBERSORT 版) ----------
library(pheatmap)
rn <- rownames(cor_m); cn <- colnames(cor_m)
lab <- matrix("", nrow(cor_m), ncol(cor_m))
lab[which(pv_m < 0.05, arr.ind = TRUE)] <- "*"
lab[which(pv_m < 0.01, arr.ind = TRUE)] <- "**"
rownames(lab) <- rn; colnames(lab) <- cn
pheatmap(cor_m, display_numbers = lab, fontsize_number = 10,
         color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         breaks = seq(-1, 1, length.out = 101),
         cluster_rows = FALSE, cluster_cols = FALSE,
         border_color = "grey90",
         filename = "figures/Fig_PCD_immune_CIBERSORT_corheatmap.png",
         width = 11, height = 4, cellwidth = 14, cellheight = 16)
cat("CIBERSORT 相关性热图已保存\n")
