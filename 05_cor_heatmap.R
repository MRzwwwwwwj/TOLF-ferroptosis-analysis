# ============================================================
# TOLF Step 05: PCD 评分 x 免疫细胞 相关性热图  【论文核心主图】
# 输入: results/PCD_scores_ssGSEA_GSE113212.csv
#       results/immune_cibersort_GSE113212.csv / immune_ssgsea_GSE113212.csv
# 输出: results/PCD_immune_cor.csv
#       figures/Fig_PCD_immune_corheatmap.pdf/png  【主图】
# 方法: Pearson 相关; n=8 样本量小 -> 报告 r 与 P, 讨论中标注探索性;
#       阈值 |r|>0.7 & P<0.05 为强相关(对标 39011365)
# ============================================================

.libPaths("D:/A生信/R-4.5.3/library")
suppressPackageStartupMessages({library(ggplot2); library(reshape2)})

wd <- "C:/Users/10215/WorkBuddy/2026-10-03-00-10-07"

pcd <- read.csv(file.path(wd, "results", "PCD_scores_ssGSEA_GSE113212.csv"),
                row.names = 1, check.names = FALSE)
group <- read.csv(file.path(wd, "results", "group_GSE113212.csv"), stringsAsFactors = FALSE)
pcd <- as.matrix(pcd[, group$sample])

read_immune <- function(f, idcol) {
  d <- read.csv(file.path(wd, "results", f), check.names = FALSE)
  m <- as.matrix(d[, -1])
  rownames(m) <- colnames(d)[-1]
  colnames(m) <- d[[1]]
  m[, group$sample, drop = FALSE]
}

imm_list <- list()
fc <- file.path(wd, "results", "immune_cibersort_GSE113212.csv")
if (file.exists(fc)) imm_list$CIBERSORT <- read_immune("immune_cibersort_GSE113212.csv")
fs <- file.path(wd, "results", "immune_ssgsea_GSE113212.csv")
if (file.exists(fs)) {
  d <- read.csv(fs, row.names = 1, check.names = FALSE)
  imm_list$ssGSEA28 <- as.matrix(d[, group$sample])
}

cor_all <- data.frame()
for (tag in names(imm_list)) {
  imm <- imm_list[[tag]]
  for (m in rownames(pcd)) {
    for (c in rownames(imm)) {
      ct <- suppressWarnings(cor.test(as.numeric(pcd[m, ]), as.numeric(imm[c, ]),
                                      method = "pearson"))
      cor_all <- rbind(cor_all, data.frame(dataset = tag, PCD_module = m,
                                           Immune_cell = c,
                                           r = unname(ct$estimate),
                                           P = ct$p.value))
    }
  }
}
cor_all$strong <- abs(cor_all$r) > 0.7 & cor_all$P < 0.05
write.csv(cor_all, file.path(wd, "results", "PCD_immune_cor.csv"), row.names = FALSE)
cat("总相关对:", nrow(cor_all), " 强相关对(|r|>0.7 & P<0.05):", sum(cor_all$strong), "\n")
print(cor_all[cor_all$strong, ][order(-abs(cor_all$r[cor_all$strong])), ])

# ---------- 主图: 相关性热图 ----------
pick <- if ("CIBERSORT" %in% names(imm_list)) "CIBERSORT" else names(imm_list)[1]
sub <- cor_all[cor_all$dataset == pick, ]
sub$label <- ifelse(sub$P < 0.05, sprintf("%.2f*", sub$r), sprintf("%.2f", sub$r))
p <- ggplot(sub, aes(Immune_cell, PCD_module, fill = r)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient2(low = "#185FA5", mid = "white", high = "#D85A30",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(title = paste0("Correlation: multi-mode PCD scores vs immune infiltration (", pick, ")"),
       x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(wd, "figures", "Fig_PCD_immune_corheatmap.pdf"), p, width = 10, height = 4.5)
ggsave(file.path(wd, "figures", "Fig_PCD_immune_corheatmap.png"), p, width = 10, height = 4.5, dpi = 300)
cat("\n=== 完成 ===\n")
