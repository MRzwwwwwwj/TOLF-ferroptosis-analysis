# ============================================================
# TOLF Step 05: PCD 评分 x 免疫细胞 相关性热图  【论文核心主图 Figure 1D】
# 输入: results/PCD_scores_ssGSEA_GSE113212.csv
#       results/immune_ssgsea_GSE113212.csv (Bindea 23 种细胞, ssGSEA28)
# 输出: results/PCD_immune_cor.csv
#       figures/Fig_PCD_immune_corheatmap.pdf/png  【主图】
# 方法: Pearson 相关; n=8 样本量小 -> 报告 r 与 P, 讨论中标注探索性;
#       阈值 |r|>0.7 & P<0.05 为强相关(对标 39011365)
# 说明: 主图固定用 ssGSEA28（论文正文报告 115 对相关即出自此矩阵）;
#       CIBERSORT (LM22) 相关矩阵由 04b/04c 单独产出
#       (results/PCD_immune_cibersort_cor*.csv)，不进本脚本，
#       避免 pick 逻辑在不同环境下切换主图数据集。
# 2026-10-15: 模块名防御性去 ".txt" 后缀; 图中模块名首字母大写、
#             细胞名下划线转空格（仅显示层，CSV 保持原始名）。
# ============================================================

.libPaths("D:/ARLARS/R-4.5.3/library")
suppressPackageStartupMessages({library(ggplot2); library(reshape2)})

wd <- "C:/Users/10215/WorkBuddy/2026-10-03-00-10-07"

pcd <- read.csv(file.path(wd, "results", "PCD_scores_ssGSEA_GSE113212.csv"),
                row.names = 1, check.names = FALSE)
rownames(pcd) <- sub("\\.txt$", "", rownames(pcd))  # 防御性清洗
group <- read.csv(file.path(wd, "results", "group_GSE113212.csv"), stringsAsFactors = FALSE)
pcd <- as.matrix(pcd[, group$sample])

d <- read.csv(file.path(wd, "results", "immune_ssgsea_GSE113212.csv"),
              row.names = 1, check.names = FALSE)
imm <- as.matrix(d[, group$sample])

cor_all <- data.frame()
for (m in rownames(pcd)) {
  for (c in rownames(imm)) {
    ct <- suppressWarnings(cor.test(as.numeric(pcd[m, ]), as.numeric(imm[c, ]),
                                    method = "pearson"))
    cor_all <- rbind(cor_all, data.frame(dataset = "ssGSEA28", PCD_module = m,
                                         Immune_cell = c,
                                         r = unname(ct$estimate),
                                         P = ct$p.value))
  }
}
cor_all$strong <- abs(cor_all$r) > 0.7 & cor_all$P < 0.05
write.csv(cor_all, file.path(wd, "results", "PCD_immune_cor.csv"), row.names = FALSE)
cat("总相关对:", nrow(cor_all), " 强相关对(|r|>0.7 & P<0.05):", sum(cor_all$strong), "\n")
print(cor_all[cor_all$strong, ][order(-abs(cor_all$r[cor_all$strong])), ])

# ---------- 主图: 相关性热图 (Figure 1D) ----------
# 仅显示层美化，CSV 保持原始命名
disp <- cor_all
cap1 <- function(s) paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
disp$PCD_module  <- cap1(disp$PCD_module)
disp$Immune_cell <- gsub("_", " ", disp$Immune_cell)
# 固定顺序：模块/细胞按原始行序
disp$PCD_module  <- factor(disp$PCD_module, levels = cap1(rownames(pcd)))
disp$Immune_cell <- factor(disp$Immune_cell, levels = gsub("_", " ", rownames(imm)))

disp$label <- ifelse(disp$P < 0.05, sprintf("%.2f*", disp$r), sprintf("%.2f", disp$r))
p <- ggplot(disp, aes(Immune_cell, PCD_module, fill = r)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient2(low = "#185FA5", mid = "white", high = "#D85A30",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(title = "Correlation: multi-mode PCD scores vs immune infiltration (ssGSEA28)",
       x = NULL, y = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(wd, "figures", "Fig_PCD_immune_corheatmap.pdf"), p, width = 10, height = 4.5)
ggsave(file.path(wd, "figures", "Fig_PCD_immune_corheatmap.png"), p, width = 10, height = 4.5, dpi = 300)
cat("\n=== 完成 ===\n")
