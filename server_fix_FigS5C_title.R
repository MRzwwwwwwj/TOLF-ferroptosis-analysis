#!/usr/bin/env Rscript
# ============================================================================
# fix_FigS5C_title.R — 只重画 FigS5C（标题 P 值指数被裁：1.78e-184 只显示到 e-1）
# 根因：03_ferroptosis_validate.R 中 ggsave width=4.5 放不下长标题
# 做法：直接读 GSE271018_percell_scores.csv（27,852 行，无需重载 9 个 RDS），
#       同样的 wilcox + ggplot 代码，width 4.5 -> 5.5
# 输出：~/GSE271018_results/FigS5C_SPP1macro.pdf（覆盖）
# ============================================================================
suppressPackageStartupMessages(library(ggplot2))

OUT_RES <- "~/GSE271018_results"
all <- read.csv(file.path(OUT_RES, "GSE271018_percell_scores.csv"),
                stringsAsFactors = FALSE)
cat("total cells:", nrow(all), "\n")

sp <- all[!is.na(all$spp1_status) & all$spp1_status != "", ]
cat("SPP1 cells:", nrow(sp), "\n")
p_sp <- wilcox.test(ferro_score ~ spp1_status, data = sp)$p.value
cat(sprintf("p_sp = %.2e\n", p_sp))

p3 <- ggplot(sp, aes(spp1_status, ferro_score, fill = spp1_status)) +
  geom_boxplot(outlier.size = 0.1) +
  scale_fill_manual(values = c("SPP1+_Macrophage" = "#D73027",
                               "SPP1-_Macrophage" = "#74ADD1")) +
  labs(x = NULL, y = "Ferroptosis score",
       title = sprintf("SPP1+ vs SPP1- macrophages (P = %.2e)", p_sp)) +
  theme_bw(base_size = 12) + theme(legend.position = "none")
ggsave(file.path(OUT_RES, "FigS5C_SPP1macro.pdf"), p3,
       width = 5.5, height = 5)
cat("=== fix_FigS5C done ===\n")
