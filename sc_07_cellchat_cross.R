# ============================================================================
# sc_07_cellchat_cross.R — CellChat 双队列交叉对照图
# 输入：results/CellChat_cross_cohort_summary.csv（由 L-R 明细表汇总）
# 输出：figures/Fig_cc_cross_cohort.png（A: SPP1 通路总强度; B: SPP1-CD44->Fibroblast）
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

BASE <- "C:/Users/10215/WorkBuddy/2026-10-03-00-10-07"
d <- fread(file.path(BASE, "results", "CellChat_cross_cohort_summary.csv"))
d$group  <- factor(d$group, levels = c("Control", "Disease"))
d$cohort <- factor(d$cohort,
                   levels = c("Cohort1_GSE294458_267819", "Cohort2_GSE271018"),
                   labels = c("Cohort 1 (GSE294458/267819)", "Cohort 2 (GSE271018)"))

cols <- c(Control = "#4DBBD5", Disease = "#E64B35")

pA <- ggplot(d, aes(group, SPP1_pathway_total_prob, fill = group)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = round(SPP1_pathway_total_prob, 3)), vjust = -0.4, size = 3.5) +
  facet_wrap(~ cohort) +
  scale_fill_manual(values = cols) +
  labs(title = "SPP1 pathway total communication probability",
       x = NULL, y = "Total probability") +
  theme_bw(base_size = 12) + theme(legend.position = "none")

pB <- ggplot(d, aes(group, SPP1_CD44_to_Fibroblast_prob, fill = group)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = round(SPP1_CD44_to_Fibroblast_prob, 4)), vjust = -0.4, size = 3.5) +
  facet_wrap(~ cohort) +
  scale_fill_manual(values = cols) +
  labs(title = "SPP1-CD44: SPP1+ Macrophage -> Fibroblast",
       x = NULL, y = "Communication probability") +
  theme_bw(base_size = 12) + theme(legend.position = "none")

png(file.path(BASE, "figures", "Fig_cc_cross_cohort.png"),
    width = 3200, height = 1800, res = 300)
gridExtra::grid.arrange(pA, pB, ncol = 2)
dev.off()
cat("=== sc_07 done ===\n")
