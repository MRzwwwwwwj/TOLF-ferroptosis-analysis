# sc_04_cross_validation.R — 两重 scRNA 验证对照汇总
# 输入：scRNA_out/（第一重 GSE294458+GSE267819）与 GSE271018_results/（第二重）
# 输出：figures/Fig_cross_scRNA_validation.png + results/scRNA_cross_validation_summary.csv
# 仅用标准包：ggplot2

suppressPackageStartupMessages(library(ggplot2))

OUT_F <- "figures"
dir.create(OUT_F, showWarnings = FALSE)

# ---------- 读入 ----------
c1_ct  <- read.csv("scRNA_out/sc_ferroptosis_by_celltype.csv")          # celltype,n,mean_score
c1_grp <- read.csv("scRNA_out/sc_ferroptosis_LSS_vs_LDH_wilcox.csv")    # celltype,mean_LSS,mean_LDH,delta,p_adj
c1_sp  <- read.csv("scRNA_out/sc_SPP1macro_ferroptosis.csv")            # celltype_detail,n,mean_score

c2_ct  <- read.csv("GSE271018_results/GSE271018_ferroptosis_by_celltype.csv")  # celltype,n,mean,median
c2_grp <- read.csv("GSE271018_results/GSE271018_deg_vs_ctrl_wilcox.csv")       # celltype,mean_deg,mean_ctrl,delta,p_adj
c2_sp  <- read.csv("GSE271018_results/GSE271018_SPP1macro_ferroptosis.csv")    # spp1_status,n,mean,p_wilcox

# ---------- 细胞类型名对齐 ----------
harm <- c("Macrophage" = "Macrophage", "Fibroblast" = "Fibroblast",
          "Endothelial" = "Endothelial", "Pericyte_SMC" = "Pericyte_SMC",
          "T_NK" = "T_NK", "B_cell" = "B_cell", "Mast_cell" = "Mast_cell",
          "T cell" = "T_NK", "NK cell" = "T_NK", "Mast cell" = "Mast_cell",
          "Pericyte/SMC" = "Pericyte_SMC", "B cell" = "B_cell",
          "Endothelial cell" = "Endothelial")
c1_ct$ct  <- harm[c1_ct$celltype];  c1_ct  <- c1_ct[!is.na(c1_ct$ct), ]
c1_grp$ct <- harm[c1_grp$celltype]; c1_grp <- c1_grp[!is.na(c1_grp$ct), ]

# ---------- Panel A：两队列细胞类型铁死亡均分 ----------
pa_df <- rbind(
  data.frame(ct = c1_ct$ct, mean = c1_ct$mean_score, cohort = "Cohort 1 (GSE294458+267819)"),
  data.frame(ct = c2_ct$celltype, mean = c2_ct$mean, cohort = "Cohort 2 (GSE271018)")
)
common <- intersect(pa_df$ct[pa_df$cohort == levels(factor(pa_df$cohort))[1]],
                    pa_df$ct[pa_df$cohort == levels(factor(pa_df$cohort))[2]])
common <- intersect(c1_ct$ct, c2_ct$celltype)
pa_df <- pa_df[pa_df$ct %in% common, ]
ord <- c2_ct$celltype[order(-c2_ct$mean)]
pa_df$ct <- factor(pa_df$ct, levels = ord)

pA <- ggplot(pa_df, aes(ct, mean, fill = cohort)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(values = c("Cohort 1 (GSE294458+267819)" = "#4393C3",
                               "Cohort 2 (GSE271018)" = "#F4A582")) +
  labs(x = NULL, y = "Mean ferroptosis module score",
       title = "A  Ferroptosis activity by cell type in two independent scRNA cohorts") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "top", legend.title = element_blank())

# ---------- Panel B：病变 vs 对照（巨噬细胞 + 全类型 delta） ----------
b1 <- data.frame(ct = c1_grp$ct, delta = c1_grp$delta,
                 cohort = "Cohort 1 (GSE294458+267819)")
b2 <- data.frame(ct = c2_grp$celltype, delta = c2_grp$delta,
                 cohort = "Cohort 2 (GSE271018)")
pb_df <- rbind(b1, b2)
pb_df$ct <- factor(pb_df$ct, levels = ord)

pB <- ggplot(pb_df, aes(ct, delta, fill = cohort)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("Cohort 1 (GSE294458+267819)" = "#4393C3",
                               "Cohort 2 (GSE271018)" = "#F4A582")) +
  labs(x = NULL, y = expression(Delta * " (diseased - control)"),
       title = "B  Disease-associated ferroptosis increase per cell type") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "top", legend.title = element_blank())

# ---------- Panel C：SPP1+ vs SPP1- 巨噬细胞 ----------
sp1 <- data.frame(status = c1_sp$celltype_detail, mean = c1_sp$mean_score,
                  cohort = "Cohort 1 (GSE294458+267819)",
                  p = c1_sp$p_wilcox[1])
sp2 <- data.frame(status = c2_sp$spp1_status, mean = c2_sp$mean,
                  cohort = "Cohort 2 (GSE271018)",
                  p = c2_sp$p_wilcox[1])
pc_df <- rbind(sp1, sp2)

pC <- ggplot(pc_df, aes(status, mean, fill = status)) +
  geom_col(width = 0.55) +
  facet_wrap(~cohort) +
  scale_fill_manual(values = c("SPP1+_Macrophage" = "#D73027",
                               "SPP1-_Macrophage" = "#74ADD1")) +
  labs(x = NULL, y = "Mean ferroptosis score",
       title = "C  SPP1+ vs SPP1- macrophages") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "none")

# ---------- 拼图输出（无 patchwork 依赖：手动 gridExtra 也不用，分三张存 + 一张组合） ----------
ggsave(file.path(OUT_F, "Fig_cross_sc_A_celltype.png"), pA, width = 8, height = 5, dpi = 300)
ggsave(file.path(OUT_F, "Fig_cross_sc_B_delta.png"),    pB, width = 8, height = 5, dpi = 300)
ggsave(file.path(OUT_F, "Fig_cross_sc_C_SPP1.png"),     pC, width = 7, height = 4.5, dpi = 300)

# ---------- 汇总表 ----------
macro1 <- c1_grp[c1_grp$ct == "Macrophage", ]
macro2 <- c2_grp[c2_grp$celltype == "Macrophage", ]
summary_tab <- data.frame(
  finding = c("Macrophage ferroptosis rank (by mean score)",
              "Macrophage ferroptosis, diseased vs control",
              "Macrophage delta (diseased - control)",
              "Macrophage P value (adj.)",
              "SPP1+ macrophage ferroptosis",
              "SPP1- macrophage ferroptosis",
              "SPP1+ vs SPP1- P value"),
  cohort1 = c(paste0("rank ", which(c1_ct$ct[order(-c1_ct$mean_score)] == "Macrophage"),
                     " of ", length(unique(c1_ct$ct)),
                     " (mean ", round(c1_ct$mean_score[c1_ct$ct == "Macrophage"], 3), ")"),
              paste(round(macro1$mean_LSS, 3), "vs", round(macro1$mean_LDH, 3)),
              round(macro1$delta, 3),
              formatC(macro1$p_adj, format = "e", digits = 1),
              round(c1_sp$mean_score[c1_sp$celltype_detail == "SPP1+_Macrophage"], 3),
              round(c1_sp$mean_score[c1_sp$celltype_detail == "SPP1-_Macrophage"], 3),
              formatC(c1_sp$p_wilcox[1], format = "e", digits = 1)),
  cohort2 = c(paste0("rank ", which(c2_ct$celltype[order(-c2_ct$mean)] == "Macrophage"),
                     " of ", nrow(c2_ct),
                     " (mean ", round(c2_ct$mean[c2_ct$celltype == "Macrophage"], 3), ")"),
              paste(round(macro2$mean_deg, 3), "vs", round(macro2$mean_ctrl, 3)),
              round(macro2$delta, 3),
              formatC(macro2$p_adj, format = "e", digits = 1),
              round(c2_sp$mean[c2_sp$spp1_status == "SPP1+_Macrophage"], 3),
              round(c2_sp$mean[c2_sp$spp1_status == "SPP1-_Macrophage"], 3),
              formatC(c2_sp$p_wilcox[1], format = "e", digits = 1))
)
write.csv(summary_tab, "results/scRNA_cross_validation_summary.csv", row.names = FALSE)
print(summary_tab)
cat("\n[DONE] 图存于 figures/Fig_cross_sc_*.png，表存于 results/scRNA_cross_validation_summary.csv\n")
