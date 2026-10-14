#!/usr/bin/env Rscript
# ============================================================================
# 08_m6a_regulators_GSE113212.R — 加分项②：m6A 调控酶在发现队列的闭环分析
# 1) 23 个 m6A 调控酶表达热图（分组注释）
# 2) 组间差异（老年肥厚 vs 年轻对照，Wilcoxon + BH）
# 3) 调控酶表达 × 铁死亡 GSVA 评分 Spearman 相关
# 4) 调控酶表达 × 巨噬细胞 ssGSEA 评分 Spearman 相关
# 输出：results/m6a_regulators_GSE113212.csv + 2 张图
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(pheatmap); library(patchwork)
})

BASE <- getwd()
R <- file.path(BASE, "results"); F <- file.path(BASE, "figures")
dir.create(F, showWarnings = FALSE, recursive = TRUE)

# ---------- 数据 ----------
expr <- fread(file.path(R, "expr_gene_GSE113212.csv"), data.table = FALSE)
rownames(expr) <- expr$Gene; expr$Gene <- NULL
grp  <- fread(file.path(R, "group_GSE113212.csv"), data.table = FALSE)
grp_vec <- setNames(grp$group, grp$sample)
grp_vec <- factor(grp_vec[colnames(expr)], levels = c("young_control", "elderly_hypertrophied"))
is_case <- grp_vec == "elderly_hypertrophied"

pcd <- fread(file.path(R, "PCD_scores_GSVA_GSE113212.csv"), data.table = FALSE)
ferr_row <- grep("ferroptosis", pcd$Module, ignore.case = TRUE)[1]
ferro_score <- as.numeric(pcd[ferr_row, -1])
names(ferro_score) <- colnames(pcd)[-1]
ferro_score <- ferro_score[colnames(expr)]

imm <- fread(file.path(R, "immune_ssgsea_GSE113212.csv"), data.table = FALSE)
mac_row <- grep("^Macrophages$", imm$CellType)[1]
mac_score <- as.numeric(imm[mac_row, -1])
names(mac_score) <- colnames(imm)[-1]
mac_score <- mac_score[colnames(expr)]

# ---------- 23 个 m6A 调控酶 ----------
regulators <- list(
  Writer = c("METTL3", "METTL14", "METTL16", "WTAP", "VIRMA", "ZC3H13",
             "RBM15", "RBM15B", "CBLL1"),
  Eraser = c("FTO", "ALKBH5", "ALKBH3"),
  Reader = c("YTHDF1", "YTHDF2", "YTHDF3", "YTHDC1", "YTHDC2",
             "IGF2BP1", "IGF2BP2", "IGF2BP3", "HNRNPA2B1", "LRPPRC", "ELAVL1")
)
all_reg <- unlist(regulators)
reg_role <- rep(names(regulators), sapply(regulators, length))
reg_avail <- all_reg %in% rownames(expr)
cat("m6A regulators covered:", sum(reg_avail), "/", length(all_reg), "\n")
reg_expr <- expr[all_reg[reg_avail], , drop = FALSE]
reg_role_df <- data.frame(Role = reg_role[reg_avail], row.names = rownames(reg_expr))

# ---------- 1. 热图 ----------
ann_col <- data.frame(Group = grp_vec, row.names = colnames(expr))
ann_row <- reg_role_df
pheatmap(as.matrix(reg_expr),
         annotation_col = ann_col, annotation_row = ann_row,
         scale = "row", cluster_cols = FALSE,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         filename = file.path(F, "Fig_m6a_regulators_heatmap_GSE113212.png"),
         width = 7, height = 6, main = "m6A regulators (GSE113212)")

# ---------- 2-4. 差异 + 双相关 ----------
res <- do.call(rbind, lapply(rownames(reg_expr), function(g) {
  v <- as.numeric(reg_expr[g, ])
  wt <- wilcox.test(v ~ is_case)
  c1 <- suppressWarnings(cor.test(v, ferro_score, method = "spearman"))
  c2 <- suppressWarnings(cor.test(v, mac_score,    method = "spearman"))
  data.frame(Gene = g, Role = reg_role_df[g, 1],
             mean_case = mean(v[is_case]), mean_ctrl = mean(v[!is_case]),
             delta_log2 = mean(v[is_case]) - mean(v[!is_case]),  # 数据本身已 log2，组间差即 logFC
             P_group = wt$p.value,
             rho_ferroptosis = unname(c1$estimate), P_ferro = c1$p.value,
             rho_macrophage  = unname(c2$estimate), P_macro = c2$p.value)
}))
res$adjP_group <- p.adjust(res$P_group, "BH")
res$adjP_ferro <- p.adjust(res$P_ferro, "BH")
res$adjP_macro <- p.adjust(res$P_macro, "BH")
fwrite(res, file.path(R, "m6a_regulators_GSE113212.csv"))
cat("=== m6A regulators: group / ferroptosis-cor / macrophage-cor ===\n")
print(res[, c("Gene", "Role", "delta_log2", "adjP_group",
              "rho_ferroptosis", "adjP_ferro", "rho_macrophage", "adjP_macro")],
      digits = 3, row.names = FALSE)

# ---------- 相关性棒棒糖图 ----------
df <- res[order(res$rho_ferroptosis), ]
df$sig_f <- ifelse(df$adjP_ferro < 0.05, "adjP<0.05", "ns")
df$sig_m <- ifelse(df$adjP_macro < 0.05, "adjP<0.05", "ns")
p1 <- ggplot(df, aes(x = rho_ferroptosis, y = reorder(Gene, rho_ferroptosis),
                     color = sig_f)) +
  geom_point(size = 3) +
  geom_segment(aes(x = 0, xend = rho_ferroptosis,
                   y = reorder(Gene, rho_ferroptosis),
                   yend = reorder(Gene, rho_ferroptosis)), linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("adjP<0.05" = "#B2182B", "ns" = "grey60")) +
  labs(x = "Spearman rho (regulator expression vs ferroptosis score)",
       y = NULL, title = "m6A regulators x Ferroptosis (GSE113212)") +
  theme_bw(base_size = 11)
p2 <- ggplot(df[order(df$rho_macrophage), ],
             aes(x = rho_macrophage, y = reorder(Gene, rho_macrophage), color = sig_m)) +
  geom_point(size = 3) +
  geom_segment(aes(x = 0, xend = rho_macrophage,
                   y = reorder(Gene, rho_macrophage),
                   yend = reorder(Gene, rho_macrophage)), linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("adjP<0.05" = "#B2182B", "ns" = "grey60")) +
  labs(x = "Spearman rho (regulator expression vs macrophage infiltration)",
       y = NULL, title = "m6A regulators x Macrophages (GSE113212)") +
  theme_bw(base_size = 11)
ggsave(file.path(F, "Fig_m6a_regulators_lollipop_GSE113212.png"),
       p1 | p2, width = 12, height = 7, dpi = 300)

cat("=== 08 done ===\n")
