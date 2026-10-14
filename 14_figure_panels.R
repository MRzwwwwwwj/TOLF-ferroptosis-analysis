# ============================================================================
# 14_figure_panels.R — 论文主图/附图拼版（magick 包，300dpi PNG）
# 输入：figures/*.png + figures/panels_src/*.png（PDF 已转 PNG）
# 输出：figures/panels/Fig1~Fig8 + FigS5~S7 拼版图
# 说明：仙桃 CLI 无拼版工具（2026-10-11 实测 tool search），按项目约定用 R 包
# ============================================================================
suppressPackageStartupMessages(library(magick))

BASE <- "C:/Users/10215/WorkBuddy/2026-10-03-00-10-07"
FIG  <- file.path(BASE, "figures")
SRC  <- file.path(FIG, "panels_src")
OUT  <- file.path(FIG, "panels")
dir.create(OUT, showWarnings = FALSE)

# ---------- 工具函数 ----------
# 读图并统一高度（保持宽高比）
load_h <- function(path, h = 1800) {
  img <- image_read(path)
  image_scale(img, paste0("x", h))
}
# 加 panel 标签（A/B/C...）
tag <- function(img, label) {
  image_annotate(img, label, size = 90, weight = 700,
                 gravity = "northwest", location = "+30+20", color = "black")
}
# 白底占位图（机制示意图预留位）
placeholder <- function(w = 1800, h = 1800, txt = "C\nMechanism diagram\n(to be inserted)") {
  img <- image_blank(w, h, color = "white")
  img <- image_border(img, "grey70", "8x8")
  image_annotate(img, txt, size = 70, gravity = "center", color = "grey40")
}
# 拼版：imgs 为 magick 向量列表，labels 为标签，ncol 列数
assemble <- function(imgs, labels, ncol, out, cell_h = 1800, pad = 40) {
  tagged <- Map(function(im, lb) tag(im, lb), imgs, labels)
  nrow_ <- ceiling(length(tagged) / ncol)
  rows <- vector("list", nrow_)
  for (r in seq_len(nrow_)) {
    idx <- ((r - 1) * ncol + 1):min(r * ncol, length(tagged))
    row_imgs <- tagged[idx]
    # 不足整行补白图，保证对齐
    if (length(row_imgs) < ncol) {
      w_ref <- image_info(row_imgs[[1]])$width
      h_ref <- image_info(row_imgs[[1]])$height
      for (k in seq_len(ncol - length(row_imgs)))
        row_imgs[[length(row_imgs) + 1]] <- image_blank(w_ref, h_ref, "white")
    }
    rows[[r]] <- image_append(do.call(c, row_imgs), stack = FALSE)
  }
  final <- image_append(do.call(c, rows), stack = TRUE)
  final <- image_border(final, "white", paste0(pad, "x", pad))
  image_write(final, out, format = "png", density = 300)
  cat("saved:", basename(out), image_info(final)$width, "x", image_info(final)$height, "\n")
}

# ---------- Fig 1：发现队列总览（火山 + PCD + 免疫 + 相关热图） ----------
assemble(
  list(load_h(file.path(FIG, "Fig_volcano_GSE113212.png")),
       load_h(file.path(FIG, "Fig_PCD_scores_boxplot.png")),
       load_h(file.path(FIG, "Fig_immune_boxplot.png")),
       load_h(file.path(FIG, "Fig_PCD_immune_corheatmap.png"))),
  c("A", "B", "C", "D"), 2,
  file.path(OUT, "Fig1_discovery_overview.png"))

# ---------- Fig 2：富集分析（GO/KEGG/Reactome 上调） ----------
assemble(
  list(load_h(file.path(FIG, "Fig_enrich_GO_BP_DEG_up.png")),
       load_h(file.path(FIG, "Fig_enrich_KEGG_DEG_up.png")),
       load_h(file.path(FIG, "Fig_enrich_REAC_DEG_up.png"))),
  c("A", "B", "C"), 3,
  file.path(OUT, "Fig2_enrichment.png"))

# ---------- Fig 3：外部验证 + m6A + 机制示意（占位） ----------
assemble(
  list(load_h(file.path(FIG, "Fig_PCD_scores_GSE206560.png")),
       load_h(file.path(FIG, "Fig_m6a_peakcounts_GSE206560.png")),
       placeholder()),
  c("A", "B", "C"), 3,
  file.path(OUT, "Fig3_validation_m6A.png"))

# ---------- Fig 4：单细胞五联 ----------
assemble(
  list(load_h(file.path(SRC, "Fig4A.png")),
       load_h(file.path(SRC, "Fig4B.png")),
       load_h(file.path(SRC, "Fig4C.png")),
       load_h(file.path(SRC, "Fig4D.png")),
       load_h(file.path(SRC, "Fig4E.png"))),
  c("A", "B", "C", "D", "E"), 3,
  file.path(OUT, "Fig4_scRNA.png"))

# ---------- Fig 5：m6A 调控酶 + GSEA ----------
assemble(
  list(load_h(file.path(FIG, "Fig_m6a_regulators_lollipop_GSE113212.png")),
       load_h(file.path(FIG, "Fig_m6a_regulators_heatmap_GSE113212.png")),
       load_h(file.path(FIG, "Fig_GSEA_ferroptosis.png")),
       load_h(file.path(FIG, "Fig_GSEA_NES_dotplot.png"))),
  c("A", "B", "C", "D"), 2,
  file.path(OUT, "Fig5_m6Areg_GSEA.png"))

# ---------- Fig 6：LASSO + 三队列 ROC ----------
assemble(
  list(load_h(file.path(FIG, "Fig_LASSO_cvcurve.png")),
       load_h(file.path(FIG, "Fig_ROC_LASSO_PCD.png")),
       load_h(file.path(FIG, "Fig_ROC_LASSO_external.png")),
       load_h(file.path(FIG, "Fig_ROC_LASSO_GSE106253.png"))),
  c("A", "B", "C", "D"), 2,
  file.path(OUT, "Fig6_LASSO_ROC.png"))

# ---------- Fig 7：检查点 + M1/M2 ----------
assemble(
  list(load_h(file.path(FIG, "Fig_ckpt_A_ckpt_corheatmap.png")),
       load_h(file.path(FIG, "Fig_ckpt_B_M1M2_cor.png"))),
  c("A", "B"), 2,
  file.path(OUT, "Fig7_checkpoint_M1M2.png"))

# ---------- Fig 8：Nomogram 三件套 ----------
assemble(
  list(load_h(file.path(FIG, "Fig_nomo_A_nomogram.png")),
       load_h(file.path(FIG, "Fig_nomo_B_calibration.png")),
       load_h(file.path(FIG, "Fig_nomo_C_DCA_train.png")),
       load_h(file.path(FIG, "Fig_nomo_D_DCA_external.png"))),
  c("A", "B", "C", "D"), 2,
  file.path(OUT, "Fig8_nomogram_DCA.png"))

# ---------- Fig S5：GSE271018 复现（A-D 服务器 + E-G 跨队列） ----------
assemble(
  list(load_h(file.path(SRC, "FigS5A.png")),
       load_h(file.path(SRC, "FigS5B.png")),
       load_h(file.path(SRC, "FigS5C.png")),
       load_h(file.path(SRC, "FigS5D.png")),
       load_h(file.path(FIG, "Fig_cross_sc_A_celltype.png")),
       load_h(file.path(FIG, "Fig_cross_sc_B_delta.png")),
       load_h(file.path(FIG, "Fig_cross_sc_C_SPP1.png"))),
  c("A", "B", "C", "D", "E", "F", "G"), 4,
  file.path(OUT, "FigS5_GSE271018.png"))

# ---------- Fig S6：CellChat 双队列（A-F 队列1 + G 队列2 + H 对照独占整行） ----------
r1 <- assemble_rows <- NULL  # 手工分行拼装
mk_row <- function(paths, labels, cell_h = 1800) {
  imgs <- Map(function(p, lb) tag(load_h(p, cell_h), lb), paths, labels)
  image_append(do.call(c, imgs), stack = FALSE)
}
row1 <- mk_row(list(file.path(FIG, "..", "scRNA_out", "cellchat", "Fig_cc_A_compare_interactions.png"),
                    file.path(FIG, "..", "scRNA_out", "cellchat", "Fig_cc_B_diff_bubble.png"),
                    file.path(FIG, "..", "scRNA_out", "cellchat", "Fig_cc_C_rankNet.png"),
                    file.path(SRC, "Fig_cc_D.png")),
               c("A", "B", "C", "D"))
row2 <- mk_row(list(file.path(SRC, "Fig_cc_E.png"),
                    file.path(FIG, "..", "scRNA_out", "cellchat", "Fig_cc_F_SPP1_LR_pairs.png"),
                    file.path(SRC, "Fig_cc_E_srv.png")),
               c("E", "F", "G"))
w_target <- max(image_info(row1)$width, image_info(row2)$width)
row1 <- image_scale(row1, paste0(w_target, "x"))
row2 <- image_scale(row2, paste0(w_target, "x"))
row3 <- tag(image_scale(image_read(file.path(FIG, "Fig_cc_cross_cohort.png")),
                        paste0(w_target, "x")), "H")
s6 <- image_append(c(row1, row2, row3), stack = TRUE)
s6 <- image_border(s6, "white", "40x40")
image_write(s6, file.path(OUT, "FigS6_CellChat.png"), format = "png", density = 300)
cat("saved: FigS6_CellChat.png", image_info(s6)$width, "x", image_info(s6)$height, "\n")

# ---------- Fig S7：monocle3 拟时序六联（2026-10-14 修正：monocle -> monocle_v2 生产路径） ----------
assemble(
  list(load_h(file.path(FIG, "..", "scRNA_out", "monocle_v2", "Fig_mono_A_trajectory.png")),
       load_h(file.path(FIG, "..", "scRNA_out", "monocle_v2", "Fig_mono_B_pseudotime.png")),
       load_h(file.path(FIG, "..", "scRNA_out", "monocle_v2", "Fig_mono_C_ferroscore.png")),
       load_h(file.path(FIG, "..", "scRNA_out", "monocle_v2", "Fig_mono_D_pt_vs_ferro.png")),
       load_h(file.path(FIG, "..", "scRNA_out", "monocle_v2", "Fig_mono_E_pt_by_celltype.png")),
       load_h(file.path(FIG, "..", "scRNA_out", "monocle_v2", "Fig_mono_F_keygenes.png"))),
  c("A", "B", "C", "D", "E", "F"), 3,
  file.path(OUT, "FigS7_pseudotime.png"))

cat("=== 14 done ===\n")
