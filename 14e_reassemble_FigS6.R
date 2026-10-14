# ============================================================================
# 14e_reassemble_FigS6.R — FigS6 重拼：H 面板换用去字母前缀的新跨队列图
# 复刻 14_figure_panels.R 的 FigS6 段（row1: A-D / row2: E-G / row3: H 独占整行）
# 输出：figures/panels/FigS6_CellChat.png（覆盖旧版）
# ============================================================================
suppressPackageStartupMessages(library(magick))

BASE <- "C:/Users/10215/WorkBuddy/2026-10-03-00-10-07"
FIG  <- file.path(BASE, "figures")
SRC  <- file.path(FIG, "panels_src")
OUT  <- file.path(FIG, "panels")

load_h <- function(path, h = 1800) image_scale(image_read(path), paste0("x", h))
tag <- function(img, label) {
  image_annotate(img, label, size = 90, weight = 700,
                 gravity = "northwest", location = "+30+20", color = "black")
}
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
