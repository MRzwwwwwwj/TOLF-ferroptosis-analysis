#!/usr/bin/env Rscript
# 14c_reassemble_FigS7.R — 用 monocle_v2（use_partition=FALSE 重跑）结果重拼 FigS7
# 替换旧 v1 分区版拼版（其 D 面板仍标 rho=0.31, P=1e-76）
suppressPackageStartupMessages({ library(magick) })

BASE <- getwd(); OUT <- file.path(BASE, "figures", "panels")
MONO <- file.path(BASE, "scRNA_out", "monocle_v2")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

load_h <- function(p, h = 1800) image_resize(image_read(p), paste0("x", h))
tag <- function(im, lab) image_annotate(im, lab, size = 64, color = "black",
                                        gravity = "northwest", location = "+20+20", weight = 700)
assemble <- function(imgs, labels, ncol, out, cell_h = 1800, pad = 40) {
  tagged <- Map(function(im, lb) tag(im, lb), imgs, labels)
  nrow_ <- ceiling(length(tagged) / ncol)
  rows <- vector("list", nrow_)
  for (r in seq_len(nrow_)) {
    idx <- ((r - 1) * ncol + 1):min(r * ncol, length(tagged))
    row_imgs <- tagged[idx]
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

assemble(
  list(load_h(file.path(MONO, "Fig_mono_A_trajectory.png")),
       load_h(file.path(MONO, "Fig_mono_B_pseudotime.png")),
       load_h(file.path(MONO, "Fig_mono_C_ferroscore.png")),
       load_h(file.path(MONO, "Fig_mono_D_pt_vs_ferro.png")),
       load_h(file.path(MONO, "Fig_mono_E_pt_by_celltype.png")),
       load_h(file.path(MONO, "Fig_mono_F_keygenes.png"))),
  c("A", "B", "C", "D", "E", "F"), 3,
  file.path(OUT, "FigS7_pseudotime.png"))
cat("=== 14c done ===\n")
