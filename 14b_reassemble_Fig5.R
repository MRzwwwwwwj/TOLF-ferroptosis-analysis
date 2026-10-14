#!/usr/bin/env Rscript
# 14b_reassemble_Fig5.R — 仅重拼 Fig5（用 09b 重绘的 GSEA 曲线替换空图）
# 拼版逻辑与 14_figure_panels.R 完全一致（magick, tag+append）
suppressPackageStartupMessages({ library(magick) })

BASE <- getwd(); FIG <- file.path(BASE, "figures"); OUT <- file.path(FIG, "panels")
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

# Fig 5：A=相关棒棒糖, B=调控酶热图, C=GSEA 曲线(09b 重绘), D=NES 点图
assemble(
  list(load_h(file.path(FIG, "Fig_m6a_regulators_lollipop_GSE113212.png")),
       load_h(file.path(FIG, "Fig_m6a_regulators_heatmap_GSE113212.png")),
       load_h(file.path(FIG, "Fig_GSEA_ferroptosis.png")),
       load_h(file.path(FIG, "Fig_GSEA_NES_dotplot.png"))),
  c("A", "B", "C", "D"), 2,
  file.path(OUT, "Fig5_m6Areg_GSEA.png"))
cat("=== 14b done ===\n")
