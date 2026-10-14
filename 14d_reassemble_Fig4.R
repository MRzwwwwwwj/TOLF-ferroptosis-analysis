# ============================================================================
# 14d_reassemble_Fig4.R — Fig4 重拼：修复 panel E 字母压标题（遗留 P3）
# 做法：Fig4E 加顶部白边（原图高度的 8%）后再统一缩放到 1800，
#       panel 字母 "E" 落在白边内，不再压标题。
# 输出：figures/panels/Fig4_scRNA.png（覆盖旧版）
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

# panel E：先加顶部白边再缩放
imgE <- image_read(file.path(SRC, "Fig4E.png"))
info <- image_info(imgE)
pad_top <- round(info$height * 0.09)          # 约 9% 高度留白
imgE <- image_border(imgE, "white", paste0("0x", pad_top))
imgE <- image_scale(imgE, "x1800")

imgs <- list(
  load_h(file.path(SRC, "Fig4A.png")),
  load_h(file.path(SRC, "Fig4B.png")),
  load_h(file.path(SRC, "Fig4C.png")),
  load_h(file.path(SRC, "Fig4D.png")),
  imgE
)
labels <- c("A", "B", "C", "D", "E")

# 复刻 14 号脚本的 assemble 逻辑（3 列）
assemble <- function(imgs, labels, ncol, out, pad = 40) {
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
    hmax <- max(sapply(row_imgs, function(im) image_info(im)$height))
    row_imgs <- lapply(row_imgs, function(im) image_extent(im, paste0(
      image_info(im)$width, "x", hmax), gravity = "north", color = "white"))
    rows[[r]] <- image_append(do.call(c, row_imgs))
  }
  wmax <- max(sapply(rows, function(im) image_info(im)$width))
  rows <- lapply(rows, function(im) image_extent(im, paste0(
    wmax, "x", image_info(im)$height), gravity = "northwest", color = "white"))
  s6 <- image_append(do.call(c, rows), stack = TRUE)
  s6 <- image_border(s6, "white", paste0(pad, "x", pad))
  image_write(s6, out, format = "png", density = 300)
  s6
}

out <- assemble(imgs, labels, 3, file.path(OUT, "Fig4_scRNA.png"))
cat("saved:", file.path(OUT, "Fig4_scRNA.png"),
    image_info(out)$width, "x", image_info(out)$height, "\n")
