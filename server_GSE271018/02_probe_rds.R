#!/usr/bin/env Rscript
# 02_probe_rds.R — 探查 GSE271018 的 9 个 RDS 对象结构
# 目的：在跑分析前回答 4 个问题：
#   1) 服务器的老 Seurat 能否读入这些 RDS（版本兼容性）？
#   2) 每个对象多大（细胞数/基因数）？
#   3) meta.data 里有没有现成的细胞注释列（有就直接用原文注释，省掉聚类）？
#   4) 有没有 SPP1 亚群/分组相关列？
# 用法：Rscript 02_probe_rds.R  →  查看 probe_report.txt，把内容发回

suppressMessages(library(Seurat))

data_dir <- "GSE271018"
rds_files <- list.files(data_dir, pattern = "\\.rds(\\.gz)?$",
                        full.names = TRUE, ignore.case = TRUE)
if (length(rds_files) == 0) stop("No .rds files found in ", data_dir)

report <- file("probe_report.txt", open = "wt")
sink(report, split = TRUE)

cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")
cat("Files found:", length(rds_files), "\n\n")

for (f in sort(rds_files)) {
  cat("======================================================\n")
  cat("FILE:", basename(f), " (", round(file.size(f)/1e9, 2), "GB )\n")
  obj <- tryCatch(readRDS(f), error = function(e) e)
  if (inherits(obj, "error")) {
    cat("  !! READ FAILED:", conditionMessage(obj), "\n")
    cat("  -> 多半是 Seurat 版本不兼容（v5 对象在老 Seurat 上读不了）。\n")
    cat("     请把此报告发回，我们走降级方案。\n\n")
    next
  }
  cat("  class:", class(obj)[1], "\n")
  if (!inherits(obj, "Seurat")) {
    cat("  !! 不是 Seurat 对象，实际结构:\n")
    print(str(obj, max.level = 1))
    cat("\n")
    next
  }
  cat("  dims:", nrow(obj), "genes x", ncol(obj), "cells\n")
  cat("  assays:", paste(names(obj@assays), collapse = ", "), "\n")
  cat("  assay class:", class(obj@assays[[1]])[1], "\n")   # Assay5 还是 Assay
  cat("  reductions:", paste(names(obj@reductions), collapse = ", "), "\n")
  cat("  default.ident levels (前10):",
      paste(head(levels(Idents(obj)), 10), collapse = ", "), "\n")

  md <- obj@meta.data
  cat("  meta.data columns (", ncol(md), "):\n")
  for (cn in colnames(md)) {
    v <- md[[cn]]
    if (is.factor(v) || is.character(v) || is.logical(v)) {
      lv <- unique(as.character(v))
      show <- if (length(lv) <= 12) paste(lv, collapse = "|") else paste0(length(lv), " levels")
      cat(sprintf("    [cat] %-24s : %s\n", cn, show))
    } else if (is.numeric(v)) {
      cat(sprintf("    [num] %-24s : range %.3g .. %.3g\n", cn, min(v, na.rm = TRUE), max(v, na.rm = TRUE)))
    }
  }

  # 注释列候选
  cand <- grep("cell|type|anno|cluster|ident|label|sub", colnames(md),
               ignore.case = TRUE, value = TRUE)
  cat("  >> annotation candidates:", paste(cand, collapse = ", "), "\n")
  rm(obj); gc(verbose = FALSE)
  cat("\n")
}

sink(); close(report)
cat("[done] probe_report.txt written — 请把此文件内容发回\n")
