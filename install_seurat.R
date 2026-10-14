#!/usr/bin/env Rscript
# 安装单细胞分析所需包（CRAN 二进制，R 4.5 Windows）
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN"))
.libPaths("D:/A生信/R-4.5.3/library")

need <- c("Seurat", "SeuratObject", "Matrix", "future", "harmony",
          "UCell", "dplyr", "ggplot2", "patchwork", "data.table", "R.utils")
have <- rownames(installed.packages())
todo <- setdiff(need, have)
cat("To install:", paste(todo, collapse = ", "), "\n")

if (length(todo) > 0) {
  install.packages(todo, dependencies = NA, quiet = TRUE,
                   lib = "D:/A生信/R-4.5.3/library", type = "binary")
}
# 验证
ok <- sapply(need, function(p) requireNamespace(p, quietly = TRUE))
cat("=== install check ===\n")
print(ok)
if (all(ok)) cat("ALL_OK\n") else cat("MISSING:", paste(names(ok)[!ok], collapse = ","), "\n")
