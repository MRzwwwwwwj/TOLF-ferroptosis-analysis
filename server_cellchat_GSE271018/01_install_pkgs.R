#!/usr/bin/env Rscript
# ============================================================================
# 01_install_pkgs.R — 服务器安装 CellChat（独立个人库，不更新/不覆盖任何系统旧包）
# 用法：Rscript 01_install_pkgs.R 2>&1 | tee install_log.txt
# 原理：所有新包装进 ~/Rlibs_new；remotes 的 upgrade="never" 保证不升级已装包
# ============================================================================
LIB <- path.expand("~/Rlibs_new")
dir.create(LIB, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(LIB, .libPaths()))
cat("安装目标库:", LIB, "\n\n")

# remotes（若系统已有则直接复用，不升级）
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes", lib = LIB, repos = "https://cloud.r-project.org")

# CellChat 本体（依赖自动装进 LIB；upgrade="never" = 不动系统旧包）
remotes::install_github("sqjin/CellChat", lib = LIB,
                        upgrade = "never", dependencies = TRUE)

# 验证
library(CellChat, lib.loc = c(LIB, .libPaths()))
cat("\n=== CellChat", as.character(packageVersion("CellChat")), "安装成功 ===\n")
cat("数据库文件位置:", system.file("extdata", package = "CellChat"), "\n")
