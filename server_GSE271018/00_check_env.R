#!/usr/bin/env Rscript
# 00_check_env.R — GSE271018 服务器环境验证
# 目的：在任何分析之前确认 R 与关键包的版本是否满足最低要求。
# 用法：Rscript 00_check_env.R  →  查看生成的 env_report.txt
# 不安装、不更新任何包，只读检查。

options(warn = 1)
report <- file("env_report.txt", open = "wt")
sink(report, split = TRUE)   # 同时输出到屏幕和文件

cat("==============================================\n")
cat(" GSE271018 environment check\n")
cat(" date:", as.character(Sys.time()), "\n")
cat("==============================================\n\n")

## ---- 1. R 版本 ----
cat("## R version\n")
cat(R.version.string, "\n")
r_ok <- getRversion() >= "3.6.0"
cat("R >= 3.6 required:", ifelse(r_ok, "PASS", "FAIL"), "\n\n")

## ---- 2. 必需包检查 ----
# 分级：必需 / 可选（没有也能跑，功能降级）
required <- c("Seurat", "ggplot2", "Matrix")
optional <- c("dplyr", "data.table", "patchwork", "future", "harmony", "RColorBrewer")

cat("## Package availability & versions\n")
pkg_info <- data.frame(package = character(), installed = logical(),
                       version = character(), level = character(),
                       stringsAsFactors = FALSE)
for (p in c(required, optional)) {
  lvl <- ifelse(p %in% required, "REQUIRED", "optional")
  ok <- requireNamespace(p, quietly = TRUE)
  ver <- if (ok) as.character(packageVersion(p)) else "NOT INSTALLED"
  pkg_info <- rbind(pkg_info, data.frame(package = p, installed = ok,
                                         version = ver, level = lvl))
  cat(sprintf("  %-14s %-9s %s  (%s)\n", p, ifelse(ok, "OK", "MISSING"), ver, lvl))
}
cat("\n")

## ---- 3. Seurat 版本分支预判 ----
cat("## Seurat API compatibility\n")
if (requireNamespace("Seurat", quietly = TRUE)) {
  sv <- packageVersion("Seurat")
  cat("Seurat version:", as.character(sv), "\n")
  major <- sv$major
  if (major >= 5) {
    cat("API path: v5 (Assay5/layers) — script will use layer= argument + JoinLayers fallback\n")
  } else if (major == 4) {
    cat("API path: v4 (Assay/slot) — fully supported\n")
  } else if (major == 3) {
    cat("API path: v3 (Assay/slot) — fully supported (minimum viable)\n")
  } else {
    cat("API path: Seurat < 3.0 — TOO OLD, AddModuleScore may misbehave; report back before proceeding\n")
  }
  has_ams <- "AddModuleScore" %in% getNamespaceExports("Seurat")
  cat("AddModuleScore available:", has_ams, "\n")
} else {
  cat("Seurat NOT INSTALLED — cannot proceed with Seurat-based analysis.\n")
  cat("Fallback: report back; we can switch to a pure-Matrix + manual scoring pipeline.\n")
}
cat("\n")

## ---- 4. 系统资源 ----
cat("## System resources\n")
mem_line <- tryCatch({
  m <- readLines("/proc/meminfo")
  total <- as.numeric(sub("\\D*(\\d+).*", "\\1", grep("MemTotal", m, value = TRUE))) / 1024^2
  avail <- as.numeric(sub("\\D*(\\d+).*", "\\1", grep("MemAvailable", m, value = TRUE))) / 1024^2
  sprintf("RAM total: %.1f GB, available: %.1f GB", total, avail)
}, error = function(e) "could not read /proc/meminfo")
cat(mem_line, "\n")
disk <- tryCatch(system("df -h . | tail -1", intern = TRUE), error = function(e) NA)
if (!is.na(disk[1])) cat("disk:", disk, "\n")
cpu <- tryCatch(length(parallel::detectCores()), error = function(e) NA)
cat("cores:", cpu, "\n\n")

## ---- 5. 结论 ----
cat("## Verdict\n")
req_ok <- all(sapply(required, requireNamespace, quietly = TRUE))
if (req_ok && r_ok) {
  cat("ENV OK — proceed to 01_download.sh\n")
} else {
  cat("ENV INCOMPLETE — send this report back before downloading data.\n")
}

sink()
close(report)
cat("\n[done] report written to env_report.txt\n")
