#!/usr/bin/env Rscript
# ============================================================================
# 02_cellchat.R — GSE271018 第二队列 CellChat 细胞通讯分析（服务器版）
# 设计：逐样本读 RDS -> 注释(原注释列优先,标志基因兜底) -> 巨噬细胞 SPP1+/-
#       -> 按组拼接稀疏矩阵(Degenerative: N0-N10 / Control: C1,C3,C5)
#       -> 每组 CellChat -> mergeCellChat 比较 -> 图 + 表 + rds
# 用法：cd ~ && nohup Rscript server_cellchat_GSE271018/02_cellchat.R > run_cellchat.log 2>&1 &
# 预计耗时：30-90 分钟（51,686 细胞）
# ============================================================================
.libPaths(c(path.expand("~/Rlibs_new"), .libPaths()))   # CellChat 装在个人库

suppressMessages({ library(Seurat); library(CellChat); library(ggplot2) })

DATA_DIR <- path.expand("~/GSE271018")          # 绝对路径，任意目录下运行均可
OUT      <- path.expand("~/rerun_2026-10-13_cellchat")
dir.create(OUT, showWarnings = FALSE)

infer_group <- function(fname) {
  if (grepl("Tendon_C", fname, ignore.case = TRUE)) return("Control")
  if (grepl("Tendon_N", fname, ignore.case = TRUE)) return("Degenerative")
  return(NA)
}
get_counts <- function(obj)
  tryCatch(GetAssayData(obj, layer = "counts"), error = function(e) GetAssayData(obj, slot = "counts"))
get_data <- function(obj)
  tryCatch(GetAssayData(obj, layer = "data"), error = function(e) GetAssayData(obj, slot = "data"))

marker_panel <- list(
  Fibroblast = c("COL1A1","COL1A2","DCN","LUM","COL3A1"),
  Macrophage = c("CD68","CD163","LYZ","C1QA","C1QB","CD14"),
  Endothelial = c("PECAM1","VWF","KDR","CLDN5"),
  Pericyte_SMC = c("ACTA2","RGS5","MYH11","TAGLN"),
  T_NK = c("CD3D","NKG7","GNLY","IL7R"),
  B_cell = c("CD79A","MS4A1","CD19"),
  Mast_cell = c("TPSAB1","TPSB2","KIT","CPA3"),
  Chondrocyte = c("ACAN","SOX9","COL2A1","CRTAC1")
)
pick_celltype_col <- function(md) {
  cand <- grep("cell.?type|annotation|celltype|sub.?type|major|label",
               colnames(md), ignore.case = TRUE, value = TRUE)
  cand <- setdiff(cand, c("nCount_RNA","nFeature_RNA"))
  for (cn in cand) {
    lv <- unique(as.character(md[[cn]]))
    if (length(lv) >= 2 && length(lv) <= 40 &&
        any(grepl("macro|fibro", lv, ignore.case = TRUE))) return(cn)
  }
  for (cn in colnames(md)) {
    v <- md[[cn]]
    if (is.character(v) || is.factor(v)) {
      lv <- unique(as.character(v))
      if (length(lv) >= 2 && length(lv) <= 40 &&
          any(grepl("macro|fibro", lv, ignore.case = TRUE))) return(cn)
    }
  }
  NA
}
annotate_by_markers <- function(obj) {
  message("    marker-panel fallback")
  obj <- tryCatch(NormalizeData(obj, verbose = FALSE), error = function(e) obj)
  obj <- tryCatch({
    obj <- FindVariableFeatures(obj, verbose = FALSE)
    obj <- ScaleData(obj, verbose = FALSE)
    obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
    obj <- FindNeighbors(obj, dims = 1:20, verbose = FALSE)
    obj <- FindClusters(obj, resolution = 0.6, verbose = FALSE)
    obj
  }, error = function(e) obj)
  avg <- tryCatch(AverageExpression(obj, assays = "RNA")$RNA,
                  error = function(e) NULL)
  if (is.null(avg)) avg <- AggregateExpression(obj, assays = "RNA")$RNA
  cl_clean <- sub("^g", "", colnames(avg))
  score <- sapply(marker_panel, function(g) {
    g <- intersect(g, rownames(avg))
    if (length(g) == 0) return(rep(0, length(cl_clean)))
    colMeans(as.matrix(avg[g, , drop = FALSE]))
  })
  best <- apply(score, 1, function(x)
    if (max(x) == 0) "Unknown" else colnames(score)[which.max(x)])
  names(best) <- cl_clean
  as.character(unname(best[as.character(obj$seurat_clusters)]))
}

# ---------- 阶段1：逐样本提取 标签 + 归一化矩阵，按组累积 ----------
group_mats <- list(); group_labels <- list()
rds_files <- sort(list.files(DATA_DIR, pattern = "\\.rds(\\.gz)?$",
                             full.names = TRUE, ignore.case = TRUE))
stopifnot(length(rds_files) > 0)
cat("=== sample -> group ===\n")
for (f in rds_files) cat(sprintf("  %-40s -> %s\n", basename(f), infer_group(basename(f))))
cat("\n")

for (f in rds_files) {
  sid <- sub("\\.rds(\\.gz)?$", "", basename(f), ignore.case = TRUE)
  grp <- infer_group(basename(f))
  cat("=====", sid, "(", grp, ") =====\n")
  obj <- tryCatch(readRDS(f), error = function(e) e)
  if (inherits(obj, "error") || !inherits(obj, "Seurat")) { cat("  !! skip\n"); next }
  cat("  cells:", ncol(obj), "\n")

  try(obj@graphs <- list(), silent = TRUE)   # v4 老对象 SNN 兼容
  ok <- tryCatch({ validObject(obj); TRUE }, error = function(e) FALSE)
  if (!ok) {
    md0 <- obj@meta.data; cm <- get_counts(obj)
    obj <- CreateSeuratObject(counts = cm, meta.data = md0)
  }
  obj <- tryCatch(NormalizeData(obj, verbose = FALSE), error = function(e) obj)

  md <- obj@meta.data
  ct_col <- pick_celltype_col(md)
  ct <- if (!is.na(ct_col)) as.character(md[[ct_col]]) else annotate_by_markers(obj)
  ct[is.na(ct)] <- "Unknown"
  if (!is.na(ct_col)) cat("  annotation col:", ct_col, "\n")

  # 巨噬细胞 SPP1+/-
  dat <- get_data(obj)
  if ("SPP1" %in% rownames(dat)) {
    sp <- as.numeric(dat["SPP1", ])
    ct <- ifelse(grepl("macro|mono", ct, ignore.case = TRUE),
                 ifelse(sp > 0, "SPP1+_Macrophage", "SPP1-_Macrophage"), ct)
  }
  keep <- ct != "Unknown"
  cells <- colnames(obj)[keep]; ct <- ct[keep]
  cat("  kept:", sum(keep), "/", length(keep), "cells\n")
  cat("  celltypes:", paste(names(table(ct)), collapse = ", "), "\n")

  m <- as(get_data(obj)[, cells], "dgCMatrix")
  group_mats[[grp]]  <- c(group_mats[[grp]],  list(m))
  group_labels[[grp]] <- c(group_labels[[grp]], list(setNames(ct, cells)))
  rm(obj, dat, m); gc(verbose = FALSE)
}

# 基因并集对齐（不同样本基因集略有差异 -> 缺失行补 0）
merge_group <- function(mats, labs) {
  genes <- sort(unique(unlist(lapply(mats, rownames))))
  # 注意：必须带 x = numeric(0) 显式创建 dgCMatrix（数值型）；
  # 缺省时 sparseMatrix 返回 ngCMatrix（逻辑型），赋值会把表达量强制转 0/1！
  mat <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                              dims = c(length(genes), sum(lengths(labs))),
                              dimnames = list(genes, unlist(lapply(labs, names))))
  off <- 0
  for (k in seq_along(mats)) {
    mk <- mats[[k]]; idx <- match(rownames(mk), genes)
    mat[idx, off + seq_len(ncol(mk))] <- mk
    off <- off + ncol(mk)
  }
  list(mat = mat, labels = unlist(labs))
}

# ---------- 阶段2：每组 CellChat（支持断点续跑：RDS 已存在则跳过） ----------
glev <- names(group_mats)
cc_list <- list()
for (g in glev) {
  rds_g <- file.path(OUT, sprintf("cellchat_%s.rds", g))
  if (file.exists(rds_g)) {
    cat("\n>>> [skip]", g, "已完成，读取", rds_g, "\n")
    cc_list[[g]] <- readRDS(rds_g); next
  }
  cat("\n>>> CellChat:", g, ":", sum(lengths(group_labels[[g]])), "cells\n")
  mg <- merge_group(group_mats[[g]], group_labels[[g]])
  cc <- createCellChat(object = mg$mat,
                       meta = data.frame(labels = as.character(mg$labels),
                                         row.names = colnames(mg$mat)),
                       group.by = "labels")
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  # projectData 为可选 PPI 平滑步骤，CellChat 官方教程默认跳过（computeCommunProb 默认用原始数据）
  cc <- computeCommunProb(cc)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)                      # CellChat 1.6.x 中整合网络（computeNetwork 不存在）
  cc <- netAnalysis_computeCentrality(cc)
  saveRDS(cc, file.path(OUT, sprintf("cellchat_%s.rds", g)))
  cc_list[[g]] <- cc
  df <- as.data.frame(subsetCommunication(cc))
  write.csv(df, file.path(OUT, sprintf("LR_pairs_all_%s.csv", g)), row.names = FALSE)
  cat(g, "显著 L-R 对:", nrow(df), "条\n")
  rm(mg); gc(verbose = FALSE)
}

# ---------- 阶段3：合并比较 + 图 ----------
cc_comp <- mergeCellChat(cc_list, add.names = glev)

p1 <- compareInteractions(cc_comp, show.legend = FALSE) + ggtitle("Interaction counts")
p2 <- compareInteractions(cc_comp, show.legend = FALSE, measure = "weight") + ggtitle("Interaction strength")
ggsave(file.path(OUT, "Fig_cc_A_compare_interactions.png"), p1 + p2,
       width = 8, height = 4.5, dpi = 300)

gg <- netVisual_bubble(cc_comp, comparison = c(1, 2), angle.x = 45, return.data = FALSE)
ggsave(file.path(OUT, "Fig_cc_B_diff_bubble.png"), gg,
       width = 12, height = 20, dpi = 300, limitsize = FALSE)

p3 <- rankNet(cc_comp, mode = "comparison", stacked = FALSE, do.stat = TRUE)
ggsave(file.path(OUT, "Fig_cc_C_rankNet.png"), p3, width = 10, height = 14, dpi = 300)

# SPP1 通路网络图
pdf(file.path(OUT, "Fig_cc_E_SPP1_network.pdf"), width = 12, height = 6)
par(mfrow = c(1, 2), xpd = TRUE)
for (g in glev) {
  cc <- cc_list[[g]]
  if ("SPP1" %in% cc@netP$pathways) {
    netVisual_aggregate(cc, signaling = "SPP1", layout = "circle", edge.weight.max = NULL)
    title(g, line = -1)
  }
}
dev.off()

cat("\n=== 02_cellchat done ===\n")
cat("取回本地图（本地 PowerShell）:\n")
cat('scp -P 25100 -r weijian@yxszxxsx.com:~/GSE271018_cellchat .\n')
