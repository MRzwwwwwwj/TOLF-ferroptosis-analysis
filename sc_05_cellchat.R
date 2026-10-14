#!/usr/bin/env Rscript
# ============================================================================
# sc_05_cellchat.R — 细胞间通讯分析（第一队列 GSE294458/267819, n=33,817）
# 目标：比较 肥厚组(LF_hypertrophy) vs 对照组 的细胞通讯强度差异，
#       聚焦 SPP1+ 巨噬细胞 -> 成纤维细胞/成软骨细胞 的 SPP1 信号轴
# 输入：scRNA_out/sc_annotated.rds
# 输出：scRNA_out/cellchat/ 目录（rds 对象 + 图 + L-R 对表格）
#
# ★ 依赖：CellChat 包（需手动安装，见 加分项分析操作指南.md）
# ★ 人工复核点见脚本内 [复核] 标记
# ============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(CellChat); library(ggplot2); library(patchwork)
})

BASE <- getwd()
OUT  <- file.path(BASE, "scRNA_out", "cellchat")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(file.path(BASE, "scRNA_out", "sc_annotated.rds"))
DefaultAssay(obj) <- "RNA"
try(obj <- JoinLayers(obj), silent = TRUE)   # Seurat v5 兼容，v4 无此函数

data.input <- tryCatch(GetAssayData(obj, layer = "data"),
                       error = function(e) GetAssayData(obj, slot = "data"))
meta <- obj@meta.data

# ---------- 细胞标签：用 celltype_detail（含 SPP1+/SPP1- 巨噬细胞区分） ----------
labels <- as.character(meta$celltype_detail)
names(labels) <- rownames(meta)
groups <- as.character(meta$group)
cat("=== celltype_detail 构成 ===\n"); print(table(labels))
cat("=== group 构成 ===\n"); print(table(groups))
# [复核] 确认 SPP1+_Macrophage 有足够细胞（几百以上）参与通讯推断

glev <- names(which(table(groups) > 100))   # 每组至少 100 细胞才做通讯
stopifnot(length(glev) == 2)

# ---------- 每组建 CellChat 对象并跑标准流程（支持断点续跑） ----------
cc_list <- list()
for (g in glev) {
  cells <- rownames(meta)[groups == g]
  rds_g <- file.path(OUT, sprintf("cellchat_%s.rds", g))
  if (file.exists(rds_g)) {
    cat(sprintf("\n>>> [skip] %s 已完成，读取 %s\n", g, rds_g))
    cc_list[[g]] <- readRDS(rds_g); next
  }
  cat(sprintf("\n>>> CellChat: %s (n=%d cells)\n", g, length(cells)))
  cc <- createCellChat(object = as(data.input[, cells], "dgCMatrix"),
                       meta = data.frame(labels = labels[cells], row.names = cells),
                       group.by = "labels")
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  # projectData 为可选的 PPI 平滑步骤（官方教程默认跳过，computeCommunProb 默认用原始数据）
  cc <- computeCommunProb(cc)                 # ~几分钟
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)                      # CellChat 1.6.1 中整合网络（旧名 computeNetwork 不存在）
  cc <- netAnalysis_computeCentrality(cc)     # 计算网络中心性（角色热图需要）
  saveRDS(cc, rds_g)
  cc_list[[g]] <- cc
}

# ---------- 合并比较 ----------
cc_comp <- mergeCellChat(cc_list, add.names = glev)

# 图1：总交互数/强度比较
p1 <- compareInteractions(cc_comp, show.legend = FALSE) + ggtitle("Interaction counts")
p2 <- compareInteractions(cc_comp, show.legend = FALSE, measure = "weight") + ggtitle("Interaction strength")
ggsave(file.path(OUT, "Fig_cc_A_compare_interactions.png"), p1 + p2,
       width = 8, height = 4.5, dpi = 300)

# 图2：差异交互气泡图（疾病组 vs 对照组）
gg <- netVisual_bubble(cc_comp, comparison = c(1, 2),
                       angle.x = 45, return.data = FALSE)
ggsave(file.path(OUT, "Fig_cc_B_diff_bubble.png"), gg,
       width = 12, height = 20, dpi = 300, limitsize = FALSE)

# 图3：信号通路总览排序（看 SPP1 排位）
p3 <- rankNet(cc_comp, mode = "comparison", stacked = FALSE, do.stat = TRUE)
ggsave(file.path(OUT, "Fig_cc_C_rankNet.png"), p3, width = 10, height = 14, dpi = 300)

# 图4：细胞角色热图（谁发出/接收信号）——两组通路数不同，分别保存
p4a <- netAnalysis_signalingRole_heatmap(cc_list[[glev[1]]], pattern = "outgoing", width = 8, height = 8, title = glev[1])
p4b <- netAnalysis_signalingRole_heatmap(cc_list[[glev[2]]], pattern = "outgoing", width = 8, height = 8, title = glev[2])
pdf(file.path(OUT, "Fig_cc_D_signalingRole_heatmap.pdf"), width = 8, height = 9)
print(p4a)
grid::grid.newpage()
print(p4b)
dev.off()

# ---------- 聚焦 SPP1 通路 ----------
sig_all <- sapply(cc_list, function(x) x@netP$pathways)
common_sig <- intersect(sig_all[[1]], sig_all[[2]])
cat("\nCommon pathways:", paste(common_sig, collapse = ", "), "\n")

if ("SPP1" %in% sig_all[[1]] | "SPP1" %in% sig_all[[2]]) {
  # 图5：SPP1 通路网络图（两组各一张）
  pdf(file.path(OUT, "Fig_cc_E_SPP1_network.pdf"), width = 12, height = 6)
  par(mfrow = c(1, 2), xpd = TRUE)
  for (g in glev) {
    cc <- cc_list[[g]]
    if ("SPP1" %in% cc@netP$pathways) {
      netVisual_aggregate(cc, signaling = "SPP1",
                          vertex.receiver = c(1, 4), layout = "circle",
                          edge.weight.max = NULL)
      title(g, line = -1)
    }
  }
  dev.off()
  # 图6：SPP1 通路贡献的 L-R 对
  ggsave(file.path(OUT, "Fig_cc_F_SPP1_LR_pairs.png"),
         netVisual_bubble(cc_list[[glev[1]]], sources.use = "SPP1+_Macrophage",
                          signaling = "SPP1", remove.isolate = FALSE, return.data = FALSE),
         width = 6, height = 5, dpi = 300)
  # [复核] 重点看 SPP1+_Macrophage -> Fibroblast 的 SPP1-(ITGA/CD44) 对是否显著
}

# ---------- 导出 L-R 对明细表 ----------
for (g in glev) {
  df <- subsetCommunication(cc_list[[g]])
  df <- as.data.frame(df)
  write.csv(df, file.path(OUT, sprintf("LR_pairs_all_%s.csv", g)), row.names = FALSE)
  key <- subset(df, source %in% c("SPP1+_Macrophage", "SPP1-_Macrophage") &
                  grepl("Fibroblast|Pericyte|Endothelial", target))
  write.csv(key, file.path(OUT, sprintf("LR_pairs_macro_to_stromal_%s.csv", g)),
            row.names = FALSE)
  cat(sprintf("\n%s: SPP1+ 巨噬细胞发出的显著 L-R 对 %d 条\n",
              g, sum(df$source == "SPP1+_Macrophage")))
}

cat("\n=== sc_05 CellChat done ===\n")
cat("输出目录:", OUT, "\n")
