#!/usr/bin/env Rscript
# ============================================================================
# sc_06c_finish.R — monocle3 收尾（从已保存的 cds_macro.rds 继续）
# 说明：本机 monocle3 版本的 graph_test 输出列为
#       status, p_value, morans_test_statistic, morans_I, gene_short_name, q_value
#       （新版本叫 moranI，此处做了两种列名兼容）
# 产出：genes_moransI_top100.csv + Fig_mono_F_keygenes.png
# ============================================================================
suppressPackageStartupMessages({
  library(monocle3); library(ggplot2)
})

BASE <- getwd()
OUT  <- file.path(BASE, "scRNA_out", "monocle")
cds  <- readRDS(file.path(OUT, "cds_macro.rds"))

# ---------- 沿轨迹的驱动基因（Moran's I） ----------
pr <- graph_test(cds, neighbor_graph = "principal_graph", cores = 1)
gd <- as.data.frame(pr)
cat("graph_test 输出列:", paste(colnames(gd), collapse = ", "), "\n")

# 列名兼容：新版 moranI / 本机版 morans_I
mi_col <- if ("moranI" %in% colnames(gd)) "moranI" else "morans_I"
q_col  <- if ("moranI_q" %in% colnames(gd)) "moranI_q" else "q_value"
gd$moranI_use <- suppressWarnings(as.numeric(gd[[mi_col]]))
gd <- gd[!is.na(gd$moranI_use), ]
gd$gene <- rownames(gd)
top <- head(gd[order(-gd$moranI_use), ], 100)
write.csv(top[, c("gene", "gene_short_name", mi_col, "p_value", q_col)],
          file.path(OUT, "genes_moransI_top100.csv"), row.names = FALSE)
cat("Top 驱动基因已保存，前 10：\n"); print(head(top[, c("gene", mi_col)], 10))

# ---------- 关键基因沿拟时序表达（铁死亡 + SPP1 相关） ----------
key_g <- intersect(c("SPP1", "SLC7A11", "STEAP3", "HMOX1", "CDKN1A", "FTL", "TFRC"),
                   rownames(cds))
cat("关键基因覆盖:", paste(key_g, collapse = ", "), "\n")
if (length(key_g) >= 3) {
  g6 <- plot_genes_in_pseudotime(cds[key_g, ], color_cells_by = "celltype_detail",
                                 min_expr = 0.5, ncol = 2)
  ggsave(file.path(OUT, "Fig_mono_F_keygenes.png"), g6,
         width = 9, height = 10, dpi = 300)
}

cat("\n=== sc_06c Monocle3 收尾 done ===\n")
