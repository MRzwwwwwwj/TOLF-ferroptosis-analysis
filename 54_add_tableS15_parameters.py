# -*- coding: utf-8 -*-
r"""
54_add_tableS15_parameters.py — 回修补充表 Table S15
回应 R1 Comment 1.5：CellChat / monocle3 / single-cell / robustness 的完整参数。
把新 sheet 追加进 supplementary/Table_S10-S14_revision.xlsx，另存为 Table_S10-S15_revision.xlsx，
并同步更新交付目录。
"""
import os
import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Font, Alignment, PatternFill

WS = os.environ.get("TOLF_WS", ".")
SRC = os.path.join(WS, "supplementary", "Table_S10-S14_revision.xlsx")
OUT = os.path.join(WS, "supplementary", "Table_S10-S15_revision.xlsx")

ROWS = [
    ("Section", "Analysis step", "Software / version", "Parameter", "Value / setting"),
    ("Single-cell QC", "Object creation", "Seurat v5", "min.cells / min.features", "3 / 200"),
    ("Single-cell QC", "Cell filtering", "Seurat v5", "nFeature_RNA range", "200-6,000"),
    ("Single-cell QC", "Cell filtering", "Seurat v5", "Mitochondrial fraction", "< 15%"),
    ("Single-cell QC", "Feature space", "Seurat v5", "Shared expressed genes across samples", "14,060"),
    ("Integration", "Batch correction", "Harmony", "Batch variable", "sample"),
    ("Clustering", "PCA", "Seurat v5 (RunPCA)", "npcs", "30"),
    ("Clustering", "Clustering", "Seurat v5 (FindClusters)", "resolution", "0.5 (21 clusters -> 7 cell types)"),
    ("Annotation", "Cell-type labels", "canonical marker panels (cluster-averaged)", "Cell types", "fibroblast, macrophage, endothelial, pericyte/SMC, NK, B, cycling"),
    ("Module scoring", "Ferroptosis score", "Seurat AddModuleScore", "Gene set / random seed", "36-gene ferroptosis set / seed = 42"),
    ("CellChat", "Database", "CellChat v1.6.1", "Ligand-receptor reference", "CellChatDB.human"),
    ("CellChat", "Input", "CellChat v1.6.1", "Expression input / PPI projection", "raw counts; no PPI projection"),
    ("CellChat", "Group labels", "CellChat v1.6.1", "Cell groups", "7 cell types with SPP1+/SPP1- macrophage subdivision"),
    ("CellChat", "Pipeline", "CellChat v1.6.1", "Steps", "subsetData -> identifyOverExpressedGenes -> identifyOverExpressedInteractions -> computeCommunProb -> filterCommunication(min.cells = 10) -> computeCommunProbPathway -> aggregateNet -> netAnalysis_computeCentrality"),
    ("CellChat", "Significance", "CellChat v1.6.1", "Test / aggregation", "default permutation test; probability aggregation by triMean (defaults)"),
    ("CellChat", "First-cohort sizes", "CellChat v1.6.1", "Cells analysed", "21,186 hypertrophic / 12,631 control (GSE294458 + GSE267819)"),
    ("CellChat", "Cross-cohort replication", "CellChat v1.6.1", "Pooling", "normalized counts pooled per group from 9 per-sample GSE271018 objects (6 degenerative / 3 traumatic)"),
    ("Pseudotime", "Preprocessing", "monocle3 (preprocess_cds)", "num_dim", "30"),
    ("Pseudotime", "Embedding", "monocle3 (reduce_dimension)", "preprocess_method", "PCA"),
    ("Pseudotime", "Trajectory", "monocle3 (learn_graph)", "use_partition", "FALSE (single connected graph over all cells)"),
    ("Pseudotime", "Root state", "monocle3 (order_cells)", "root_cells", "cell with the lowest ferroptosis module score (caveat stated in main text)"),
    ("Pseudotime", "Trajectory genes", "monocle3 (graph_test)", "Test", "Moran's I graph-autocorrelation (root-free; independent of rooting)"),
    ("Cell-level tests", "Group comparison", "stats", "Test / correction", "two-sided Wilcoxon rank-sum per cell type; Benjamini-Hochberg"),
    ("Pseudobulk (revision)", "Sample-level validation", "custom (aggregate per sample x cell type)", "Design", "6 degenerative vs 3 traumatic samples (GSE271018)"),
    ("Correlation pipelines", "Three pipelines", "R", "Scoring x method", "GSVA + Pearson; ssGSEA + Spearman; ssGSEA + Pearson (labelled per Methods 4.5; see Table S13B)"),
    ("Nested CV (revision)", "Classifier robustness", "glmnet", "Outer / inner loop", "LOOCV outer; feature selection (limma + cv.glmnet) re-done inside each outer fold"),
    ("Permutation (revision)", "Null distribution", "custom", "Label shuffles (B)", "200 (permutation P = 0.234; Fig. S8A)"),
    ("Scoring sensitivity (revision)", "Ferroptosis variants", "GSVA / ssGSEA", "Gene-set definitions", "full 36-gene; pro-ferroptosis drivers; suppressors; PCD-overlap-removed; iron-metabolism-removed (Table S11, Fig. S8B)"),
    ("m6A depth (revision)", "Depth sensitivity", "custom", "Proxy", "paired Input library read counts; 2.49 -> 2.23 fold after adjustment (Table S13A)"),
    ("Environment", "Computing", "R 4.5.3", "Key packages", "limma, GSVA, edgeR, fgsea, glmnet, pROC, Seurat v5, Harmony, CellChat v1.6.1, monocle3, rms, dcurves, Biostrings, CIBERSORT v1.03"),
]

wb = load_workbook(SRC)
if "TableS15_parameters" in wb.sheetnames:
    del wb["TableS15_parameters"]
ws = wb.create_sheet("TableS15_parameters")

hdr_font = Font(bold=True, color="FFFFFF", size=11)
hdr_fill = PatternFill("solid", fgColor="1F4E79")
sec_font = Font(bold=True, color="1F4E79", size=10)
wrap = Alignment(wrap_text=True, vertical="top")

widths = [20, 26, 26, 22, 60]
for c, wdt in enumerate(widths, 1):
    ws.column_dimensions[chr(64 + c)].width = wdt

for r, row in enumerate(ROWS, 1):
    for c, val in enumerate(row, 1):
        cell = ws.cell(row=r, column=c, value=val)
        cell.alignment = wrap
        if r == 1:
            cell.font = hdr_font
            cell.fill = hdr_fill
        elif c == 1:
            cell.font = sec_font

ws.freeze_panes = "A2"
wb.save(OUT)

xl = pd.ExcelFile(OUT)
print("已生成:", OUT)
print("sheets:", xl.sheet_names)

# 同步交付目录（覆盖旧 S10-S14 文件）
deliv = os.path.join(WS, "回修_v2_2026-11-02", "补充材料_v2", "Table_S10-S15.xlsx")
import shutil
old = os.path.join(WS, "回修_v2_2026-11-02", "补充材料_v2", "Table_S10-S14.xlsx")
shutil.copy(OUT, deliv)
if os.path.exists(old):
    os.remove(old)
print("交付目录已更新 ->", os.path.basename(deliv))
