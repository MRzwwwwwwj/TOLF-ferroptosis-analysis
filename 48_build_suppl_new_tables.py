# -*- coding: utf-8 -*-
r"""
48_build_suppl_new_tables.py — 2026-11-02 回修
把 D1/D2/D3 新分析整理为新增补充表 Table S10–S14（单独成册，不动原 S1–S9 存档）。
输出 supplementary/Table_S10-S14_revision.xlsx
"""
import os
import pandas as pd

WS = os.environ.get("TOLF_WS", ".")
RES = os.path.join(WS, "results")
OUT = os.path.join(WS, "supplementary", "Table_S10-S14_revision.xlsx")

def rd(f):
    return pd.read_csv(os.path.join(RES, f))

def block(writer, sheet, title, df, row):
    """在同一 sheet 写一个带标题的数据块，返回下一可用行。"""
    pd.DataFrame({sheet: [title]}).to_excel(writer, sheet_name=sheet,
                                            index=False, header=False, startrow=row)
    row += 1
    df.to_excel(writer, sheet_name=sheet, index=False, startrow=row)
    return row + len(df) + 2  # 空一行

with pd.ExcelWriter(OUT, engine="openpyxl") as w:
    # ---- Table S10: signature robustness ----
    s10a = rd("nested_CV_summary.csv")
    s10a.columns = ["Analysis", "Value"]
    s10b = rd("nested_CV_gene_stability.csv")
    s10b.columns = ["Gene", "Fold selection frequency (of 8 outer folds)"]
    s10c = rd("lambda_grid_extension.csv")
    r = 0
    r = block(w, "TableS10_signature_robustness",
              "A. Nested cross-validation and permutation test (six-gene signature)", s10a, r)
    r = block(w, "TableS10_signature_robustness",
              "B. Gene selection stability across nested-CV outer folds", s10b, r)
    r = block(w, "TableS10_signature_robustness",
              "C. LASSO lambda-grid extension (edge behaviour)", s10c, r)

    # ---- Table S11: ferroptosis scoring sensitivity ----
    s11 = rd("ferroptosis_sensitivity_immune_cor.csv")
    s11m = s11[s11["immune_cell"] == "Macrophages"].copy()
    s11m = s11m[["set", "method", "immune_cell", "r", "P"]]
    s11m.columns = ["Gene set", "Scoring method", "Immune cell", "Pearson r", "P value"]
    s11g = rd("ferroptosis_sensitivity_groupcompare.csv")
    s11g.columns = ["Gene set", "Method", "Module", "Mean young", "Mean elderly",
                    "Delta", "Direction", "Wilcoxon P"]
    r = 0
    r = block(w, "TableS11_scoring_sensitivity",
              "A. Ferroptosis-score vs macrophage infiltration across 5 gene-set definitions", s11m, r)
    r = block(w, "TableS11_scoring_sensitivity",
              "B. Group comparison (hypertrophy vs control) under each scoring variant", s11g, r)

    # ---- Table S12: pseudobulk + BH ----
    s12p = rd("pseudobulk_GSE271018_wilcox.csv")
    s12b = rd("BH_correction_summary.csv")
    r = 0
    r = block(w, "TableS12_pseudobulk_BH",
              "A. Sample-level (pseudobulk) ferroptosis comparison in GSE271018", s12p, r)
    r = block(w, "TableS12_pseudobulk_BH",
              "B. Benjamini-Hochberg correction across all correlation matrices", s12b, r)

    # ---- Table S13: m6A depth + H2 provenance ----
    s13d = rd("m6A_peak_depth_groupcompare.csv")
    s13h = rd("H2_cuproptosis_M1_trace.csv")
    r = 0
    r = block(w, "TableS13_m6A_depth_provenance",
              "A. m6A peak-count sensitivity to sequencing-depth adjustment (GSE206560)", s13d, r)
    r = block(w, "TableS13_m6A_depth_provenance",
              "B. Provenance of the cuproptosis-Macrophage correlation (three pipelines)", s13h, r)

    # ---- Table S14: SLC7A11 cell-type expression ----
    s14 = rd("SLC7A11_keygenes_by_celltype.csv")
    s14.columns = ["Gene", "Cell type", "Mean expression", "% cells expressing", "N cells"]
    block(w, "TableS14_SLC7A11_celltype",
          "Key-gene expression by cell type (single-cell, GSE294458)", s14, 0)

print("已生成:", OUT, "%.1f KB" % (os.path.getsize(OUT) / 1024))
# 校验读回
xl = pd.ExcelFile(OUT)
print("sheets:", xl.sheet_names)
