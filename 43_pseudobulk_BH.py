# -*- coding: utf-8 -*-
"""
D2-F：GSE271018 样本级 pseudobulk 验证 + 相关分析 BH 校正
=========================================================
回应 R4「细胞级 P 值当独立样本（伪重复）」与「多重比较未校正」。

输出（全部写 results/，不覆盖任何原文件）：
  - pseudobulk_GSE271018_by_celltype.csv   每样本×每细胞类型 mean ferro_score
  - pseudobulk_GSE271018_wilcox.csv        样本级 6v3 Wilcoxon + BH
  - pseudobulk_GSE271018_SPP1macro.csv     样本级 SPP1+/SPP1- 巨噬细胞
  - BH_correction_summary.csv              三套相关矩阵 + 合并 BH 存活情况
  - BH_corrected_PCD_immune_cor.csv        带 q 值的 115 对
  - BH_corrected_checkpoint_cor.csv        带 q 值的 115 对
  - BH_corrected_M1M2_PCD_cor.csv          带 q 值的 30 对
"""
import os, itertools
import numpy as np
import pandas as pd
from scipy import stats

WD = os.environ.get("TOLF_WS", ".")
RES = os.path.join(WD, "results")


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


# ---------------------------------------------------------------- 1. pseudobulk
pc = pd.read_csv(os.path.join(WD, "GSE271018_results", "GSE271018_percell_scores.csv"))
# 样本级聚合：每样本 × 每细胞类型 的 mean ferro_score（= pseudobulk 的评分版本）
pb = (pc.groupby(["sample", "group", "celltype"], as_index=False)
        .agg(mean_ferro=("ferro_score", "mean"),
             median_ferro=("ferro_score", "median"),
             n_cells=("ferro_score", "size")))
pb.to_csv(os.path.join(RES, "pseudobulk_GSE271018_by_celltype.csv"), index=False)

rows = []
for ct, d in pb.groupby("celltype"):
    a = d.loc[d.group == "Degenerative", "mean_ferro"].values
    b = d.loc[d.group == "Control", "mean_ferro"].values
    if len(a) >= 2 and len(b) >= 2:
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        rows.append(dict(celltype=ct, n_deg=len(a), n_ctrl=len(b),
                         mean_deg=a.mean(), mean_ctrl=b.mean(),
                         delta=a.mean() - b.mean(),
                         median_deg=np.median(a), median_ctrl=np.median(b),
                         U=u, P_sample_level=p))
pbw = pd.DataFrame(rows)
pbw["q_BH"] = bh(pbw["P_sample_level"].values)
pbw["sig_q05"] = pbw["q_BH"] < 0.05
pbw = pbw.sort_values("P_sample_level")
pbw.to_csv(os.path.join(RES, "pseudobulk_GSE271018_wilcox.csv"), index=False)

# SPP1+/- 巨噬细胞：样本级
mac = pc[pc.celltype == "Macrophage"].copy()
sp = (mac.dropna(subset=["spp1_status"])
        .groupby(["sample", "group", "spp1_status"], as_index=False)
        .agg(mean_ferro=("ferro_score", "mean"), n_cells=("ferro_score", "size")))
sp.to_csv(os.path.join(RES, "pseudobulk_GSE271018_SPP1macro.csv"), index=False)

sp_rows = []
for st, d in sp.groupby("spp1_status"):
    a = d.loc[d.group == "Degenerative", "mean_ferro"].values
    b = d.loc[d.group == "Control", "mean_ferro"].values
    if len(a) >= 2 and len(b) >= 2:
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        sp_rows.append(dict(spp1_status=st, n_deg=len(a), n_ctrl=len(b),
                            mean_deg=a.mean(), mean_ctrl=b.mean(),
                            delta=a.mean() - b.mean(), P_sample_level=p))
spw = pd.DataFrame(sp_rows)
spw.to_csv(os.path.join(RES, "pseudobulk_GSE271018_SPP1macro_wilcox.csv"), index=False)

# ---------------------------------------------------------------- 2. BH 校正
specs = [
    ("PCD_immune_cor.csv", "P", "r", "PCD_module", "Immune_cell", "PCD x immune (ssGSEA28, 115 pairs)"),
    ("checkpoint_cor.csv", "P", "rho", "module", "checkpoint", "PCD x checkpoint (115 pairs)"),
    ("M1M2_PCD_cor.csv", "P", "rho", "module", "marker", "PCD x M1M2 (30 pairs)"),
]
summ, allp = [], []
for fn, pcol, rcol, c1, c2, label in specs:
    d = pd.read_csv(os.path.join(RES, fn))
    d["q_BH"] = bh(d[pcol].values)
    d["sig_q05"] = d["q_BH"] < 0.05
    d["nominal_p05"] = d[pcol] < 0.05
    out = os.path.join(RES, "BH_corrected_" + fn)
    d.to_csv(out, index=False)
    summ.append(dict(matrix=label, n=len(d),
                     nominal_p05=int(d["nominal_p05"].sum()),
                     BH_q05=int(d["sig_q05"].sum()),
                     min_p=round(float(d[pcol].min()), 6),
                     min_q=round(float(d["q_BH"].min()), 6)))
    allp.extend(d[pcol].tolist())

# 最保守：三套合并（260 对）统一 BH
summ.append(dict(matrix="ALL COMBINED (most conservative)", n=len(allp),
                 nominal_p05=int(sum(np.array(allp) < 0.05)),
                 BH_q05=int((bh(np.array(allp)) < 0.05).sum()),
                 min_p=round(float(np.min(allp)), 6),
                 min_q=round(float(np.min(bh(np.array(allp)))), 6)))
pd.DataFrame(summ).to_csv(os.path.join(RES, "BH_correction_summary.csv"), index=False)

# ---------------------------------------------------------------- 报告
print("=== pseudobulk (样本级 6 Degenerative vs 3 Control) ===")
print(pbw.to_string(index=False))
print("\n=== SPP1+/- 巨噬细胞（样本级）===")
print(spw.to_string(index=False))
print("\n=== BH 校正汇总 ===")
print(pd.DataFrame(summ).to_string(index=False))
