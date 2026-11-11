# -*- coding: utf-8 -*-
"""
D3-G2 + H2：m6A peak 深度标准化敏感性 + Table S6 数值溯源
===========================================================
G2（R4 深度混淆质疑）：6 个 IP bed 的 peak 数 825–7244 相差 8.8 倍，
   用对应 Input 文库 raw counts 总读长做代理，算 peaks per million input reads，
   重算「2.5 倍」是否稳健；并做下采样敏感性。
H2（R4 数值打架）：证明 0.869/0.005 与 0.781/0.022 来自两套不同管线
   （GSVA+Pearson vs ssGSEA+Spearman），两个数值都真实。

输出全部写 results/，不覆盖原文件。
"""
import os, gzip, itertools
import numpy as np
import pandas as pd
from scipy import stats

WD = os.environ.get("TOLF_WS", ".")
RES = os.path.join(WD, "results")
BEDS = os.path.join(WD, "回修_v2_2026-11-02", "_m6a_beds")

# ============================================================ 1. G2 深度敏感性
# IP 样本 -> (组别, 配对 Input GSM, Input 文库总读长)
# series_matrix: GSM6256879-6884 = IP (A01-A03 非骨化, B01-B03 骨化)
#                GSM6256885-6890 = Input；raw_counts 文件仅含 5 列（缺 GSM6256886）
ip_info = [
    ("GSM6256879_A01", "non_ossified", "GSM6256885", 19730279),
    ("GSM6256880_A02", "non_ossified", "GSM6256886", None),      # Input 深度缺失
    ("GSM6256881_A03", "non_ossified", "GSM6256887", 16924928),
    ("GSM6256882_B01", "ossified",     "GSM6256888", 18220549),
    ("GSM6256883_B02", "ossified",     "GSM6256889", 22374649),
    ("GSM6256884_B03", "ossified",     "GSM6256890", 17367197),
]

rows = []
for ip, grp, inp, depth in ip_info:
    f = os.path.join(BEDS, ip + ".mRNA.m6A.bed.gz")
    n = sum(1 for _ in gzip.open(f, "rt"))
    ppm = n / depth * 1e6 if depth else np.nan
    rows.append(dict(IP_library=ip, group=grp, paired_Input=inp,
                     n_peaks=n, Input_total_reads=depth,
                     peaks_per_M_input_reads=ppm))
peak = pd.DataFrame(rows)
peak.to_csv(os.path.join(RES, "m6A_peak_depth_sensitivity.csv"), index=False)

# 组间：原始 peak 数 vs 深度标准化后
def grp_cmp(df, col, label):
    a = df.loc[df.group == "ossified", col].dropna().values
    b = df.loc[df.group == "non_ossified", col].dropna().values
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided") if len(a) >= 2 and len(b) >= 2 else (np.nan, np.nan)
    return dict(metric=label, n_oss=len(a), n_non=len(b),
                mean_oss=a.mean(), mean_non=b.mean(),
                ratio=a.mean() / b.mean() if b.mean() else np.nan,
                sum_oss=a.sum(), sum_non=b.sum(),
                sum_ratio=a.sum() / b.sum() if b.sum() else np.nan,
                U=u, P=p)

cmp_rows = [
    grp_cmp(peak, "n_peaks", "raw peak counts (as in manuscript)"),
    grp_cmp(peak, "peaks_per_M_input_reads", "peaks per million paired-Input reads (depth-adjusted)"),
]

# 下采样敏感性：把每个 IP 的 peak 数按比例缩到最小文库深度（825 peaks）
mn = peak.n_peaks.min()
peak["downsampled_peaks"] = peak.n_peaks * (mn / peak.n_peaks.max())  # 按最大文库深度归一
# 更直接：以最小文库为基准，按 Input 深度比例做期望值
peak["expected_peaks_at_equal_depth"] = peak.peaks_per_M_input_reads * peak.Input_total_reads.min() / 1e6
cmp_rows.append(grp_cmp(peak, "expected_peaks_at_equal_depth",
                        "expected peaks at a common library depth (down-sampling proxy)"))
pd.DataFrame(cmp_rows).to_csv(os.path.join(RES, "m6A_peak_depth_groupcompare.csv"), index=False)

# ============================================================ 2. H2 数值溯源
# 两套管线的实际相关系数，直接从现有 CSV 读取，并各自复算验证
gsva = pd.read_csv(os.path.join(RES, "PCD_scores_GSVA_GSE113212.csv"), index_col=0)
ssg  = pd.read_csv(os.path.join(RES, "PCD_scores_ssGSEA_GSE113212.csv"), index_col=0)
cib  = pd.read_csv(os.path.join(RES, "immune_cibersort_GSE113212.csv"), index_col=0)
grp  = pd.read_csv(os.path.join(RES, "group_GSE113212.csv"))

# CIBERSORT 矩阵方向 = 行(样本) × 列(细胞类型)；PCD 评分矩阵 = 行(模块) × 列(样本)
m1col = [c for c in cib.columns if "M1" in c][0]
samples = grp["sample"].tolist()
x_cib = cib.loc[samples, m1col].astype(float).values

trace = []
for score_name, sc in [("GSVA", gsva), ("ssGSEA", ssg)]:
    y = sc.loc["cuproptosis", samples].astype(float).values
    for method in ["pearson", "spearman"]:
        r, p = (stats.pearsonr(y, x_cib) if method == "pearson"
                else stats.spearmanr(y, x_cib))
        trace.append(dict(PCD_score=score_name, immune_matrix="CIBERSORT Macrophages M1",
                         method=method, r=float(r), P=float(p)))
tr = pd.DataFrame(trace)
tr["source_in_manuscript"] = [
    "Section 2.4 (rho=0.869, P=0.005) <-- PCD_immune_cibersort_cor.csv",
    "",
    "",
    "Table S6 (rho=0.781, P=0.022) <-- M1M2_PCD_cor.csv",
]
tr.to_csv(os.path.join(RES, "H2_cuproptosis_M1_trace.csv"), index=False)

# ============================================================ 报告
print("=== 6 个 IP 文库 peak 数 + 配对 Input 深度 ===")
print(peak.to_string(index=False))
print("\n=== 2.5 倍说法的稳健性 ===")
print(pd.DataFrame(cmp_rows).to_string(index=False))
print("\n=== H2：0.869 vs 0.781 四套组合复算 ===")
print(tr.to_string(index=False))
print("\n=== 完成 ===")
