"""
生成补充表 S1-S9 汇总 xlsx（2026-10-14）
全部数据来自已核实的结果 CSV，不做任何重新计算（除 S5 的 GSE271018 汇总统计——直接由已核实的 percell 评分表聚合）。
"""
import os
import pandas as pd
from scipy import stats

# 路径可用环境变量覆盖：TOLF_WS（工作区）/ TOLF_RES（结果 CSV 目录）
WD = os.environ.get("TOLF_WS", r"C:/Users/10215/WorkBuddy/2026-10-03-00-10-07")
RES = os.environ.get("TOLF_RES", os.path.join(WD, "results"))
OUT = os.path.join(WD, "supplementary")
os.makedirs(OUT, exist_ok=True)

sheets = {}

# ---------- S1: PCD gene sets ----------
gs_dir = os.path.join(WD, "pcd_genesets")
sets = {}
for f in ["gs_ferroptosis.txt", "gs_pyroptosis.txt", "gs_cuproptosis.txt",
          "gs_disulfidptosis.txt", "gs_panoptosis.txt"]:
    name = f.replace("gs_", "").replace(".txt", "")
    genes = [l.strip() for l in open(os.path.join(gs_dir, f), encoding="utf-8") if l.strip()]
    sets[name] = genes
s1 = pd.DataFrame({k: pd.Series(v) for k, v in sets.items()})
sheets["TableS1_PCD_genesets"] = s1
print("S1:", {k: len(v) for k, v in sets.items()})

# ---------- S2: strong PCD x immune pairs ----------
cor = pd.read_csv(os.path.join(RES, "PCD_immune_cor.csv"))
s2 = cor[(cor["r"].abs() > 0.7) & (cor["P"] < 0.05)].sort_values("r", key=abs, ascending=False)
sheets["TableS2_PCDxImmune_strong"] = s2
print("S2:", len(s2))

# ---------- S3: ossification-hypermethylated PCD genes ----------
s3 = pd.read_csv(os.path.join(RES, "m6a_PCD_hyper_genes_GSE206560.csv"))
sheets["TableS3_m6A_hyper_PCD"] = s3
print("S3:", len(s3), list(s3.columns))

# ---------- S4: m6A regulators ----------
s4 = pd.read_csv(os.path.join(RES, "m6a_regulators_GSE113212.csv"))
sheets["TableS4_m6A_regulators"] = s4
print("S4:", len(s4))

# ---------- S5: single-cell ferroptosis scores (both cohorts) ----------
c1 = pd.read_csv(os.path.join(WD, "scRNA_out", "sc_ferroptosis_by_celltype.csv"))
c1.insert(0, "cohort", "Cohort1_GSE294458+GSE267819")
c1b = pd.read_csv(os.path.join(WD, "scRNA_out", "sc_ferroptosis_LSS_vs_LDH_wilcox.csv"))
c1b.insert(0, "cohort", "Cohort1_GSE294458_LSSvsLDH")

d = pd.read_csv(os.path.join(WD, "GSE271018_results", "GSE271018_percell_scores.csv"))
rows = []
for ct in sorted(d.celltype.unique()):
    dd = d[d.celltype == ct]
    rows.append(dict(cohort="Cohort2_GSE271018", celltype=ct, n=len(dd),
                     mean_score=round(dd.ferro_score.mean(), 4)))
for ct in sorted(d.celltype.unique()):
    dd = d[d.celltype == ct]
    de = dd[dd.group == "Degenerative"].ferro_score
    co = dd[dd.group == "Control"].ferro_score
    pu = stats.mannwhitneyu(de, co).pvalue
    rows.append(dict(cohort="Cohort2_GSE271018_DegVsCtrl", celltype=ct,
                     n_Degenerative=len(de), n_Control=len(co),
                     mean_Degenerative=round(de.mean(), 4), mean_Control=round(co.mean(), 4),
                     delta=round(de.mean() - co.mean(), 4), p_wilcox="%.3g" % pu))
mac = d[d.celltype == "Macrophage"].dropna(subset=["spp1_status"])
for s in ["SPP1+_Macrophage", "SPP1-_Macrophage"]:
    m = mac[mac.spp1_status == s]
    rows.append(dict(cohort="Cohort2_GSE271018_SPP1subset", celltype=s, n=len(m),
                     mean_score=round(m.ferro_score.mean(), 4)))
s5 = pd.DataFrame(rows)
s5_full = pd.concat([
    c1.rename(columns={"celltype": "celltype"}),
    c1b.rename(columns={"celltype": "celltype"}).assign(
        mean_case=lambda x: x["mean_LSS"], mean_ctrl=lambda x: x["mean_LDH"]),
    s5], ignore_index=True, sort=False)
sheets["TableS5_sc_ferroptosis"] = s5_full
print("S5:", len(s5_full))

# ---------- S6: checkpoint + M1/M2 ----------
ck = pd.read_csv(os.path.join(RES, "checkpoint_cor.csv"))
mm = pd.read_csv(os.path.join(RES, "M1M2_PCD_cor.csv"))
s6 = pd.concat([
    ck.assign(analysis="immune_checkpoint").rename(columns={"module": "PCD_module", "checkpoint": "partner"}),
    mm.assign(analysis="M1M2_polarization").rename(columns={"module": "PCD_module", "marker": "partner"})
], ignore_index=True, sort=False)
sheets["TableS6_checkpoint_M1M2"] = s6
print("S6:", len(s6))

# ---------- S7: CellChat LR pairs (4 tables concatenated) ----------
lr = []
for f, lab in [("scRNA_out/cellchat/LR_pairs_all_Control_LDH.csv", "Cohort1_Control"),
               ("scRNA_out/cellchat/LR_pairs_all_LF_hypertrophy.csv", "Cohort1_Hypertrophy"),
               ("GSE271018_cellchat/LR_pairs_all_Control.csv", "Cohort2_Control"),
               ("GSE271018_cellchat/LR_pairs_all_Degenerative.csv", "Cohort2_Degenerative")]:
    x = pd.read_csv(os.path.join(WD, f))
    x.insert(0, "cohort_group", lab)
    lr.append(x)
s7 = pd.concat(lr, ignore_index=True)
sheets["TableS7_CellChat_LR"] = s7
print("S7:", len(s7))

# ---------- S8: Moran's I top 100 ----------
s8 = pd.read_csv(os.path.join(WD, "scRNA_out", "monocle_v2", "genes_moransI_top100.csv"))
sheets["TableS8_pseudotime_genes"] = s8
print("S8:", len(s8))

# ---------- S9: probe re-annotation ----------
s9 = pd.read_csv(os.path.join(RES, "probe2gene_GPL21827_targeted.csv"))
sheets["TableS9_probe2gene"] = s9
print("S9:", len(s9))

# ---------- 写出 xlsx ----------
out = os.path.join(OUT, "Table_S1-S9_supplementary.xlsx")
with pd.ExcelWriter(out, engine="openpyxl") as w:
    for name, df in sheets.items():
        df.to_excel(w, sheet_name=name[:31], index=False)
print("\nSaved:", out, round(os.path.getsize(out) / 1024, 1), "KB")
