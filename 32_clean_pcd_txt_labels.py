# -*- coding: utf-8 -*-
r"""
路径可用环境变量覆盖：TOLF_WS（工作区）/ TOLF_RES（结果 CSV 目录）。
32_clean_pcd_txt_labels.py — 清洗 PCD 模块名 ".txt" 污染 + 修正 R 脚本旧库路径
背景：results/PCD_scores_{GSVA,ssGSEA}_GSE113212.csv 等由旧版 03 脚本生成，
模块名残留 ".txt"（文件名残留），污染 Figure 1B/1D 与 Table S2。
本脚本只改标签字符串，不动任何数值。
"""
import os, re, shutil, csv

WS = os.environ.get("TOLF_WS", r"C:/Users/10215/WorkBuddy/2026-10-03-00-10-07")
RES = os.path.join(WS, "results")
DESK_RES = os.environ.get("TOLF_RES", os.path.join(WS, "results"))
BAK = os.path.join(WS, "_backup", "pcd_txt_cleanup_2026-10-15")
os.makedirs(BAK, exist_ok=True)

# ---------- 1. CSV 清洗（只改 Module/PCD_module 标签列，数值原样保留） ----------
TARGETS = [
    (RES, "PCD_scores_GSVA_GSE113212.csv", "Module"),
    (RES, "PCD_scores_ssGSEA_GSE113212.csv", "Module"),
    (RES, "PCD_scores_groupcompare.csv", "Module"),
    (DESK_RES, "PCD_scores_GSVA_GSE113212.csv", "Module"),
    (DESK_RES, "PCD_scores_ssGSEA_GSE113212.csv", "Module"),
    (DESK_RES, "PCD_scores_groupcompare.csv", "Module"),
]
# PCD_immune_cor.csv 由 R 重跑 05 时重建，这里也顺手清洗旧文件（桌面副本稍后由重建结果覆盖）
TARGETS.append((RES, "PCD_immune_cor.csv", "PCD_module"))
TARGETS.append((DESK_RES, "PCD_immune_cor.csv", "PCD_module"))

n_fixed = 0
for d, fn, col in TARGETS:
    p = os.path.join(d, fn)
    if not os.path.exists(p):
        print("  跳过(不存在):", p)
        continue
    raw = open(p, "rb").read()
    if b".txt" not in raw:
        print("  无需清洗:", fn)
        continue
    # 备份（带来源目录标识）
    bak_name = fn.replace(".csv", "") + "__" + ("desktop" if "Desktop" in d else "ws") + ".csv"
    shutil.copy2(p, os.path.join(BAK, bak_name))
    # 精确替换：仅第一列（模块名列）的 ".txt" 后缀；行首或引号包裹两种情况
    text = raw.decode("utf-8")
    def repl(m):
        return m.group(1)  # 去掉 .txt
    new = re.sub(r'(?m)^("?)(ferroptosis|pyroptosis|cuproptosis|disulfidptosis|panoptosis)\.txt("?)(?=,)',
                 lambda m: m.group(1) + m.group(2) + m.group(3), text)
    # PCD_immune_cor.csv 的 PCD_module 在中间列（"dataset","PCD_module","Immune_cell",...）
    new = re.sub(r'(?<=")(ferroptosis|pyroptosis|cuproptosis|disulfidptosis|panoptosis)\.txt(?=")',
                 lambda m: m.group(1), new)
    open(p, "wb").write(new.encode("utf-8"))
    assert ".txt" not in new or new.count(".txt") == 0, fn
    n_fixed += 1
    print("  已清洗:", p)
print("CSV 清洗完成: %d 个文件, 备份在 %s" % (n_fixed, BAK))

# ---------- 2. R 脚本 .libPaths 旧路径修正 ----------
OLD1 = "D:/A生信/R-4.5.3/library"
NEW1 = "D:/ARLARS/R-4.5.3/library"
OLD2 = r"D:\A生信\R-4.5.3\bin\Rscript.exe"
NEW2 = r"D:\ARLARS\R-4.5.3\bin\Rscript.exe"

rdirs = [os.path.join(WS, "scripts"), os.path.join(WS, "github_zenodo_upload", "03_代码")]
n_path = 0
for rd in rdirs:
    if not os.path.isdir(rd):
        continue
    for fn in os.listdir(rd):
        if not fn.endswith(".R"):
            continue
        p = os.path.join(rd, fn)
        raw = open(p, "rb").read()
        text = raw.decode("utf-8")
        new = text.replace(OLD1, NEW1).replace(OLD2, NEW2)
        if new != text:
            open(p, "wb").write(new.encode("utf-8"))
            n_path += 1
            print("  路径修正:", p)
print("R 脚本路径修正完成: %d 个文件" % n_path)
