# -*- coding: utf-8 -*-
# 07b_prep_m6a.py: 从 GSE206560_Methylated_RNA_sites.mRNA.xlsx 提取 6 个 IP 样本的 m6A 位点
# sheet: A01-A03 = 非骨化(non-ossified) IP, B01-B03 = 骨化(ossified) IP
# 表头在第 25 行(0基24), 19 列
import openpyxl, csv, os

WD = r"C:/Users/10215/WorkBuddy/2026-10-03-00-10-07"
XLSX = os.path.join(WD, "GSE206560", "GSE206560_Methylated_RNA_sites.mRNA.xlsx")
OUT = os.path.join(WD, "results", "m6a_sites_all.csv")

wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
rows_out = []
summary = []
for sn in wb.sheetnames:
    group = "Ossified" if sn.startswith("B") else "NonOssified"
    ws = wb[sn]
    n = 0
    for i, r in enumerate(ws.iter_rows(values_only=True)):
        if i <= 24 or r[0] is None:
            continue
        # 0chrom 1start 2end 3peakID 4score 5len 6txid 7gene 8m10logp 9fold 10FDR 11GeneID...
        gene = r[7]
        if not gene:
            continue
        rows_out.append([sn, group, r[0], r[1], r[2], r[3], gene,
                         r[6] if r[6] is not None else "",
                         round(r[9], 3) if isinstance(r[9], (int, float)) else "",
                         round(r[10], 6) if isinstance(r[10], (int, float)) else "",
                         r[5] if r[5] is not None else ""])
        n += 1
    summary.append((sn, group, n))

with open(OUT, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["sample", "group", "chrom", "start", "end", "peakID", "gene",
                "transcript_id", "fold_enrichment", "FDR", "peak_length"])
    w.writerows(rows_out)

print("sample,group,n_sites")
for s in summary:
    print(",".join(map(str, s)))
print("total rows:", len(rows_out), "->", OUT)
