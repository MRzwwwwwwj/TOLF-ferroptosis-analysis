# GSE271018 第二重单细胞验证 —— 服务器操作指令

> 目标：在 Linux 服务器上完成 GSE271018（9 个样本：3 对照 + 6 退变，6.3GB RDS）的铁死亡验证分析。
> 约束：R 包为老版本、**不可安装/更新任何包**。本套脚本已按此设计：
> - 全程只用 `Seurat + ggplot2 + Matrix`（几乎任何老服务器都有）
> - 自动适配 Seurat v3 / v4 / v5（脚本内部版本分支，无需人工改代码）
> - 逐样本处理、不合并大对象（内存友好 + 绕开老版本 merge/整合 API 差异）
> - 不需要 harmony（不合并就不需要批次校正）、不需要联网装包

---

## 第 0 步：上传本文件夹到服务器

在**本地**（Windows PowerShell 或 Git Bash）执行：

```bash
# 把 server_GSE271018 整个文件夹上传到服务器家目录（替换成你的账号和地址）
scp -r "C:\Users\10215\WorkBuddy\2026-10-03-00-10-07\server_GSE271018" 用户名@服务器IP:~/
```

文件夹内容（共 5 个文件，几十 KB）：

| 文件 | 作用 |
|---|---|
| `00_check_env.R` | 环境验证（R/包版本/内存/磁盘） |
| `01_download.sh` | 下载 6.3GB 数据 + 解包 + 显示分组 |
| `02_probe_rds.R` | 探查 9 个 RDS 的结构和注释列 |
| `03_ferroptosis_validate.R` | 主分析（打分 + 统计 + 图） |
| `gs_ferroptosis.txt` | 铁死亡 36 基因集（与论文其余部分同源） |

---

## 第 1 步：环境验证（必须先做，不可跳过）

```bash
ssh 用户名@服务器IP
cd ~/server_GSE271018
Rscript 00_check_env.R
cat env_report.txt
```

**判断标准**（报告里也有自动结论）：

| 检查项 | 通过标准 | 不通过怎么办 |
|---|---|---|
| R 版本 | ≥ 3.6 | 发回报告，换方案 |
| Seurat | 已安装，≥ 3.0 | 未装 → 发回，我走纯 Matrix 降级管线；< 3.0 → 发回 |
| ggplot2 | 已安装 | 发回（图可改 base R） |
| 内存 | ≥ 32 GB 总内存 | < 16 GB 发回（C3 文件加载需 8-15GB） |
| 磁盘 | ≥ 20 GB 空闲 | 清理后再来 |

✅ 报告结尾显示 `ENV OK` → 进第 2 步
⚠️ 显示 `ENV INCOMPLETE` 或任何拿不准的 → **把 `env_report.txt` 内容原样发给我**，我给降级方案，不要硬跑

---

## 第 2 步：下载数据（约 6.3GB，放后台跑）

```bash
nohup bash 01_download.sh > download.log 2>&1 &
# 查看进度：
tail -f download.log
```

下完（log 里出现"下载完成"）后**核对分组**：

```bash
# log 最后会打印每个样本的 source_name，确认：
#   Tendon_C1/C3/C5 → 对照（control / non-degenerative）
#   Tendon_N0/N1/N2/N3/N5/N10 → 退变（degenerative ligament）
```

❗ 如果实际分组和上面不一致（比如某个 C 其实是病变），**先告诉我再往下跑**——分组错则全错。

---

## 第 3 步：探查 RDS 结构（几分钟）

```bash
Rscript 02_probe_rds.R
cat probe_report.txt
```

**把 `probe_report.txt` 内容发给我**（这一步很重要）。我看三样东西：

1. **9 个文件是否都能被老 Seurat 读入**——若某文件报 READ FAILED（通常是 v5 对象撞上老 Seurat），我给你针对该文件的降级处理
2. **有没有现成的细胞注释列**（原文已发表，大概率有 `cell_type` / `seurat_clusters` 之类）——有的话我会让你把列名填进 03 脚本第 25 行的 `CELLTYPE_COL_OVERRIDE`，**直接用原文注释**，比我们自己注释更权威
3. 细胞总数和 meta.data 列名，确认分析设计

如果你想直接跑不等我确认：03 脚本的自动识别逻辑会自己找注释列（找含 macro/fibro 关键词的低基数列），找不到才走标志基因面板兜底。直接跑也可以，**但发我确认一遍更稳**。

---

## 第 4 步：主分析（约 0.5–2 小时，放后台）

```bash
# （可选）如果第 3 步确认了注释列名，先编辑这一行：
#   CELLTYPE_COL_OVERRIDE <- "实际的列名"
# 用 vi/nano 改 03_ferroptosis_validate.R 第 25 行

nohup Rscript 03_ferroptosis_validate.R > log03.txt 2>&1 &
tail -f log03.txt
```

脚本设计是**断点续跑**的：每个样本的打分结果单独存 `GSE271018_scores/<样本>.csv`，已存在的样本会跳过——中途被杀（内存不足等）直接重新执行同一命令即可。

出现 `[ALL DONE]` 即完成。产物在 `GSE271018_results/`：

| 文件 | 内容 |
|---|---|
| `GSE271018_percell_scores.csv` | 每细胞一行（细胞类型/分组/铁死亡评分/SPP1 状态/关键基因） |
| `GSE271018_ferroptosis_by_celltype.csv` | 铁死亡评分细胞类型排序（**验证点①：巨噬细胞是否第 1**） |
| `GSE271018_deg_vs_ctrl_wilcox.csv` | 退变 vs 对照逐细胞类型统计（**验证点②：巨噬细胞增幅**） |
| `GSE271018_SPP1macro_ferroptosis.csv` | SPP1+ vs SPP1- 巨噬细胞（**验证点③：与原文 SPP1+ 亚群对接**） |
| `FigS5A–D` 四张 PDF | 细胞类型排序图、分组对比图、SPP1 对比图、关键基因热图 |

---

## 第 5 步：结果传回本地

```bash
# 服务器上先打包
cd ~/server_GSE271018
tar -czf GSE271018_results.tar.gz GSE271018_results/ env_report.txt probe_report.txt log03.txt

# 本地 PowerShell/Git Bash 拉回
scp 用户名@服务器IP:~/server_GSE271018/GSE271018_results.tar.gz "C:\Users\10215\WorkBuddy\2026-10-03-00-10-07\"
```

拉回来后**发我一声**，我负责：解包核对 → 与第一重验证（GSE294458/267819）对照写论文段落 → 更新 manuscript v3。

---

## 预期结果（对照第一重验证）

| 验证点 | 第一重（腰椎 LSS/LFH） | 第二重预期（GSE271018） |
|---|---|---|
| ① 铁死亡最活跃的细胞类型 | 巨噬细胞第 1（0.365） | 巨噬细胞仍居前 |
| ② 病 vs 对照 | 7 类细胞全升，巨噬增幅最大 | 退变组升高、巨噬显著 |
| ③ SPP1+ 亚群 | SPP1+ > SPP1-（P=2.9e-15） | 同向（该数据集正是 SPP1+ 巨噬的原始出处，**若命中是最强背书**） |

⚠️ 若结果不命中预期（比如该数据集中巨噬细胞不是第一），也不要紧——如实写进论文（"在两个独立单细胞数据集中一个复现"），反而显得客观。**不要挑数据。**

---

## 常见问题

| 问题 | 处理 |
|---|---|
| `readRDS` 报 "unknown class Assay5" | v5 对象撞老 Seurat → 发我 probe 报告，我给该文件单独的矩阵提取脚本 |
| 跑 03 时内存爆（killed） | 逐样本设计下只可能死在最大的 C3；重跑一次（断点续跑），仍死则发我，我给 C3 的降采样方案 |
| 没有注释列、兜底聚类太慢 | 正常，C3 大文件聚类可能要 20+ 分钟，耐心等；或发我报告，按原文补充材料手工映射 |
| wget 断流 | `wget -c` 支持断点续传，重跑 `bash 01_download.sh` 即可 |
| 服务器 R 版本极低（< 3.6） | 发我 env_report，别硬跑 |
