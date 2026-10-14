# 服务器 CellChat 复现指令（GSE271018 第二队列）

> 目的：在同学服务器上对 **GSE271018（51,686 细胞，3 创伤对照 vs 6 退变）** 复现 CellChat
> 细胞通讯分析，验证本地第一队列（GSE294458/267819）发现的 **SPP1⁺巨噬细胞 → 基质细胞**
> 信号轴，构成两队列交叉证据。
>
> 本地第一队列结果已出（SPP1 信号完全特异于 SPP1⁺巨噬细胞；肥厚组 SPP1⁺发出 L-R 对 152
> 条 vs 对照 96 条），服务器跑完即可对照。

---

## 0. 前置条件（已具备）

- ✅ 服务器 `~/GSE271018/` 已有 9 个 RDS 文件（Tendon_C1/C3/C5 = 对照，Tendon_N0-N10 = 退变）
- ✅ 服务器可联网装新包（允许装新包，**不能更新旧包**——本包全部遵守）
- ✅ SSH 端口 25100（`ssh -p 25100 weijian@yxszxxsx.com`，注意 scp 用大写 `-P`）
- ⚠️ **用户名是 `weijian`**（2026-10-11 更正）：`weijian@shpc` 中的 `shpc` 是服务器主机名而非用户名一部分，误用会被服务器当不存在用户拒绝（Permission denied）

## 1. 上传分析包（本地 PowerShell）

```powershell
scp -P 25100 -r "C:\Users\10215\WorkBuddy\2026-10-03-00-10-07\server_cellchat_GSE271018" "weijian@yxszxxsx.com:~/"
```

> ⚠️ 必须在 PowerShell 里用**完整路径**（如上），或先 `cd C:\Users\10215\WorkBuddy\2026-10-03-00-10-07`
> 再用相对路径——在 `C:\Users\10215` 主目录直接 scp 相对路径会报 "No such file or directory"。

## 2. 安装 CellChat（服务器，约 10-30 分钟）

```bash
ssh -p 25100 weijian@yxszxxsx.com     # 登录服务器
cd ~/server_cellchat_GSE271018
Rscript 01_install_pkgs.R 2>&1 | tee install_log.txt
```

**要点（脚本已内置，无需手动操作）**：
- 所有新包装进独立个人库 `~/Rlibs_new`，**完全不碰系统库**
- `remotes::install_github(..., upgrade = "never")` 保证不升级任何已有包
- 若 GitHub 直连慢/失败：重跑一次即可（已装依赖不会重复装）

**成功标志**：日志末尾出现 `=== CellChat 1.x.x 安装成功 ===`

## 3. 跑 CellChat 主分析（服务器，约 30-90 分钟）

```bash
cd ~                                                    # 必须在主目录运行（脚本内部为绝对路径）
nohup Rscript server_cellchat_GSE271018/02_cellchat.R > run_cellchat.log 2>&1 &
tail -f run_cellchat.log        # 看进度；退出 tail 用 Ctrl+C
```

> ⚠️ 不要 `cd ~/server_cellchat_GSE271018` 后再跑——早期版本用相对路径会在错误位置找数据，报 `length(rds_files) > 0 is not TRUE`。2026-10-11 起脚本已改绝对路径，但保持从主目录运行的习惯最稳。

**流程**（脚本自动完成）：
1. 逐样本读 RDS → 注释（原文注释列优先，标志基因面板兜底）→ 巨噬细胞 SPP1⁺/⁻ 细分
2. 按组拼接（Control 3 样本 / Degenerative 6 样本）→ 每组独立 CellChat 推断
3. 合并比较 → 差异气泡图 + 通路排序 + SPP1 网络图 + L-R 明细表

**注意**：本机 CellChat 1.6.1 踩过并已修复的坑（服务器脚本同样已修）：
- `projectData` 是可选步骤，已删（`computeCommunProb` 默认用原始数据）
- 1.6.x 无 `computeNetwork`，已改 `aggregateNet` + `netAnalysis_computeCentrality`

**2026-10-11 上传前核查记录（已通过）**：
- ✅ 两个脚本 R 语法解析通过
- ✅ 分组正则 `Tendon_C/N` 与已在服务器跑通的 03 脚本逐字一致
- ✅ `merge_group()` 稀疏矩阵对齐经合成数据 5 项断言实测通过；**修复了一个严重 bug**：`sparseMatrix()` 缺省创建逻辑型 ngCMatrix 会把表达量强制转 0/1，已改为显式 `x = numeric(0)`（dgCMatrix）——**请务必上传修复后的版本**（02_cellchat.R 文件日期 ≥ 2026-10-11）

**成功标志**：日志出现 `=== 02_cellchat done ===`，且 `~/GSE271018_cellchat/` 里有图和 CSV。

## 4. 取回结果（本地 PowerShell）

```powershell
cd C:\Users\10215\WorkBuddy\2026-10-03-00-10-07
scp -P 25100 -r "weijian@yxszxxsx.com:~/GSE271018_cellchat" .
```

## 5. 复核要点（结果回来后我帮你做）

| 检查项 | 预期 |
|---|---|
| SPP1 通路发送者 | 应**完全/主要**来自 SPP1⁺巨噬细胞（与第一队列一致） |
| Degenerative vs Control 的 SPP1⁺ L-R 对数 | 退变组应**多于**对照组（第一队列 152 vs 96） |
| SPP1-CD44 / SPP1-整合素 → Fibroblast | 两组均显著 |
| 上传回来的图 A/B/C | 与本地 `scRNA_out/cellchat/` 同名图可并排放 Figure |

## 常见问题

- **装包失败提示缺 Bio\* 依赖**：`BiocManager::install("包名", lib="~/Rlibs_new")` 补装后重跑 01
- **02 中途断掉**：回到主目录重新 `nohup Rscript server_cellchat_GSE271018/02_cellchat.R > run_cellchat.log 2>&1 &`——**有断点续跑**（已完成组会秒过），只重算未完成的组
- **不影响服务器其他用户**：所有计算单核、内存峰值 <50GB，包装在个人库
