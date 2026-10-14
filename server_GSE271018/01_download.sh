#!/bin/bash
# 01_download.sh — 下载 GSE271018 全部数据（约 6.3 GB）
# 用法：bash 01_download.sh
# 建议放在 screen/tmux 或 nohup 里跑：
#   nohup bash 01_download.sh > download.log 2>&1 &
set -e

mkdir -p GSE271018
cd GSE271018

echo "==> [1/3] series_matrix（样本分组信息，几 KB）"
wget -c "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE271nnn/GSE271018/matrix/GSE271018_series_matrix.txt.gz"

echo "==> [2/3] RAW.tar（9 个 RDS，6.3 GB，耗时取决于带宽）"
wget -c "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE271nnn/GSE271018/suppl/GSE271018_RAW.tar"

echo "==> [3/3] 解包"
tar -xf GSE271018_RAW.tar

echo ""
echo "==> 下载完成，文件清单："
ls -lh

echo ""
echo "==> 样本分组（从 series_matrix 提取，请核对你的分组假设）:"
zcat GSE271018_series_matrix.txt.gz | grep -E "!Sample_title|!Sample_source_name_ch1" | cut -c1-400
