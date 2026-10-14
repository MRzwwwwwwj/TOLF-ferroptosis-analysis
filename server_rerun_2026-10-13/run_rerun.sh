#!/bin/bash
cd ~
Rscript rerun_2026-10-13/rerun_03_ferroptosis_validate.R > rerun_2026-10-13/log_03_ferroptosis.txt 2>&1
S1=$?
Rscript rerun_2026-10-13/rerun_02_cellchat.R > rerun_2026-10-13/log_02_cellchat.txt 2>&1
S2=$?
echo "s1=$S1 s2=$S2" > rerun_2026-10-13/STATUS.txt
