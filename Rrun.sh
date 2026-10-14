#!/bin/bash
# 干净环境运行用户的 R 4.5.3（规避 Git Bash 环境变量导致的 R 启动失败）
# 用法: bash scripts/Rrun.sh <script.R> [args...]
env -i \
  PATH="/d/A生信/R-4.5.3/bin:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin:/c/Windows/System32:/c/Windows" \
  SystemRoot="C:\\Windows" \
  TEMP="C:\\Users\\10215\\AppData\\Local\\Temp" \
  TMP="C:\\Users\\10215\\AppData\\Local\\Temp" \
  USERPROFILE="C:\\Users\\10215" \
  HOMEDRIVE="C:" HOMEPATH="\\Users\\10215" \
  APPDATA="C:\\Users\\10215\\AppData\\Roaming" \
  LOCALAPPDATA="C:\\Users\\10215\\AppData\\Local" \
  ComSpec="C:\\Windows\\System32\\cmd.exe" \
  "/d/A生信/R-4.5.3/bin/Rscript.exe" "$@"
