#!/usr/bin/env bash
set -euo pipefail
# Does not install or replace QIIME 2 environments. Use your existing compatible environment.
if [[ "$EUID" -eq 0 ]]; then APT=(apt-get); else APT=(sudo apt-get); fi
"${APT[@]}" update
"${APT[@]}" install -y python3 r-base r-base-dev libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev libpng-dev libtiff5-dev libjpeg-dev fonts-dejavu-core
Rscript - <<'R'
local_lib <- path.expand('~/R/16S-workstation-library')
dir.create(local_lib,recursive=TRUE,showWarnings=FALSE)
.libPaths(c(local_lib,.libPaths()))
pkgs <- c('ggplot2','patchwork','svglite','xml2','jsonlite','systemfonts','multcompView')
missing <- pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) install.packages(missing,repos='https://cloud.r-project.org',lib=local_lib)
if(any(!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE))) stop('R包安装未完成')
R
# Reuse legally installed Windows Times New Roman; do not download fonts from third parties.
if [[ -f /mnt/c/Windows/Fonts/times.ttf ]]; then
  mkdir -p "$HOME/.local/share/fonts/16S-workstation"
  for font in /mnt/c/Windows/Fonts/times*.ttf; do cp -n "$font" "$HOME/.local/share/fonts/16S-workstation/"; done
  fc-cache -f "$HOME/.local/share/fonts/16S-workstation"
fi
printf '%s\n' 'R/字体依赖安装完成。请确认已有 QIIME 2 环境与兼容的分类器；本脚本不会自动下载它们。'
