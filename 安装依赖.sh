#!/usr/bin/env bash
set -euo pipefail
# Dedicated binary R packages: do not install/upgrade global R or mix Ubuntu/CRAN builds.
cd -- "$(dirname -- "$0")"
r_env_name="${WORKSTATION_R_ENV:-16s-workstation-r}"
[[ "$r_env_name" =~ ^[A-Za-z0-9_][A-Za-z0-9_.-]*$ ]] || { echo 'R环境名无效' >&2; exit 1; }
conda_bin="${CONDA_EXE:-}"
if [[ ! -x "$conda_bin" ]]; then conda_bin="$(type -P conda || true)"; fi
if [[ ! -x "$conda_bin" ]]; then
  for candidate in "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda" /opt/conda/bin/conda; do
    if [[ -x "$candidate" ]]; then conda_bin="$candidate"; break; fi
  done
fi
[[ -x "$conda_bin" ]] || { echo '找不到conda，请从已安装Miniconda的Ubuntu终端执行此脚本' >&2; exit 1; }
if [[ "$EUID" -eq 0 ]]; then APT=(apt-get); else APT=(sudo apt-get); fi
"${APT[@]}" update
"${APT[@]}" install -y python3 ca-certificates fontconfig fonts-dejavu-core
# Binary packages are resolved together from one channel, with normal TLS/hash verification.
if "$conda_bin" env list --json | python3 -c 'import json,os,sys; name=sys.argv[1]; sys.exit(0 if any(os.path.basename(p)==name for p in json.load(sys.stdin)["envs"]) else 1)' "$r_env_name"; then
  operation=install
else
  operation=create
fi
"$conda_bin" "$operation" --name "$r_env_name" --override-channels -c conda-forge --strict-channel-priority -y \
  'r-base=4.4' r-ggplot2 r-patchwork r-svglite r-xml2 r-jsonlite r-systemfonts r-multcompview r-mass cairo fontconfig
# Make full DejaVu Serif and existing Windows Times available to either fontconfig.
font_dir="$HOME/.local/share/fonts/16S-workstation"
mkdir -p "$font_dir"
for font in /usr/share/fonts/truetype/dejavu/*.ttf /mnt/c/Windows/Fonts/times*.ttf; do
  [[ -f "$font" ]] || continue
  [[ -e "$font_dir/$(basename "$font")" ]] || cp "$font" "$font_dir/"
done
fc-cache -f "$font_dir"
# Drop inherited user/system R search paths and compiler/library overrides for checks.
env -u R_HOME -u LD_LIBRARY_PATH -u LD_PRELOAD -u PKG_CONFIG_PATH -u LIBRARY_PATH -u CPATH -u C_INCLUDE_PATH -u CPLUS_INCLUDE_PATH -u FONTCONFIG_PATH -u FONTCONFIG_FILE \
  R_LIBS=/dev/null R_LIBS_USER=/dev/null R_LIBS_SITE=/dev/null R_PROFILE=/dev/null R_PROFILE_USER=/dev/null R_ENVIRON=/dev/null R_ENVIRON_USER=/dev/null \
  "$conda_bin" run --no-capture-output -n "$r_env_name" /usr/bin/env -u R_HOME -u LD_LIBRARY_PATH -u LD_PRELOAD R_LIBS=/dev/null R_LIBS_USER=/dev/null R_LIBS_SITE=/dev/null R_PROFILE=/dev/null R_PROFILE_USER=/dev/null R_ENVIRON=/dev/null R_ENVIRON_USER=/dev/null Rscript --vanilla check_r.R 'DejaVu Serif'
printf '%s\n' "专用R环境 $r_env_name 安装及检查完成。请运行 bash 启动.sh。"
