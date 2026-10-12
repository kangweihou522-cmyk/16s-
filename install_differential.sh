#!/usr/bin/env bash
set -euo pipefail
# Separate binary R environment; does not alter QIIME or the plotting R environment.
ancom_env_name="${WORKSTATION_ANCOM_ENV:-16s-workstation-ancom}"
[[ "$ancom_env_name" =~ ^[A-Za-z0-9_][A-Za-z0-9_.-]*$ ]] || { echo 'ANCOM环境名无效' >&2; exit 1; }
conda_bin="${CONDA_EXE:-}"
if [[ ! -x "$conda_bin" ]]; then conda_bin="$(type -P conda || true)"; fi
if [[ ! -x "$conda_bin" ]]; then
  for candidate in "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda" /opt/conda/bin/conda; do
    if [[ -x "$candidate" ]]; then conda_bin="$candidate"; break; fi
  done
fi
[[ -x "$conda_bin" ]] || { echo '找不到conda，请从已安装Miniconda的Ubuntu终端执行此脚本' >&2; exit 1; }
if "$conda_bin" env list --json | python3 -c 'import json,os,sys; name=sys.argv[1]; sys.exit(0 if any(os.path.basename(p)==name for p in json.load(sys.stdin)["envs"]) else 1)' "$ancom_env_name"; then
  operation=install
else
  operation=create
fi
"$conda_bin" "$operation" --name "$ancom_env_name" --override-channels -c conda-forge -c bioconda --strict-channel-priority -y \
  'r-base>=4.4,<4.6' 'bioconductor-ancombc>=2.8,<3' bioconductor-phyloseq r-jsonlite r-ggplot2 r-svglite
"$conda_bin" run --no-capture-output -n "$ancom_env_name" /usr/bin/env -u R_HOME -u LD_LIBRARY_PATH -u LD_PRELOAD \
  R_LIBS=/dev/null R_LIBS_USER=/dev/null R_LIBS_SITE=/dev/null R_PROFILE=/dev/null R_PROFILE_USER=/dev/null \
  R_ENVIRON=/dev/null R_ENVIRON_USER=/dev/null Rscript --vanilla -e \
  'for(p in c("ANCOMBC","phyloseq","jsonlite","ggplot2","svglite")){stopifnot(requireNamespace(p,quietly=TRUE));cat(p,as.character(packageVersion(p)),"\n")};stopifnot("ancombc2" %in% getNamespaceExports("ANCOMBC"));cat(R.version.string,"\n")'
