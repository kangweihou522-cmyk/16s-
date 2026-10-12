#!/usr/bin/env bash
set -euo pipefail
# Install the official Bioconda package in the selected LEfSe environment only.
lefse_env_name="${WORKSTATION_LEFSE_ENV:-lefse}"
[[ "$lefse_env_name" =~ ^[A-Za-z0-9_][A-Za-z0-9_.-]*$ ]] || { echo 'LEfSe环境名无效' >&2; exit 1; }
conda_bin="${CONDA_EXE:-}"
if [[ ! -x "$conda_bin" ]]; then conda_bin="$(type -P conda || true)"; fi
if [[ ! -x "$conda_bin" ]]; then
  for candidate in "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda" /opt/conda/bin/conda; do
    if [[ -x "$candidate" ]]; then conda_bin="$candidate"; break; fi
  done
fi
[[ -x "$conda_bin" ]] || { echo '找不到conda，请从已安装Miniconda的Ubuntu终端执行此脚本' >&2; exit 1; }
if "$conda_bin" env list --json | python3 -c 'import json,os,sys; name=sys.argv[1]; sys.exit(0 if any(os.path.basename(p)==name for p in json.load(sys.stdin)["envs"]) else 1)' "$lefse_env_name"; then
  operation=install
else
  operation=create
fi
"$conda_bin" "$operation" --name "$lefse_env_name" --override-channels -c conda-forge -c bioconda --strict-channel-priority -y 'lefse=1.1.2'
"$conda_bin" run --no-capture-output -n "$lefse_env_name" python -c '
import importlib.metadata,os,shutil,subprocess,sys
commands=("lefse_format_input.py","lefse_run.py","lefse_plot_res.py","lefse_plot_cladogram.py")
print("LEfSe",importlib.metadata.version("lefse"),"environment:",sys.prefix,flush=True)
for name in commands:
    path=shutil.which(name)
    if not path or os.path.commonpath([os.path.abspath(path),os.path.abspath(sys.prefix)])!=os.path.abspath(sys.prefix):
        raise SystemExit("环境内缺少命令："+name)
    check=subprocess.run([path,"-h"],capture_output=True,text=True,timeout=90)
    if check.returncode:
        raise SystemExit(name+" 启动失败：\n"+check.stdout+check.stderr)
    print("OK:",path,flush=True)
from lefse.lefse import init
init()
print("OK: R统计依赖已加载",flush=True)
'
echo "LEfSe安装与启动检查通过。网页LEfSe conda环境填写：$lefse_env_name"
