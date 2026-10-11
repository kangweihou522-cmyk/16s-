#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
export R_LIBS_USER="$HOME/R/16S-workstation-library${R_LIBS_USER:+:$R_LIBS_USER}"
exec python3 launch.py "$@"
