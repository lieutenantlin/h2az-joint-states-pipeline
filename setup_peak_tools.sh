#!/usr/bin/env bash
# Provision or update peak-calling tools and the external SEACR 1.3 scripts.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ENV_PREFIX=${H2AZ_TOOLS_ENV:-"${HOME}/.local/share/conda/envs/h2az-joint-states"}
TOOL_HOME=${H2AZ_TOOL_HOME:-"${HOME}/.local/share/h2az-joint-states"}
CONDA=${CONDA_EXE:-}
if [[ -z "$CONDA" ]]; then CONDA=$(command -v conda || true); fi
[[ -n "$CONDA" ]] || { printf 'Conda is required to provision samtools, bedtools, and MACS3.\n' >&2; exit 1; }
mkdir -p "$(dirname -- "$ENV_PREFIX")" "$TOOL_HOME"
if [[ -x "$ENV_PREFIX/bin/samtools" && -x "$ENV_PREFIX/bin/bedtools" && -x "$ENV_PREFIX/bin/macs3" ]]; then
  printf 'Using existing peak-tools environment at %s\n' "$ENV_PREFIX"
elif [[ -d "$ENV_PREFIX" ]]; then
  "$CONDA" -y env update -p "$ENV_PREFIX" -f "$ROOT/setup/peak-tools.yml"
else
  "$CONDA" -y env create -p "$ENV_PREFIX" -f "$ROOT/setup/peak-tools.yml"
fi
# Keep the official SEACR 1.3 files outside the source tree and avoid needless downloads on rerun.
if [[ ! -s "$TOOL_HOME/SEACR_1.3.sh" ]]; then
  curl -fL --retry 3 https://raw.githubusercontent.com/FredHutch/SEACR/master/SEACR_1.3.sh -o "$TOOL_HOME/SEACR_1.3.sh"
fi
if [[ ! -s "$TOOL_HOME/SEACR_1.3.R" ]]; then
  curl -fL --retry 3 https://raw.githubusercontent.com/FredHutch/SEACR/master/SEACR_1.3.R -o "$TOOL_HOME/SEACR_1.3.R"
fi
chmod 755 "$TOOL_HOME/SEACR_1.3.sh"
cat <<EOF
Peak-calling tools installed in:
  $ENV_PREFIX
SEACR 1.3 installed in:
  $TOOL_HOME/SEACR_1.3.sh
Use them with:
  export PATH="$ENV_PREFIX/bin:\$PATH"
  export SEACR="$TOOL_HOME/SEACR_1.3.sh"
  export MACS3="$ENV_PREFIX/bin/macs3"
  export THREADS=4
EOF
"$ENV_PREFIX/bin/samtools" --version | head -n 1
"$ENV_PREFIX/bin/bedtools" --version
"$ENV_PREFIX/bin/macs3" --version
