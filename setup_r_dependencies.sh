#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ENV_PREFIX=${H2AZ_R_ENV:-"${HOME}/.local/share/conda/envs/h2az-joint-states-r"}
CONDA=${CONDA_EXE:-}
if [[ -z "$CONDA" ]]; then CONDA=$(command -v conda || true); fi
[[ -n "$CONDA" ]] || { printf 'Conda is required for the isolated R environment.\n' >&2; exit 1; }
mkdir -p "$(dirname -- "$ENV_PREFIX")"
R_LIB="$ENV_PREFIX/lib/R/library"
if [[ -x "$ENV_PREFIX/bin/Rscript" ]] && \
   R_LIBS="$R_LIB" R_LIBS_USER="$R_LIB" "$ENV_PREFIX/bin/Rscript" -e \
   'quit(status=if (all(vapply(c("DiffBind","DESeq2","GenomicRanges","ggplot2","Rsamtools","csaw"), requireNamespace, logical(1), quietly=TRUE))) 0 else 1)' >/dev/null 2>&1; then
  printf 'Using existing R environment at %s\n' "$ENV_PREFIX"
elif [[ -d "$ENV_PREFIX" ]]; then
  "$CONDA" -y env update -p "$ENV_PREFIX" -f "$ROOT/setup/r-dependencies.yml"
else
  "$CONDA" -y env create -p "$ENV_PREFIX" -f "$ROOT/setup/r-dependencies.yml"
fi
R_LIBS="$R_LIB" R_LIBS_USER="$R_LIB" "$ENV_PREFIX/bin/Rscript" -e \
  'for(x in c("DiffBind","DESeq2","GenomicRanges","ggplot2","Rsamtools","csaw")) { stopifnot(requireNamespace(x, quietly=TRUE)); cat(x, as.character(packageVersion(x)), "\n") }'
printf 'Rscript installed at %s/bin/Rscript\n' "$ENV_PREFIX"
printf 'Select it explicitly with H2AZ_USE_CONDA_R=1 bash run_pipeline.sh.\n'
