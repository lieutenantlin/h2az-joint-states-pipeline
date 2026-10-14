# Runtime tool and optional Conda-R selection for the pipeline entry points.
THREADS=${THREADS:-4}
H2AZ_TOOLS_ENV=${H2AZ_TOOLS_ENV:-"${HOME}/.local/share/conda/envs/h2az-joint-states"}
H2AZ_R_ENV=${H2AZ_R_ENV:-"${HOME}/.local/share/conda/envs/h2az-joint-states-r"}
SEACR=${SEACR:-"${HOME}/.local/share/h2az-joint-states/SEACR_1.3.sh"}
MACS3=${MACS3:-"${H2AZ_TOOLS_ENV}/bin/macs3"}
H2AZ_RSCRIPT=${H2AZ_RSCRIPT:-Rscript}
if [[ "${H2AZ_USE_CONDA_R:-0}" == 1 ]]; then
  [[ -x "${H2AZ_R_ENV}/bin/Rscript" ]] || { printf 'Conda Rscript not found: %s/bin/Rscript\n' "$H2AZ_R_ENV" >&2; return 1 2>/dev/null || exit 1; }
  PATH="${H2AZ_R_ENV}/bin:${PATH}"
  H2AZ_RSCRIPT="${H2AZ_R_ENV}/bin/Rscript"
  R_LIBS="${H2AZ_R_ENV}/lib/R/library"
  R_LIBS_USER="$R_LIBS"
  export R_LIBS R_LIBS_USER
fi
if [[ -d "${H2AZ_TOOLS_ENV}/bin" ]]; then
  PATH="${H2AZ_TOOLS_ENV}/bin:${PATH}"
fi
export THREADS SEACR MACS3 H2AZ_TOOLS_ENV H2AZ_R_ENV H2AZ_RSCRIPT PATH
