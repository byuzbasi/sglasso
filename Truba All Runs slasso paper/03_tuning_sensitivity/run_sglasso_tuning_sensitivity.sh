#!/bin/bash
#SBATCH --job-name=sgl_tuning
#SBATCH --partition=hamsi
#SBATCH --clusters=arf

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=00-12:00:00

#SBATCH --output=logs/sgl_tuning_%j.out
#SBATCH --error=logs/sgl_tuning_%j.err

set -euo pipefail

log_msg() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

log_msg "=================================================="
log_msg "SGLASSO tuning sensitivity job started"
log_msg "Job ID          : ${SLURM_JOB_ID:-NA}"
log_msg "Job Name        : ${SLURM_JOB_NAME:-NA}"
log_msg "Cluster         : ${SLURM_CLUSTER_NAME:-NA}"
log_msg "Partition       : ${SLURM_JOB_PARTITION:-NA}"
log_msg "Hostname        : $(hostname)"
log_msg "Node list       : ${SLURM_NODELIST:-NA}"
log_msg "Submit directory: ${SLURM_SUBMIT_DIR:-NA}"
log_msg "Allocated CPUs  : ${SLURM_CPUS_PER_TASK:-NA}"
log_msg "Start time      : $(date)"
log_msg "=================================================="

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR="${WORKDIR:-$PWD}"

mkdir -p "${WORKDIR}/logs"
mkdir -p "${WORKDIR}/results/tuning_sensitivity"
mkdir -p "${WORKDIR}/tmp"

cd "$WORKDIR"

export TMPDIR="${WORKDIR}/tmp"
export TEMP="${WORKDIR}/tmp"
export TMP="${WORKDIR}/tmp"

MODULE_AVAILABLE=0
if command -v module >/dev/null 2>&1; then
  MODULE_AVAILABLE=1
  module purge 2>/dev/null || true
  module load "${R_MODULE:-apps/R/4.3.0-gcc-11.3.1}"
fi

if [ "$MODULE_AVAILABLE" -eq 1 ]; then
  log_msg "Loaded modules:"
  module list 2>&1
else
  log_msg "Environment modules are not available; using current R environment."
fi

export R_LIBS_USER="${R_LIBS_USER:-/arf/home/byuzbasi/R/x86_64-pc-linux-gnu-library/4.3}"
export R_LIBS="${R_LIBS:-$R_LIBS_USER}"

# Keep BLAS/OpenMP single-threaded while foreach uses SLURM workers.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

NREP="${NREP:-5}"
SEED="${SEED:-2026}"
CORES="${CORES:-${SLURM_CPUS_PER_TASK:-1}}"
OUTDIR="${OUTDIR:-results/tuning_sensitivity}"

if [ -f "${SCRIPT_DIR}/run_sglasso_tuning_sensitivity.R" ]; then
  RSCRIPT_FILE="${SCRIPT_DIR}/run_sglasso_tuning_sensitivity.R"
elif [ -f "run_sglasso_tuning_sensitivity.R" ]; then
  RSCRIPT_FILE="run_sglasso_tuning_sensitivity.R"
else
  echo "Cannot find run_sglasso_tuning_sensitivity.R." >&2
  exit 1
fi

log_msg "Working directory      : $(pwd)"
log_msg "R_LIBS_USER            : ${R_LIBS_USER}"
log_msg "R_LIBS                 : ${R_LIBS}"
log_msg "NREP                   : ${NREP}"
log_msg "SEED                   : ${SEED}"
log_msg "CORES                  : ${CORES}"
log_msg "OUTDIR                 : ${OUTDIR}"
log_msg "OMP_NUM_THREADS        : ${OMP_NUM_THREADS}"
log_msg "OPENBLAS_NUM_THREADS   : ${OPENBLAS_NUM_THREADS}"
log_msg "MKL_NUM_THREADS        : ${MKL_NUM_THREADS}"
log_msg "R executable           : $(which Rscript)"
log_msg "R version              : $(Rscript --version 2>&1)"
log_msg "R library paths        : $(Rscript -e 'cat(paste(.libPaths(), collapse=\" | \"))')"
log_msg "R script               : ${RSCRIPT_FILE}"

log_msg "Starting tuning sensitivity R script..."
Rscript --vanilla "$RSCRIPT_FILE" \
  --nrep="${NREP}" \
  --seed="${SEED}" \
  --cores="${CORES}" \
  --outdir="${OUTDIR}"
log_msg "Tuning sensitivity R script finished."

log_msg "=================================================="
log_msg "SGLASSO tuning sensitivity job finished"
log_msg "End time: $(date)"
log_msg "=================================================="
