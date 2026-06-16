#!/bin/bash
#SBATCH --job-name=real_slasso
#SBATCH --partition=hamsi
#SBATCH --clusters=arf

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=00-03:00:00

#SBATCH --output=logs/real_slasso_%j.out
#SBATCH --error=logs/real_slasso_%j.err

set -euo pipefail

log_msg() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

log_msg "=================================================="
log_msg "Real-data SGLASSO benchmark job started"
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
DEFAULT_WORKDIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
WORKDIR="${WORKDIR:-$DEFAULT_WORKDIR}"

mkdir -p "${WORKDIR}/logs"
mkdir -p "${WORKDIR}/results/real_data"

cd "$WORKDIR"

mkdir -p "${WORKDIR}/tmp"
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

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

OUTDIR="${OUTDIR:-results/real_data}"

if [ -f "${SCRIPT_DIR}/run_real_slasso_all.R" ]; then
  RSCRIPT_FILE="${SCRIPT_DIR}/run_real_slasso_all.R"
  RFUNCTION_FILE="${SCRIPT_DIR}/sglasso_real_function.R"
elif [ -f "run_real_slasso_all.R" ]; then
  RSCRIPT_FILE="run_real_slasso_all.R"
  RFUNCTION_FILE="sglasso_real_function.R"
else
  echo "Cannot find run_real_slasso_all.R." >&2
  exit 1
fi

log_msg "Working directory      : $(pwd)"
log_msg "R_LIBS_USER            : ${R_LIBS_USER}"
log_msg "R_LIBS                 : ${R_LIBS}"
log_msg "OUTDIR                 : ${OUTDIR}"
log_msg "OMP_NUM_THREADS        : ${OMP_NUM_THREADS}"
log_msg "OPENBLAS_NUM_THREADS   : ${OPENBLAS_NUM_THREADS}"
log_msg "MKL_NUM_THREADS        : ${MKL_NUM_THREADS}"
log_msg "R executable           : $(which Rscript)"
log_msg "R version              : $(Rscript --version 2>&1)"
log_msg "R script               : ${RSCRIPT_FILE}"

log_msg "Checking required files..."
test -f "$RSCRIPT_FILE"
test -f "$RFUNCTION_FILE"
log_msg "Required files found."

log_msg "Starting R script..."
Rscript --vanilla "$RSCRIPT_FILE" --outdir="${OUTDIR}"
log_msg "R script finished."

log_msg "=================================================="
log_msg "Real-data sparse group lasso job finished"
log_msg "End time: $(date)"
log_msg "=================================================="
