#!/bin/bash
#SBATCH --job-name=sim_slasso
#SBATCH --account=byuzbasi
#SBATCH --partition=hamsi
#SBATCH --clusters=arf
#SBATCH --chdir=..

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=02-00:00:00

#SBATCH --mail-user=b.yzb@hotmail.com
#SBATCH --mail-type=END,FAIL

#SBATCH --output=logs/sim_slasso_%j.out
#SBATCH --error=logs/sim_slasso_%j.err

set -euo pipefail

log_msg() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

log_msg "=================================================="
log_msg "SGLASSO simulation job started"
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
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
if [ "$(basename "$SUBMIT_DIR")" = "paper_codes" ]; then
  DEFAULT_WORKDIR="$(cd "${SUBMIT_DIR}/.." && pwd)"
elif [ -d "${SUBMIT_DIR}/paper_codes" ]; then
  DEFAULT_WORKDIR="$SUBMIT_DIR"
else
  DEFAULT_WORKDIR="$(pwd)"
fi
WORKDIR="${WORKDIR:-$DEFAULT_WORKDIR}"
OUTDIR="${OUTDIR:-results/main_simulation}"

mkdir -p "${WORKDIR}/logs"
mkdir -p "${WORKDIR}/${OUTDIR}/rds"
mkdir -p "${WORKDIR}/${OUTDIR}/tables"
mkdir -p "${WORKDIR}/tmp"

cd "$WORKDIR"

export TMPDIR="${WORKDIR}/tmp"
export TEMP="${WORKDIR}/tmp"
export TMP="${WORKDIR}/tmp"

export MODULEPATH="${MODULEPATH:-}:/arf/sw/ek-moduller/modules/all"

module purge 2>/dev/null || true
module load apps/R/4.3.0-gcc-11.3.1

log_msg "Loaded modules:"
module list 2>&1

R_LIBRARY_DIR="${R_LIBRARY_DIR:-/arf/home/byuzbasi/R/x86_64-pc-linux-gnu-library/4.3}"
export R_LIBS_USER="${R_LIBRARY_DIR}"
export R_LIBS="${R_LIBRARY_DIR}"

# foreach uses SLURM_CPUS_PER_TASK workers. Keep BLAS/OpenMP single-threaded
# to avoid CPU oversubscription.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

REPEATNUM="${REPEATNUM:-100}"

if [ -f "${SCRIPT_DIR}/02_main_simulation/run_sim_slasso.all.R" ]; then
  RSCRIPT_FILE="${SCRIPT_DIR}/02_main_simulation/run_sim_slasso.all.R"
elif [ -f "${SCRIPT_DIR}/run_sim_slasso.all.R" ]; then
  RSCRIPT_FILE="${SCRIPT_DIR}/run_sim_slasso.all.R"
elif [ -f "paper_codes/02_main_simulation/run_sim_slasso.all.R" ]; then
  RSCRIPT_FILE="paper_codes/02_main_simulation/run_sim_slasso.all.R"
elif [ -f "02_main_simulation/run_sim_slasso.all.R" ]; then
  RSCRIPT_FILE="02_main_simulation/run_sim_slasso.all.R"
elif [ -f "run_sim_slasso.all.R" ]; then
  RSCRIPT_FILE="run_sim_slasso.all.R"
else
  echo "Cannot find run_sim_slasso.all.R." >&2
  exit 1
fi

if [ ! -f "${SCRIPT_DIR}/02_main_simulation/sglasso_sim_function_full_tuning.R" ] &&
   [ ! -f "${SCRIPT_DIR}/sglasso_sim_function_full_tuning.R" ] &&
   [ ! -f "paper_codes/02_main_simulation/sglasso_sim_function_full_tuning.R" ] &&
   [ ! -f "02_main_simulation/sglasso_sim_function_full_tuning.R" ] &&
   [ ! -f "sglasso_sim_function_full_tuning.R" ]; then
  echo "Cannot find sglasso_sim_function_full_tuning.R." >&2
  exit 1
fi

log_msg "Working directory      : $(pwd)"
log_msg "R_LIBS_USER            : ${R_LIBS_USER}"
log_msg "R_LIBS                 : ${R_LIBS}"
log_msg "REPEATNUM              : ${REPEATNUM}"
log_msg "OUTDIR                 : ${OUTDIR}"
log_msg "SLURM_CPUS_PER_TASK    : ${SLURM_CPUS_PER_TASK:-NA}"
log_msg "OMP_NUM_THREADS        : ${OMP_NUM_THREADS}"
log_msg "OPENBLAS_NUM_THREADS   : ${OPENBLAS_NUM_THREADS}"
log_msg "MKL_NUM_THREADS        : ${MKL_NUM_THREADS}"
log_msg "R executable           : $(which Rscript)"
log_msg "R version              : $(Rscript --version 2>&1)"
R_LIB_PATHS="$(Rscript --vanilla -e 'cat(paste(.libPaths(), collapse = " | "))')"
log_msg "R library paths        : ${R_LIB_PATHS}"
log_msg "R script               : ${RSCRIPT_FILE}"

log_msg "Starting R simulation script..."
REPEATNUM="${REPEATNUM}" OUTDIR="${OUTDIR}" Rscript --vanilla "$RSCRIPT_FILE"
log_msg "R simulation script finished."

log_msg "=================================================="
log_msg "SGLASSO simulation job finished"
log_msg "End time: $(date)"
log_msg "=================================================="
