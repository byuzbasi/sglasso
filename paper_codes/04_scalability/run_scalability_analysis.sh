#!/bin/bash
#SBATCH --job-name=sgl_scale
#SBATCH --partition=hamsi
#SBATCH --clusters=arf
#SBATCH --chdir=..

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=03-00:00:00

#SBATCH --output=logs/sgl_scale_%j.out
#SBATCH --error=logs/sgl_scale_%j.err

set -euo pipefail

log_msg() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

log_msg "=================================================="
log_msg "SGLASSO scalability analysis job started"
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

mkdir -p "${WORKDIR}/logs"
mkdir -p "${WORKDIR}/results/scalability"
mkdir -p "${WORKDIR}/results/scalability/checkpoints"
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

R_LIBRARY_DIR="${R_LIBRARY_DIR:-/arf/home/byuzbasi/R/x86_64-pc-linux-gnu-library/4.3}"
export R_LIBS_USER="${R_LIBRARY_DIR}"
export R_LIBS="${R_LIBRARY_DIR}"

# Keep BLAS/OpenMP single-threaded for fair timing.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

P_GRID="${P_GRID:-3000,10000,20000,50000}"
J_GRID="${J_GRID:-1000,2000,5000,10000}"
PJ="${PJ:-3}"
ACTIVE_GROUPS="${ACTIVE_GROUPS:-10}"
STRONG_FRAC="${STRONG_FRAC:-0.10}"
NREP="${NREP:-5}"
SEED="${SEED:-2026}"
ALPHA="${ALPHA:-0.5}"
GRPREG_ALPHA="${GRPREG_ALPHA:-1}"
D_VALUE="${D_VALUE:-0.5}"
NLAMBDA="${NLAMBDA:-30}"
STANDARDIZE="${STANDARDIZE:-FALSE}"
LAMBDA_MIN_RATIO="${LAMBDA_MIN_RATIO:-0.01}"
EPS="${EPS:-1e-4}"
MAXIT="${MAXIT:-1e5}"
SGLASSO_SCREEN="${SGLASSO_SCREEN:-SSR_fast}"
SGLASSO_TRANSFORM="${SGLASSO_TRANSFORM:-lazy}"
ADELIE_EARLY_EXIT="${ADELIE_EARLY_EXIT:-FALSE}"
CORES="${CORES:-${SLURM_CPUS_PER_TASK:-56}}"
FAST_DATA="${FAST_DATA:-TRUE}"
OUTDIR="${OUTDIR:-results/scalability}"

if [ -f "${SCRIPT_DIR}/run_scalability_analysis.R" ]; then
  RSCRIPT_FILE="${SCRIPT_DIR}/run_scalability_analysis.R"
elif [ -f "paper_codes/04_scalability/run_scalability_analysis.R" ]; then
  RSCRIPT_FILE="paper_codes/04_scalability/run_scalability_analysis.R"
elif [ -f "run_scalability_analysis.R" ]; then
  RSCRIPT_FILE="run_scalability_analysis.R"
else
  echo "Cannot find run_scalability_analysis.R." >&2
  exit 1
fi

log_msg "Working directory      : $(pwd)"
log_msg "R_LIBS_USER            : ${R_LIBS_USER}"
log_msg "R_LIBS                 : ${R_LIBS}"
log_msg "P_GRID                 : ${P_GRID}"
log_msg "J_GRID                 : ${J_GRID}"
log_msg "PJ                     : ${PJ}"
log_msg "ACTIVE_GROUPS          : ${ACTIVE_GROUPS}"
log_msg "STRONG_FRAC            : ${STRONG_FRAC}"
log_msg "NREP                   : ${NREP}"
log_msg "SEED                   : ${SEED}"
log_msg "ALPHA                  : ${ALPHA}"
log_msg "GRPREG_ALPHA           : ${GRPREG_ALPHA}"
log_msg "D_VALUE                : ${D_VALUE}"
log_msg "NLAMBDA                : ${NLAMBDA}"
log_msg "STANDARDIZE            : ${STANDARDIZE}"
log_msg "LAMBDA_MIN_RATIO       : ${LAMBDA_MIN_RATIO}"
log_msg "EPS                    : ${EPS}"
log_msg "MAXIT                  : ${MAXIT}"
log_msg "SGLASSO_SCREEN         : ${SGLASSO_SCREEN}"
log_msg "SGLASSO_TRANSFORM      : ${SGLASSO_TRANSFORM}"
log_msg "ADELIE_EARLY_EXIT      : ${ADELIE_EARLY_EXIT}"
log_msg "CORES                  : ${CORES}"
log_msg "FAST_DATA              : ${FAST_DATA}"
log_msg "OUTDIR                 : ${OUTDIR}"
log_msg "OMP_NUM_THREADS        : ${OMP_NUM_THREADS}"
log_msg "OPENBLAS_NUM_THREADS   : ${OPENBLAS_NUM_THREADS}"
log_msg "MKL_NUM_THREADS        : ${MKL_NUM_THREADS}"
log_msg "R executable           : $(which Rscript)"
log_msg "R version              : $(Rscript --version 2>&1)"
R_LIB_PATHS="$(Rscript --vanilla -e 'cat(paste(.libPaths(), collapse = " | "))')"
log_msg "R library paths        : ${R_LIB_PATHS}"
log_msg "R script               : ${RSCRIPT_FILE}"

log_msg "Starting scalability R script..."
Rscript --vanilla "$RSCRIPT_FILE" \
  --p_grid="${P_GRID}" \
  --J_grid="${J_GRID}" \
  --pj="${PJ}" \
  --active_groups="${ACTIVE_GROUPS}" \
  --strong_frac="${STRONG_FRAC}" \
  --nrep="${NREP}" \
  --seed="${SEED}" \
  --alpha="${ALPHA}" \
  --grpreg_alpha="${GRPREG_ALPHA}" \
  --d="${D_VALUE}" \
  --nlambda="${NLAMBDA}" \
  --standardize="${STANDARDIZE}" \
  --lambda_min_ratio="${LAMBDA_MIN_RATIO}" \
  --eps="${EPS}" \
  --maxit="${MAXIT}" \
  --sglasso_screen="${SGLASSO_SCREEN}" \
  --sglasso_transform="${SGLASSO_TRANSFORM}" \
  --adelie_early_exit="${ADELIE_EARLY_EXIT}" \
  --cores="${CORES}" \
  --fast_data="${FAST_DATA}" \
  --outdir="${OUTDIR}"
log_msg "Scalability R script finished."

log_msg "=================================================="
log_msg "SGLASSO scalability analysis job finished"
log_msg "End time: $(date)"
log_msg "=================================================="
