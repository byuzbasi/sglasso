#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR="${WORKDIR:-$(cd "${SCRIPT_DIR}/.." && pwd)}"
cd "$SCRIPT_DIR"

mkdir -p "${WORKDIR}/logs"
mkdir -p "${WORKDIR}/results/real_data"
mkdir -p "${WORKDIR}/results/main_simulation"
mkdir -p "${WORKDIR}/results/main_simulation/rds"
mkdir -p "${WORKDIR}/results/tuning_sensitivity"
mkdir -p "${WORKDIR}/results/tuning_sensitivity/task_rds"
mkdir -p "${WORKDIR}/results/scalability"
mkdir -p "${WORKDIR}/results/scalability/checkpoints"
mkdir -p "${WORKDIR}/tmp"

echo "Submitting real-data run..."
sbatch 01_real_data/run_real_slasso_all.sh

echo "Submitting main simulation run..."
sbatch 02_main_simulation/run_sim_slasso_all.sh

echo "Submitting tuning sensitivity run..."
sbatch 03_tuning_sensitivity/run_sglasso_tuning_sensitivity.sh

echo "Submitting scalability run..."
sbatch 04_scalability/run_scalability_analysis.sh

echo "All jobs submitted."
