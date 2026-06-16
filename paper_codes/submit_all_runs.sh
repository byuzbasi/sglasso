#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR="${WORKDIR:-$(cd "${SCRIPT_DIR}/.." && pwd)}"
cd "$WORKDIR"

mkdir -p logs
mkdir -p results/real_data
mkdir -p results/main_simulation
mkdir -p results/tuning_sensitivity
mkdir -p results/scalability
mkdir -p tmp

echo "Submitting real-data run..."
sbatch paper_codes/01_real_data/run_real_slasso_all.sh

echo "Submitting main simulation run..."
sbatch paper_codes/02_main_simulation/run_sim_slasso_all.sh

echo "Submitting tuning sensitivity run..."
sbatch paper_codes/03_tuning_sensitivity/run_sglasso_tuning_sensitivity.sh

echo "Submitting scalability run..."
sbatch paper_codes/04_scalability/run_scalability_analysis.sh

echo "All jobs submitted."
