#!/bin/bash
#SBATCH --job-name=bacdive_test_eval
#SBATCH --output=eval_test_%j.out
#SBATCH --error=eval_test_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --partition=gpu
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:h100:1
#SBATCH --mem=80G
#SBATCH --nodelist=bioinf040

# BacDive test set inference: gen from FASTA + model predict, 800 batches
# Model uses ~79GB VRAM (fits on 1x H100 80GB)
# CPU memory: ~30GB for FASTA I/O + data gen + prediction storage

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -f "${SCRIPT_DIR}/.env" ]]; then
  set -a
  # shellcheck disable=SC1091
  source "${SCRIPT_DIR}/.env"
  set +a
fi

BACDIVE_CONDA_SH="${BACDIVE_CONDA_SH:-${HOME}/miniconda3/etc/profile.d/conda.sh}"
BACDIVE_TF_ENV="${BACDIVE_TF_ENV:-tf}"
EVAL_SCRIPT_PATH="${EVAL_SCRIPT_PATH:-${SCRIPT_DIR}/R/eval_test.R}"

source "${BACDIVE_CONDA_SH}"
conda activate "${BACDIVE_TF_ENV}"

echo "=== BacDive Test Set Evaluation ==="
echo "Node: $(hostname)"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null)"
echo "Date: $(date)"
echo ""

Rscript "${EVAL_SCRIPT_PATH}"

echo ""
echo "=== Done ==="
echo "Date: $(date)"
