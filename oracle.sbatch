#!/bin/bash
#SBATCH --job-name=oracle
#SBATCH --partition=SallesInfo
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
#SBATCH --time=02:00:00
#SBATCH --array=1-45
#SBATCH --output=oracle_%A_%a.out
# Self-locating: works from any checkout. This script lives in <ROOT>/bayesian_oracle/.
# Config knobs (paths, MAX_2B_N, ...) live in bayesian_oracle/config.py and are env-overridable;
# the array size must match the number of tumours with N<=MAX_2B_N (45 in the paper run).
HERE="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
ROOT="$(dirname "$HERE")"
cd "$ROOT"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-2} MAX_2B_N=${MAX_2B_N:-20000}
mkdir -p logs "${ORACLE_OUT:-scoring/oracle_scores}"
echo "### task=$SLURM_ARRAY_TASK_ID node=$(hostname) root=$ROOT start=$(date)"
python3 bayesian_oracle/oracle_score_one.py
echo "### task=$SLURM_ARRAY_TASK_ID DONE=$(date)"
