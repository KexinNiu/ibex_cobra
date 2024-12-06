#!/bin/sh
#SBATCH --job-name=bigg_meta
#SBATCH --time=08:00:00
#SBATCH --mem=20Gb
#SBATCH --cpus-per-task=1
#SBATCH --partition=workq
#SBATCH -e getbigg_$SLURM_ARRAY_TASK_ID.%J.err
#SBATCH -o getbigg_$SLURM_ARRAY_TASK_ID.%J.out
#SBATCH --array=0-3

python get_bigg.py --part $SLURM_ARRAY_TASK_ID