#!/bin/sh
#SBATCH --job-name=rm_all7
#SBATCH --time=08:00:00
#SBATCH --mem=20Gb
#SBATCH --cpus-per-task=1
#SBATCH --partition=workq
#SBATCH -e rm_all7.%A.%j.%a.err
#SBATCH -o rm_all7.%A.%j.%a.out
#SBATCH --array=1-6

seednum=$1
cleanf=$2
unreviewedf=$3
model=$4
# meth=$2
# Correct floating-point arithmetic using awk
# PRED=$(awk "BEGIN {print 0.05 * $SLURM_ARRAY_TASK_ID}")
# python 01remove_rxn.py -- --prec $PRED --seednum $seednum --meth $meth


# Run the Python script with the calculated precision value
echo 'running for original clean df'
# echo seednum $seednum
# python test.py --prec $PRED --seednum $seednum
# python 01remove_rxn.py --prec 0.05  --seednum $seednum --meth $SLURM_ARRAY_TASK_ID --cleanf  $cleanf --unreviewedf $unreviewedf --model $model --ori $5
python 01remove_rxn.py --prec 0.1  --seednum $seednum --meth $SLURM_ARRAY_TASK_ID --cleanf  $cleanf --unreviewedf $unreviewedf --model $model --ori $5
# python 01remove_rxn.py --prec 0.15  --seednum $seednum --meth $SLURM_ARRAY_TASK_ID --cleanf  $cleanf --unreviewedf $unreviewedf --model $model --ori $5
python 01remove_rxn.py --prec 0.20  --seednum $seednum --meth $SLURM_ARRAY_TASK_ID --cleanf  $cleanf --unreviewedf $unreviewedf --model $model --ori $5
# python 01remove_rxn.py --prec 0.25  --seednum $seednum --meth $SLURM_ARRAY_TASK_ID --cleanf  $cleanf --unreviewedf $unreviewedf --model $model --ori $5
python 01remove_rxn.py --prec 0.30  --seednum $seednum --meth $SLURM_ARRAY_TASK_ID --cleanf  $cleanf --unreviewedf $unreviewedf --model $model --ori $5
# python 01remove_rxn.py --prec 0.35  --seednum $seednum --meth $SLURM_ARRAY_TASK_ID --cleanf  $cleanf --unreviewedf $unreviewedf --model $model --ori $5
python 01remove_rxn.py --prec 0.40  --seednum $seednum --meth $SLURM_ARRAY_TASK_ID --cleanf  $cleanf --unreviewedf $unreviewedf --model $model --ori $5
# python 01remove_rxn.py --prec 0.45  --seednum $seednum --meth $SLURM_ARRAY_TASK_ID --cleanf  $cleanf --unreviewedf $unreviewedf --model $model --ori $5
python 01remove_rxn.py --prec 0.50  --seednum $seednum --meth $SLURM_ARRAY_TASK_ID --cleanf  $cleanf --unreviewedf $unreviewedf --model $model --ori $5


# python test.py --prec 0.05
# python test.py --prec 0.1
# python test.py --prec 0.15
# python test.py --prec 0.2
# python test.py --prec 0.25
# python test.py --prec 0.3
# python test.py --prec 0.35
# python test.py --prec 0.4
# python test.py --prec 0.45
# python test.py --prec 0.5
