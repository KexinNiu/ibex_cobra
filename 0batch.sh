#!/bin/sh
#SBATCH --job-name=0_batch_eval
#SBATCH --time=00:10:00
#SBATCH --mem=5Gb
#SBATCH --cpus-per-task=1
#SBATCH --partition=workq
#SBATCH -e 0_batch_eval.%J.err
#SBATCH -o 0_batch_eval.%J.out

# sbatch 0batch.sh 00_eval.py /ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_243274.tsv /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1  iLJ478 1
# python 00_eval.py --file_name /ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_243274.tsv --pref /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1  --modelname iLJ478 --index 1
python 00_eval.py --file_name $1 --pref $2 --modelname $3 --index $4