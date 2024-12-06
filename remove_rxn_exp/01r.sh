#!/bin/sh
#SBATCH --job-name=testgf
#SBATCH --time=00:30:00
#SBATCH --mem=20Gb
#SBATCH --cpus-per-task=1
#SBATCH --partition=workq
#SBATCH -e test_se.%j.err
#SBATCH -o test_se.%j.out


# Correct floating-point arithmetic using awk
# PRED=$(awk "BEGIN {print 0.05 * $SLURM_ARRAY_TASK_ID}")

# Run the Python script with the calculated precision value
seednum=$1
# echo seednum $seednum
# python test.py --prec $PRED --seednum $seednum
python 01remove_rxn.py --prec 0.05  --seednum $seednum
python 01remove_rxn.py --prec 0.1  --seednum $seednum
python 01remove_rxn.py --prec 0.15  --seednum $seednum
python 01remove_rxn.py --prec 0.2  --seednum $seednum
python 01remove_rxn.py --prec 0.25  --seednum $seednum
python 01remove_rxn.py --prec 0.3  --seednum $seednum
python 01remove_rxn.py --prec 0.35  --seednum $seednum
python 01remove_rxn.py --prec 0.4  --seednum $seednum
python 01remove_rxn.py --prec 0.45  --seednum $seednum
python 01remove_rxn.py --prec 0.5  --seednum $seednum

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
