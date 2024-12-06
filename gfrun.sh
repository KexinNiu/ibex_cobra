#!/bin/sh
#SBATCH --job-name=gftest 
#SBATCH --time=00:20:00
#SBATCH --mem=20Gb
#SBATCH --cpus-per-task=1
#SBATCH --partition=workq
#SBATCH -e gftest_1228.%A.%j.%a.err
#SBATCH -o gftest_1228.%A.%j.%a.out



# Correct floating-point arithmetic using awk
# PRED=$(awk "BEGIN {print 0.05 * $SLURM_ARRAY_TASK_ID}")

# # Run the Python script with the calculated precision value
# seednum=$1
# echo seednum $seednum
# python test.py --prec $PRED --seednum $seednum
# python gapfill.py -dm -type gf 
# echo "#######################"
# echo "Running gapfill.py type gf"
# python gapfill.py -dm /ibex/user/niuk0a/funcarve/cobra/iAF987_555_0.05_rm.xml --type gf -cf /ibex/user/niuk0a/funcarve/cobra/aaseq_CP000148.1_maxsep.csv
echo "#######################"
# echo "Running gapfill.py type gfw"
# python gapfill.py -dm /ibex/user/niuk0a/funcarve/cobra/iAF987_555_0.05_rm.xml --type gf -cf /ibex/user/niuk0a/funcarve/cobra/aaseq_CP000148.1_maxsep.csv
python iyl1228.py
echo "#######################"


