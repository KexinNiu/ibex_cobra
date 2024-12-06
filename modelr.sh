#!/bin/sh
#SBATCH --job-name=comv0
#SBATCH --time=02:00:00
#SBATCH --mem=20Gb
#SBATCH --cpus-per-task=1
#SBATCH --partition=workq
#SBATCH -e modelcomv0.%J.err
#SBATCH -o modelcomv0.%J.out


echo "v0"
python model_compare_v0_final.py --prec 0.1 --seednum 555 --gf w
python model_compare_v0_final.py --prec 0.1 --seednum 555 --gf ori
# echo "v1"
# # python model_compare_v1.py --prec 0.1 --seednum 555 --model /ibex/user/niuk0a/funcarve/cobra/iAF987.xml --type gfw --cleanfile /ibex/user/niuk0a/funcarve/cobra/aaseq_CP000148.1_maxsep.csv
# # python model_compare_v1.py --prec 0.25 --seednum 555 --model /ibex/user/niuk0a/funcarve/cobra/iAF987.xml --type gf 
# echo model iYL1228
# python model_compare_v1.py --prec 0.1 --seednum 555 --model /ibex/user/niuk0a/funcarve/cobra/iYL1228.xml --type gfw --cleanfile /ibex/user/niuk0a/funcarve/cobra/CP000647.1_sequence_iyj1228.csv