#!/bin/sh
#SBATCH --job-name=iyl1228noai
#SBATCH --time=10:00:00
#SBATCH --mem=50Gb
#SBATCH --cpus-per-task=8
#SBATCH --partition=workq
#SBATCH -e iyl1228noai.%J.err
#SBATCH -o iyl1228noai.%J.out


python iyl1228.py