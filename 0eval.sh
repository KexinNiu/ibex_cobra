# t1	iLJ478	NC_000853.1	243274
# t2	iJN1463	NC_002947.4	160488
# t2	iCN900	AM180355.1	272563
# t4	iAF1260	NC_000913.3	83333
# t4	iAF987	CP000148.1	269799
# t5	iNJ661	NC_000962.3	83332
# t4	iYO844	AL009126.3	224308
# sbatch 0batch.sh /ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_243274.tsv /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1  iLJ478 1
# sbatch 0batch.sh /ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_160488.tsv /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4 iJN1463 2
# sbatch 0batch.sh /ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_272563.tsv /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1 iCN900 2
# sbatch 0batch.sh /ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_83333.tsv /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000913.3 iAF1260 4
# sbatch 0batch.sh /ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_269799.tsv /ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1 iAF987 4
# sbatch 0batch.sh /ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_83332.tsv /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3 iNJ661 5
# sbatch 0batch.sh /ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_224308.tsv /ibex/user/niuk0a/CLEAN/app/results/inputs/AL009126.3 iYO844 4

## only ori 
sbatch 0batch.sh /ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_224308.tsv /ibex/user/niuk0a/CLEAN/app/results/inputs/AL009126.3 iYO844 4
sbatch 0batch.sh /ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_272563.tsv /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1 iCN900 2
sbatch 0batch.sh /ibex/user/niuk0a/funcarve/cobra/uniprot/uniprot_83333.tsv /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000913.3 iAF1260 4

