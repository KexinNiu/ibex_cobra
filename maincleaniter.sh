#!/bin/sh
#SBATCH --job-name=maiter_R_tainn
#SBATCH --time=10:50:00
#SBATCH --mem=200Gb
#SBATCH --cpus-per-task=8
#SBATCH --partition=workq
#SBATCH -e maiter_R_tainn.%J.err
#SBATCH -o maiter_R_tainn.%J.out


python mainclean_interation_v1ecoli.py --input_file $1 --cleanfile $2 --file_type 5 --reward 1 --iter 10 --cpu 4 --gram $3 --name $4 --threshold 8 --media $5
python mainclean_interation_v2.py \
    --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1_t1.fasta \
    --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl \
    --file_type 7 --reward 0.1 --iter 10 --cpu 8 --gram negative  --block_flage 0 --flux_flage 0 --name iLJ478_allec
# media="'cpd00001_e','cpd00007_e','cpd00009_e','cpd00013_e','cpd00018_e','cpd00023_e','cpd00027_e','cpd00028_e','cpd00030_e','cpd00033_e','cpd00034_e','cpd00035_e','cpd00039_e','cpd00041_e','cpd00046_e','cpd00048_e','cpd00051_e','cpd00054_e','cpd00058_e','cpd00060_e','cpd00063_e','cpd00065_e','cpd00066_e','cpd00067_e','cpd00069_e','cpd00091_e','cpd00092_e','cpd00099_e','cpd00107_e','cpd00119_e','cpd00126_e','cpd00129_e','cpd00149_e','cpd00156_e','cpd00161_e','cpd00182_e','cpd00184_e','cpd00205_e','cpd00215_e','cpd00218_e','cpd00220_e','cpd00226_e','cpd00239_e','cpd00244_e','cpd00246_e','cpd00249_e','cpd00254_e','cpd00311_e','cpd00322_e','cpd00381_e','cpd00393_e','cpd00438_e','cpd00531_e','cpd00541_e','cpd00644_e','cpd00654_e','cpd00793_e','cpd00971_e','cpd01012_e','cpd01048_e','cpd03424_e','cpd10515_e','cpd10516_e','cpd11574_e','cpd09695_e'"
# echo media: $media
# python mainclean_interation_v1ecoli.py --input_file $1 --cleanfile $2 --file_type 5 --reward 1 --iter 10 --cpu 4 --gram $3 --name $4 --threshold 0.8 --media $media



# python mainclean_interation_v1ecoli.py --input_file $1 --cleanfile $2 --file_type 5 --reward 0.1 --iter 10 --cpu 4 --gram $3 --name $4 --threshold 0.3
# python mainclean_interation_v1ecoli.py --input_file $1 --cleanfile $2 --file_type 5 --reward 0.1 --iter 10 --cpu 4 --gram $3 --name $4 --threshold 0.5

# python mainclean_interation_v1ecoli_conti.py --input_file $1 --cleanfile $2 --file_type 5 --reward 0.1 --iter 10 --cpu 4 --gram $3 --name $4 --threshold 0.3 --conti $5
# python mainclean_interation_v1ecoli_conti.py --input_file $1 --cleanfile $2 --file_type 5 --reward 0.1 --iter 10 --cpu 4 --gram $3 --name $4 --threshold 0.5 --conti $5

# echo try for iaf987 reverse
# python mainclean_interation_v0.py --input_file /ibex/user/niuk0a/funcarve/cobra/aaseq_CP000148.1.fasta --cleanfile /ibex/user/niuk0a/funcarve/cobra/aaseq_CP000148.1_maxsep.csv  --file_type 5 --reward 0.1 --iter 10 --cpu 8 --gram negative --name iaf987INTER_RR5
                              



# echo for ecoli  reverse
# python mainclean_interation_v1ecoli.py --input_file /ibex/user/niuk0a/funcarve/cobra/NC_000913.3_ecoli.fasta  --cleanfile /ibex/user/niuk0a/funcarve/cobra/NC_000913.3_ecoli_maxsep.csv  --file_type 5 --reward 0.1 --iter 10 --cpu 4 --gram negative --name ecoliINTER1_RR5

# echo for 

# iLJ478	NC_000853.1	243274	Thermotoga maritima MSB8		243274 negative
# iJN678	BA000022.2	1148	Synechocystis sp. PCC 6803		1148    negative
# iJN1463	NC_002947.4	160488	Pseudomonas putida KT2440		160488  negative
# iCN900	AM180355.1	272563	Clostridioides difficile 630		272563  postive
# iHN637	NC_014328.1	748727	Clostridium ljungdahlii DSM 13528		748727  postive
# iAF1260	NC_000913.3	511145/83333	Escherichia coli str. K-12 substr. MG1655		511145  negative
# iAF987	CP000148.1	269799	Geobacter metallireducens GS-15		269799  negative
# iCN718	NC_010410.1	509173	Acinetobacter baumannii AYE		509173  negative

# NOT USEFUL FROM NOW ON echo for iLJ478 reverse
# python mainclean_interation_v0.py --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1.fasta --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_maxsep.csv  --file_type 5 --reward 0.1 --iter 10 --cpu 8 --gram negative --name iLJ478_R


# echo for iJN678 reverse
# python mainclean_interation_v0.py --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/BA000022.2.fasta --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/BA000022.2_maxsep.csv  --file_type 5 --reward 0.1 --iter 10 --cpu 8 --gram negative --name iJN678_R_test
# echo for iJN1463 reverse
# python mainclean_interation_v0.py --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_002947.4.fasta --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_maxsep.csv  --file_type 5 --reward 0.1 --iter 10 --cpu 8 --gram negative --name iJN1463_R
# echo for iCN900 reverse
# python mainclean_interation_v0.py --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/AM180355.1.fasta --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_maxsep.csv  --file_type 5 --reward 0.1 --iter 10 --cpu 8 --gram positive --name iCN900_R
# echo for iHN637 reverse
# python mainclean_interation_v0.py --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_014328.1.fasta --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_014328.1_maxsep.csv  --file_type 5 --reward 0.1 --iter 10 --cpu 8 --gram positive --name iHN637_R
# echo for iAF1260 reverse
# python mainclean_interation_v0.py --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000913.3.fasta --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000913.3_maxsep.csv  --file_type 5 --reward 0.1 --iter 10 --cpu 8 --gram negative --name iAF1260_R

# echo for iCN718 reverse
# python mainclean_interation_v0.py --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_010410.1.fasta --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_010410.1_maxsep.csv  --file_type 5 --reward 0.1 --iter 10 --cpu 8 --gram negative --name iCN718_R
# echo done all reverse

# wait for clean result finished run 
# 9	iNJ661	NC_000962.3	83332	Mycobacterium tuberculosis H37Rv		83332	pos
# echo for inj661 reverse
# python mainclean_interation_v0.py --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000962.3.fasta --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3_maxsep.csv  --file_type 5 --reward 0.1 --iter 10 --cpu 8 --gram positive --name iNJ661_R

# echo for iAF987 reverse
# python mainclean_interation_v0.py --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/CP000148.1.fasta --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1_maxsep.csv  --file_type 5 --reward 0.1 --iter 10 --cpu 8 --gram negative --name iAF987_R

## allec reward =0.1
# 		genome 	ORGID	ORGNASIM	
# t1	iLJ478	NC_000853.1	243274	Thermotoga maritima MSB8	neg
# t2	iJN1463	NC_002947.4	160488	Pseudomonas putida KT2440	neg
# t2	iCN900	AM180355.1	272563	Clostridioides difficile 630	pos
# t4	iAF1260	NC_000913.3	83333	Escherichia coli str. K-12 substr. MG1655	neg
# t4	iAF987	CP000148.1	269799	Geobacter metallireducens GS-15	neg
# t5	iNJ661	NC_000962.3	83332	Mycobacterium tuberculosis H37Rv	pos
# t4	iYO844	AL009126.3	224308	Bacillus subtilis subsp. subtilis str. 168	pos
python mainclean_interation_v2.py \
    --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1_t1.fasta \
    --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl \
    --file_type 7 --reward 0.1 --iter 10 --cpu 8 --gram negative  --block_flage 0 --flux_flage 0 --name iLJ478_allec

python mainclean_interation_v2.py \
    --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_002947.4_t2.fasta \
    --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_t2_maxsep_df.pkl \
    --file_type 7 --reward 0.1 --iter 10 --cpu 8 --gram negative  --block_flage 0 --flux_flage 0 --name iJN1463_allec

# python mainclean_interation_v2.py \
#     --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/AM180355.1_t2.fasta \
#     --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_t2_maxsep_df.pkl \
#     --file_type 7 --reward 0.1 --iter 10 --cpu 8 --gram positive  --block_flage 0 --flux_flage 0 --name iCN900_allec

# python mainclean_interation_v2.py \
#     --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000913.3_t4.fasta \
#     --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000913.3_t4_maxsep_df.pkl \
#     --file_type 7 --reward 0.1 --iter 10 --cpu 8 --gram negative  --block_flage 0 --flux_flage 0 --name iAF1260_allec

# python mainclean_interation_v2.py \
#     --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/CP000148.1_t4.fasta \
#     --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1_t4_maxsep_df.pkl \
#     --file_type 7 --reward 0.1 --iter 10 --cpu 8 --gram negative  --block_flage 0 --flux_flage 0 --name iAF987_allec

python mainclean_interation_v2.py \
    --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000962.3_t5.fasta \
    --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3_t5_maxsep_df.pkl \
    --file_type 7 --reward 0.1 --iter 10 --cpu 8 --gram positive  --block_flage 0 --flux_flage 0 --name iNJ661_allec

python mainclean_interation_v2.py \
    --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/AL009126.3_t4.fasta \
    --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/AL009126.3_t4_maxsep_df.pkl \
    --file_type 7 --reward 0.1 --iter 10 --cpu 8 --gram positive  --block_flage 0 --flux_flage 0 --name iYO844_allec
# AL009126.3_t4_maxsep_df.pkl
# AM180355.1_t2_maxsep_df.pkl
# CP000148.1_t4_maxsep_df.pkl
# NC_000853.1_t1_maxsep_df.pkl
# NC_000913.3_t4_maxsep_df.pkl
# NC_000962.3_t5_maxsep_df.pkl
# NC_002947.4_t2_maxsep_df.pkl
# AL009126.3_t4.fasta
# AM180355.1_t2.fasta
# CP000148.1_t4.fasta
# NC_000853.1_t1.fasta
# NC_000913.3_t4.fasta
# NC_000962.3_t5.fasta
# NC_002947.4_t2.fasta