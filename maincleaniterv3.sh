#!/bin/sh
#SBATCH --job-name=v3_t7rewards
#SBATCH --time=10:50:00
#SBATCH --mem=40Gb
#SBATCH --cpus-per-task=1
#SBATCH --partition=workq
#SBATCH -e v3_t7rewards.%J.err
#SBATCH -o v3_t7rewards.%J.out

# echo "python mainclean_interation_v2.py --input_file $1 --cleanfile $2 --file_type 7 --reward $8 --iter 10 --cpu 8 --gram $3 --block_flage $4 --flux_flage $5 --name $6 --threshold 8 --media $7"
# python mainclean_interation_v2.py --input_file $1 --cleanfile $2 --file_type 7 --reward $8 --iter 4 --cpu 8 --gram $3 --block_flage $4 --flux_flage $5 --name $6 --threshold 8 --media $7
echo "python mainclean_interation_v3.py --input_file $1 --cleanfile $2 --file_type 7 --reward $8 --iter 10 --cpu 8 --gram $3 --block_flage $4 --flux_flage $5 --name $6 --threshold 8 --media $7"
python mainclean_interation_v3.py --input_file $1 --cleanfile $2 --file_type 7 --reward $8 --iter 4 --cpu 8 --gram $3 --block_flage $4 --flux_flage $5 --name $6 --threshold 6 --media $7
## allec reward =0.1
# 		genome 	ORGID	ORGNASIM	
# t1	iLJ478	NC_000853.1	243274	Thermotoga maritima MSB8	neg
# t2	iJN1463	NC_002947.4	160488	Pseudomonas putida KT2440	neg
# t2	iCN900	AM180355.1	272563	Clostridioides difficile 630	pos
# t4	iAF1260	NC_000913.3	83333	Escherichia coli str. K-12 substr. MG1655	neg
# t4	iAF987	CP000148.1	269799	Geobacter metallireducens GS-15	neg
# t5	iNJ661	NC_000962.3	83332	Mycobacterium tuberculosis H37Rv	pos
# t4	iYO844	AL009126.3	224308	Bacillus subtilis subsp. subtilis str. 168	pos
