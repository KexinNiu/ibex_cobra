# python mainclean_interation_v0.py --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1.fasta --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_maxsep.csv  --file_type 5 --reward 0.1 --iter 10 --cpu 8 --gram negative --name iLJ478_R
# echo for iJN678 reverse
# python mainclean_interation_v0.py --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/BA000022.2.fasta --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/BA000022.2_maxsep.csv  --file_type 5 --reward 0.1 --iter 10 --cpu 8 --gram negative --name iJN678_R
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
# 0.3
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_maxsep.csv negative iLJ478_R 6
# # sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/BA000022.2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/BA000022.2_maxsep.csv negative iJN678_R no
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_002947.4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_maxsep.csv negative iJN1463_R 7
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AM180355.1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_maxsep.csv positive iCN900_R  7
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_014328.1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_014328.1_maxsep.csv positive iHN637_R 8
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000913.3.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000913.3_maxsep.csv negative iAF1260_R 6
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_010410.1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_010410.1_maxsep.csv negative iCN718_R 7
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000962.3.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3_maxsep.csv positive iNJ661_R 7
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/CP000148.1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1_maxsep.csv negative iAF987_R 7

# 0.5
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_maxsep.csv negative iLJ478_R 7
# # sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/BA000022.2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/BA000022.2_maxsep.csv negative iJN678_R no
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_002947.4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_maxsep.csv negative iJN1463_R 8
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AM180355.1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_maxsep.csv positive iCN900_R  9
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_014328.1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_014328.1_maxsep.csv positive iHN637_R 9
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000913.3.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000913.3_maxsep.csv negative iAF1260_R 7
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_010410.1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_010410.1_maxsep.csv negative iCN718_R 8
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000962.3.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3_maxsep.csv positive iNJ661_R 8
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/CP000148.1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1_maxsep.csv negative iAF987_R 8


# /ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1_t4_maxsep.csv
# reward=1 threshold=8
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/CP000148.1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1_t4_maxsep.csv negative iAF987_Rt4 
# python CLEAN_infer_fasta_train.py --fasta_data NC000913.3_t4 --traindata split70_bactaxo_train4 
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC000913.3.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC000913.3_t4_maxsep.csv negative iAF1260_Rt4
# python CLEAN_infer_fasta_train.py --fasta_data AL009126.3_t4 --traindata split70_bactaxo_train4 
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AL009126.3.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AL009126.3_t4_maxsep.csv negative iYO844_Rt4
# # NC_000853.1_t1
# python CLEAN_infer_fasta_train.py --fasta_data NC_000853.1_t1 --traindata split70_bactaxo_train1
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep.csv negative iLJ478_Rt1
# # # NC_000962.3_t5
# # python CLEAN_infer_fasta_train.py --fasta_data NC_000962.3_t5 --traindata split70_bactaxo_train5
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000962.3.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3_t5_maxsep.csv positive iNJ661_Rt5
# # # NC_002947.4_t2
# # python CLEAN_infer_fasta_train.py --fasta_data NC_002947.4_t2 --traindata split70_bactaxo_train2
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_002947.4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_t2_maxsep.csv negative iJN1463_Rt2

# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000913.3.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000913.3_t4_maxsep.csv negative iAF1260_Rt4
# AM180355.1
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AM180355.1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_t2_maxsep.csv positive iCN900_Rt2

# # /ibex/user/niuk0a/CLEAN/app/results/inputs/S1_5_aminoacid_maxsep.csv
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/S1_5_aminoacid.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/S1_5_aminoacid_maxsep.csv negative S1_5_aminoacid_Rt4
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/L11_genomwhole.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/L11_genomwhole_maxsep.csv negative L11_genomwhole_Rt4
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/tiannv/sample1_proteins.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/sample1_proteins_maxsep.csv negative sample1_proteins_Rt4
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/tiannv/sample11_proteins.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/sample11_proteins_maxsep.csv negative sample11_proteins
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/tiannv/sample12_proteins.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/sample12_proteins_maxsep.csv negative sample12_proteins
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/tiannv/sample14_proteins.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/sample14_proteins_maxsep.csv positive sample14_proteins
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/tiannv/sample31_proteins.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/sample31_proteins_maxsep.csv positive sample13_proteins

# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/krish/D21_C_firmus_B_proteins.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/D21_C_firmus_B_proteins_maxsep.csv positive D21_C_firmus_B_proteins
# # /ibex/user/niuk0a/CLEAN/app/results/inputs/D21_C_israelensis_proteins_maxsep.csv
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/krish/D21_C_israelensis_proteins.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/D21_C_israelensis_proteins_maxsep.csv positive D21_C_israelensis_proteins
# # /ibex/user/niuk0a/CLEAN/app/results/inputs/D21_H_caseinilytica_proteins_maxsep.csv
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/krish/D21_H_caseinilytica_proteins.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/D21_H_caseinilytica_proteins_maxsep.csv negative D21_H_caseinilytica_proteins


# /ibex/user/niuk0a/CLEAN/app/results/inputs/UP000009073_595494_maxsep.csv 
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/UP000009073_595494.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/UP000009073_595494_maxsep.csv positive UP000009073_595494
# UP000000238_349521.fasta
# UP000000248_246195.fasta
# UP000000268_329726.fasta
# UP000000269_350688.fasta
# UP000000366_420662.fasta
# UP000000422_273121.fasta
# UP000000432_220668.fasta
# UP000000528_272635.fasta
# UP000000556_160488.fasta
# UP000000577_243231.fasta
# UP000000636_292414.fasta
# UP000000640_196162.fasta
# UP000000657_326424.fasta
# UP000000719_373903.fasta
# UP000000788_93059.fasta
# UP000000813_176299.fasta
# UP000001007_194439.fasta
# UP000001010_190485.fasta
# UP000001114_339671.fasta
# UP000001168_66692.fasta
# UP000001208_517418.fasta
# UP000001349_394503.fasta
# UP000001362_243159.fasta
# UP000001383_458233.fasta
# UP000001418_710127.fasta
# UP000001938_321332.fasta
# UP000001946_138119.fasta
# UP000001964_394221.fasta
# UP000002156_340099.fasta
# UP000002402_246197.fasta
# UP000002495_235279.fasta
# UP000002601_526222.fasta
# UP000002671_227377.fasta
# UP000002707_314315.fasta
# UP000002709_319225.fasta
# UP000006637_266117.fasta
# UP000006647_265311.fasta
# UP000006820_247156.fasta
# UP000007966_218491.fasta
# UP000008080_264462.fasta
# UP000008561_96561.fasta
# UP000008701_290317.fasta
# UP000008808_314225.fasta
# UP000008871_393595.fasta
# UP000009073_595494.fasta


## version 2 for test genomes
# 		genome 	ORGID	ORGNASIM	
# t1	iLJ478	NC_000853.1	243274	Thermotoga maritima MSB8	neg
# t2	iJN1463	NC_002947.4	160488	Pseudomonas putida KT2440	neg
# t2	iCN900	AM180355.1	272563	Clostridioides difficile 630	pos
# t4	iAF1260	NC_000913.3	83333	Escherichia coli str. K-12 substr. MG1655	neg
# t4	iAF987	CP000148.1	269799	Geobacter metallireducens GS-15	neg
# t5	iNJ661	NC_000962.3	83332	Mycobacterium tuberculosis H37Rv	pos
# t4	iYO844	AL009126.3	224308	Bacillus subtilis subsp. subtilis str. 168	pos 

# all ec reward =0.1

sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1_t1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl negative 0 0 iLJ478_allec default
sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_002947.4_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_t2_maxsep_df.pkl negative 0 0 iJN1463_allec default
sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AM180355.1_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_t2_maxsep_df.pkl positive 0 0 iCN900_allec default
sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000913.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000913.3_t4_maxsep_df.pkl negative 0 0 iAF1260_allec default
sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/CP000148.1_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1_t4_maxsep_df.pkl negative 0 0 iAF987_allec default
sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000962.3_t5.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3_t5_maxsep_df.pkl positive 0 0 iNJ661_allec default
sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AL009126.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AL009126.3_t4_maxsep_df.pkl positive 0 0 iYO844_allec default

##  block only 
sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1_t1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl negative 1 0 iLJ478_block default
sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_002947.4_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_t2_maxsep_df.pkl negative 1 0 iJN1463_block default
sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AM180355.1_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_t2_maxsep_df.pkl positive 1 0 iCN900_block default
sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000913.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000913.3_t4_maxsep_df.pkl negative 1 0 iAF1260_block default
sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/CP000148.1_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1_t4_maxsep_df.pkl negative 1 0 iAF987_block default
sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000962.3_t5.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3_t5_maxsep_df.pkl positive 1 0 iNJ661_block default
sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AL009126.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AL009126.3_t4_maxsep_df.pkl positive 1 0 iYO844_block default

## flux only
sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1_t1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl negative 0 1 iLJ478_flux default
sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_002947.4_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_t2_maxsep_df.pkl negative 0 1 iJN1463_flux default
sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AM180355.1_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_t2_maxsep_df.pkl positive 0 1 iCN900_flux default
sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000913.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000913.3_t4_maxsep_df.pkl negative 0 1 iAF1260_flux default
sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/CP000148.1_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1_t4_maxsep_df.pkl negative 0 1 iAF987_flux default
sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000962.3_t5.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3_t5_maxsep_df.pkl positive 0 1 iNJ661_flux default
sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AL009126.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AL009126.3_t4_maxsep_df.pkl positive 0 1 iYO844_flux default

## block and flux
sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1_t1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl negative 1 1 iLJ478_fluxblock default
sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_002947.4_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_t2_maxsep_df.pkl negative 1 1 iJN1463_fluxblock default
sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AM180355.1_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_t2_maxsep_df.pkl positive 1 1 iCN900_fluxblock default
sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000913.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000913.3_t4_maxsep_df.pkl negative 1 1 iAF1260_fluxblock default
sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/CP000148.1_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1_t4_maxsep_df.pkl negative 1 1 iAF987_fluxblock default
sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000962.3_t5.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3_t5_maxsep_df.pkl positive 1 1 iNJ661_fluxblock default
sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AL009126.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AL009126.3_t4_maxsep_df.pkl positive 1 1 iYO844_fluxblock default

# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/S1_5_aminoacid.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/S1_5_aminoacid_maxsep.csv negative S1_5_aminoacid_Rt4
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/L11_genomwhole.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/L11_genomwhole_maxsep.csv negative L11_genomwhole_Rt4
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/tiannv/sample1_proteins.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/sample1_proteins_maxsep.csv negative sample1_proteins_Rt4
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/tiannv/sample11_proteins.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/sample11_proteins_maxsep.csv negative sample11_proteins
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/tiannv/sample12_proteins.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/sample12_proteins_maxsep.csv negative sample12_proteins
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/tiannv/sample14_proteins.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/sample14_proteins_maxsep.csv positive sample14_proteins
# sbatch maincleaniter.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/tiannv/sample31_proteins.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/sample31_proteins_maxsep.csv positive sample13_proteins

sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/tiannv/sample1_proteins.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/tiannv/sample1_proteins_maxsep_df.pkl negative 0 0 sample1_proteins ????media