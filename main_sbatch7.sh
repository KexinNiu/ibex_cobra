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

# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1_t1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl negative 0 0 iLJ478_allec default
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_002947.4_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_t2_maxsep_df.pkl negative 0 0 iJN1463_allec default
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AM180355.1_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_t2_maxsep_df.pkl positive 0 0 iCN900_allec default
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000913.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000913.3_t4_maxsep_df.pkl negative 0 0 iAF1260_allec default
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/CP000148.1_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1_t4_maxsep_df.pkl negative 0 0 iAF987_allec default
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000962.3_t5.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3_t5_maxsep_df.pkl positive 0 0 iNJ661_allec default
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AL009126.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AL009126.3_t4_maxsep_df.pkl positive 0 0 iYO844_allec default

# ##  block only 
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1_t1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl negative 1 0 iLJ478_block_$rname default
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_002947.4_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_t2_maxsep_df.pkl negative 1 0 iJN1463_block_$rname default
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AM180355.1_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_t2_maxsep_df.pkl positive 1 0 iCN900_block_$rname default
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000913.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000913.3_t4_maxsep_df.pkl negative 1 0 iAF1260_block_$rname default
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/CP000148.1_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1_t4_maxsep_df.pkl negative 1 0 iAF987_block_$rname default
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000962.3_t5.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3_t5_maxsep_df.pkl positive 1 0 iNJ661_block_$rname default
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AL009126.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AL009126.3_t4_maxsep_df.pkl positive 1 0 iYO844_block_$rname default

# ## flux only
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1_t1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl negative 0 1 iLJ478_flux_$rname default
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_002947.4_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_t2_maxsep_df.pkl negative 0 1 iJN1463_flux_$rname default
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AM180355.1_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_t2_maxsep_df.pkl positive 0 1 iCN900_flux_$rname default
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000913.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000913.3_t4_maxsep_df.pkl negative 0 1 iAF1260_flux_$rname default
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/CP000148.1_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1_t4_maxsep_df.pkl negative 0 1 iAF987_flux_$rname default
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000962.3_t5.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3_t5_maxsep_df.pkl positive 0 1 iNJ661_flux_$rname default
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AL009126.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AL009126.3_t4_maxsep_df.pkl positive 0 1 iYO844_flux_$rname default

# ## block and flux
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1_t1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl negative 1 1 iLJ478_fluxblock_$rname  default
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_002947.4_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_t2_maxsep_df.pkl negative 1 1 iJN1463_fluxblock_$rname  default
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AM180355.1_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_t2_maxsep_df.pkl positive 1 1 iCN900_fluxblock_$rname  default
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000913.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000913.3_t4_maxsep_df.pkl negative 1 1 iAF1260_fluxblock_$rname  default
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/CP000148.1_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1_t4_maxsep_df.pkl negative 1 1 iAF987_fluxblock_$rname  default
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000962.3_t5.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3_t5_maxsep_df.pkl positive 1 1 iNJ661_fluxblock_$rname  default
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AL009126.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AL009126.3_t4_maxsep_df.pkl positive 1 1 iYO844_fluxblock_$rname  default


## test reward scale from 0.1 to 1.0

# Test reward scale from 0.1 to 1.0
echo "Starting reward tests"
# rewards=(0.1 0.3 0.5 0.7 1.0)
# rewards=(0.9)
# rewards=(0.05 1.5 2.0 3.0 )
# rewards=(0.05 0.2 5.0)
typrds=(allec block flux fluxblock)

# Iterate over reward scales
for reward in "${rewards[@]}"; do
    if (( $(echo "$reward < 1.0" | bc -l) )); then
        rname="p$(echo "$reward" | sed 's/\.//g')"
    else
        rname=$(printf "%.0f" "$reward")
    fi

    # Iterate over types
    for typrd in "${typrds[@]}"; do
        if [[ "$typrd" == "allec" ]]; then
            sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1_t1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl negative 0 0 iLJ478_allec_"$rname" default "$reward"
            sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_002947.4_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_t2_maxsep_df.pkl negative 0 0 iJN1463_allec_"$rname" default "$reward"
            sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AM180355.1_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_t2_maxsep_df.pkl positive 0 0 iCN900_allec_"$rname" default "$reward"
            sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000913.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000913.3_t4_maxsep_df.pkl negative 0 0 iAF1260_allec_"$rname" default "$reward"
            sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/CP000148.1_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1_t4_maxsep_df.pkl negative 0 0 iAF987_allec_"$rname" default "$reward"
            sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000962.3_t5.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3_t5_maxsep_df.pkl positive 0 0 iNJ661_allec_"$rname" default "$reward"
            sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AL009126.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AL009126.3_t4_maxsep_df.pkl positive 0 0 iYO844_allec_"$rname" default "$reward"
        elif [[ "$typrd" == "block" ]]; then
            sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1_t1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl negative 1 0 iLJ478_block_"$rname" default "$reward"
            sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_002947.4_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_t2_maxsep_df.pkl negative 1 0 iJN1463_block_"$rname" default "$reward"
            sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AM180355.1_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_t2_maxsep_df.pkl positive 1 0 iCN900_block_"$rname" default "$reward"
            sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000913.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000913.3_t4_maxsep_df.pkl negative 1 0 iAF1260_block_"$rname" default "$reward"
            sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/CP000148.1_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1_t4_maxsep_df.pkl negative 1 0 iAF987_block_"$rname" default "$reward"
            sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000962.3_t5.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3_t5_maxsep_df.pkl positive 1 0 iNJ661_block_"$rname" default "$reward"
            sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AL009126.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AL009126.3_t4_maxsep_df.pkl positive 1 0 iYO844_block_"$rname" default "$reward"
        elif [[ "$typrd" == "flux" ]]; then
            sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1_t1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl negative 0 1 iLJ478_flux_$rname default $reward
            sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_002947.4_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_t2_maxsep_df.pkl negative 0 1 iJN1463_flux_$rname default $reward
            sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AM180355.1_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_t2_maxsep_df.pkl positive 0 1 iCN900_flux_$rname default $reward
            sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000913.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000913.3_t4_maxsep_df.pkl negative 0 1 iAF1260_flux_$rname default $reward
            sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/CP000148.1_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1_t4_maxsep_df.pkl negative 0 1 iAF987_flux_$rname default $reward
            sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000962.3_t5.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3_t5_maxsep_df.pkl positive 0 1 iNJ661_flux_$rname default $reward
            sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AL009126.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AL009126.3_t4_maxsep_df.pkl positive 0 1 iYO844_flux_$rname default $reward
        elif [[ "$typrd" == "fluxblock" ]]; then
            sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1_t1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl negative 1 1 iLJ478_fluxblock_$rname  default $reward
            sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_002947.4_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_t2_maxsep_df.pkl negative 1 1 iJN1463_fluxblock_$rname  default $reward
            sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AM180355.1_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_t2_maxsep_df.pkl positive 1 1 iCN900_fluxblock_$rname  default $reward
            sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000913.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000913.3_t4_maxsep_df.pkl negative 1 1 iAF1260_fluxblock_$rname  default $reward
            sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/CP000148.1_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1_t4_maxsep_df.pkl negative 1 1 iAF987_fluxblock_$rname  default $reward
            sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000962.3_t5.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3_t5_maxsep_df.pkl positive 1 1 iNJ661_fluxblock_$rname  default $reward
            sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AL009126.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AL009126.3_t4_maxsep_df.pkl positive 1 1 iYO844_fluxblock_$rname  default $reward

        fi
    done
done

# python mainclean_interation_v2.py --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1_t1.fasta --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl --file_type 7 --reward 0.5 --iter 10 --cpu 8 --gram negative --block_flage 0 --flux_flage 1 --name iLJ478_flux_p05 --threshold 8 --media default
# python mainclean_interation_v2.py --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_002947.4_t2.fasta --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_t2_maxsep_df.pkl --file_type 7 --reward 0.5 --iter 10 --cpu 8 --gram negative --block_flage 0 --flux_flage 1 --name iJN1463_flux_p05 --threshold 8 --media default
# python mainclean_interation_v2.py --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1_t1.fasta --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl --file_type 7 --reward 0.5 --iter 10 --cpu 8 --gram negative --block_flage 1 --flux_flage 1 --name iLJ478_fluxblock_p05 --threshold 8 --media default
# python mainclean_interation_v2.py --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000913.3_t4.fasta --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000913.3_t4_maxsep_df.pkl --file_type 7 --reward 0.5 --iter 10 --cpu 8 --gram negative --block_flage 1 --flux_flage 1 --name iAF1260_fluxblock_p05 --threshold 8 --media default
# python mainclean_interation_v2.py --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000962.3_t5.fasta --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3_t5_maxsep_df.pkl --file_type 7 --reward 0.5 --iter 10 --cpu 8 --gram positive --block_flage 1 --flux_flage 1 --name iNJ661_fluxblock_p05 --threshold 8 --media default
# python mainclean_interation_v2.py --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/AL009126.3_t4.fasta --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/AL009126.3_t4_maxsep_df.pkl --file_type 7 --reward 0.5 --iter 10 --cpu 8 --gram positive --block_flage 1 --flux_flage 1 --name iYO844_fluxblock_p05 --threshold 8 --media default
# python mainclean_interation_v2.py --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/AM180355.1_t2.fasta --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_t2_maxsep_df.pkl --file_type 7 --reward 0.7 --iter 10 --cpu 8 --gram positive --block_flage 0 --flux_flage 0 --name iCN900_allec_p07 --threshold 8 --media default
# python mainclean_interation_v2.py --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1_t1.fasta --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl --file_type 7 --reward 0.7 --iter 10 --cpu 8 --gram negative --block_flage 1 --flux_flage 0 --name iLJ478_block_p07 --threshold 8 --media default
# python mainclean_interation_v2.py --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/AM180355.1_t2.fasta --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_t2_maxsep_df.pkl --file_type 7 --reward 0.7 --iter 10 --cpu 8 --gram positive --block_flage 1 --flux_flage 0 --name iCN900_block_p07 --threshold 8 --media default
# python mainclean_interation_v2.py --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_002947.4_t2.fasta --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_t2_maxsep_df.pkl --file_type 7 --reward 0.7 --iter 10 --cpu 8 --gram negative --block_flage 1 --flux_flage 1 --name iJN1463_fluxblock_p07 --threshold 8 --media default
# python mainclean_interation_v2.py --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/AM180355.1_t2.fasta --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_t2_maxsep_df.pkl --file_type 7 --reward 0.7 --iter 10 --cpu 8 --gram positive --block_flage 1 --flux_flage 1 --name iCN900_fluxblock_p07 --threshold 8 --media default
# python mainclean_interation_v2.py --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1_t1.fasta --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl --file_type 7 --reward 1.0 --iter 10 --cpu 8 --gram negative --block_flage 0 --flux_flage 0 --name iLJ478_allec_1 --threshold 8 --media default
# python mainclean_interation_v2.py --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/AM180355.1_t2.fasta --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_t2_maxsep_df.pkl --file_type 7 --reward 1.0 --iter 10 --cpu 8 --gram positive --block_flage 0 --flux_flage 0 --name iCN900_allec_1 --threshold 8 --media default
# python mainclean_interation_v2.py --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000962.3_t5.fasta --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3_t5_maxsep_df.pkl --file_type 7 --reward 1.0 --iter 10 --cpu 8 --gram positive --block_flage 1 --flux_flage 0 --name iNJ661_block_1 --threshold 8 --media default
# python mainclean_interation_v2.py --input_file /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1_t1.fasta --cleanfile /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl --file_type 7 --reward 1.0 --iter 10 --cpu 8 --gram negative --block_flage 0 --flux_flage 1 --name iLJ478_flux_1 --threshold 8 --media default

# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1_t1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl negative 0 0 iLJ478_flux_p05 default 0.5
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_002947.4_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_t2_maxsep_df.pkl negative 0 0 iJN1463_flux_p05 default 0.5
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1_t1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl negative 1 0 iLJ478_fluxblock_p05 default 0.5
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000913.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000913.3_t4_maxsep_df.pkl negative 1 0 iAF1260_fluxblock_p05 default 0.5
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000962.3_t5.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3_t5_maxsep_df.pkl positive 1 0 iNJ661_fluxblock_p05 default 0.5
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AL009126.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AL009126.3_t4_maxsep_df.pkl positive 1 0 iYO844_fluxblock_p05 default 0.5
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AM180355.1_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_t2_maxsep_df.pkl positive 0 0 iCN900_allec_p07 default 0.7
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1_t1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl negative 1 0 iLJ478_block_p07 default 0.7
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AM180355.1_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_t2_maxsep_df.pkl positive 1 0 iCN900_block_p07 default 0.7
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_002947.4_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_t2_maxsep_df.pkl negative 1 1 iJN1463_fluxblock_p07 default 0.7
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AM180355.1_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_t2_maxsep_df.pkl positive 1 1 iCN900_fluxblock_p07 default 0.7
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1_t1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl negative 0 0 iLJ478_allec_1 default 1.0
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AM180355.1_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_t2_maxsep_df.pkl positive 0 0 iCN900_allec_1 default 1.0
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000962.3_t5.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3_t5_maxsep_df.pkl positive 1 0 iNJ661_block_1 default 1.0
# sbatch maincleaniterv2.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1_t1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl negative 0 1 iLJ478_flux_1 default 1.0

