
# Test reward scale from 0.1 to 1.0
echo "Starting reward tests"
# rewards=(0.05 0.1 0.2 0.4 0.6 0.8 1.0 2.0)
# rewards=(0.9)
# rewards=(0.05 1.5 2.0 3.0 )
# rewards=(0.05 0.2 5.0)
rewards=( 0.1 0.4 0.8 1.0 3.0)
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
            sbatch maincleaniterv3.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1_t1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl negative 0 0 iLJ478_allec_"$rname" default "$reward"
            sbatch maincleaniterv3.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_002947.4_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_t2_maxsep_df.pkl negative 0 0 iJN1463_allec_"$rname" default "$reward"
            sbatch maincleaniterv3.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AM180355.1_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_t2_maxsep_df.pkl positive 0 0 iCN900_allec_"$rname" default "$reward"
            sbatch maincleaniterv3.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000913.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000913.3_t4_maxsep_df.pkl negative 0 0 iAF1260_allec_"$rname" default "$reward"
            sbatch maincleaniterv3.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/CP000148.1_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1_t4_maxsep_df.pkl negative 0 0 iAF987_allec_"$rname" default "$reward"
            sbatch maincleaniterv3.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000962.3_t5.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3_t5_maxsep_df.pkl positive 0 0 iNJ661_allec_"$rname" default "$reward"
            sbatch maincleaniterv3.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AL009126.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AL009126.3_t4_maxsep_df.pkl positive 0 0 iYO844_allec_"$rname" default "$reward"
        elif [[ "$typrd" == "block" ]]; then
            sbatch maincleaniterv3.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1_t1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl negative 1 0 iLJ478_block_"$rname" default "$reward"
            sbatch maincleaniterv3.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_002947.4_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_t2_maxsep_df.pkl negative 1 0 iJN1463_block_"$rname" default "$reward"
            sbatch maincleaniterv3.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AM180355.1_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_t2_maxsep_df.pkl positive 1 0 iCN900_block_"$rname" default "$reward"
            sbatch maincleaniterv3.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000913.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000913.3_t4_maxsep_df.pkl negative 1 0 iAF1260_block_"$rname" default "$reward"
            sbatch maincleaniterv3.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/CP000148.1_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1_t4_maxsep_df.pkl negative 1 0 iAF987_block_"$rname" default "$reward"
            sbatch maincleaniterv3.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000962.3_t5.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3_t5_maxsep_df.pkl positive 1 0 iNJ661_block_"$rname" default "$reward"
            sbatch maincleaniterv3.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AL009126.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AL009126.3_t4_maxsep_df.pkl positive 1 0 iYO844_block_"$rname" default "$reward"
        elif [[ "$typrd" == "flux" ]]; then
            sbatch maincleaniterv3.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1_t1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl negative 0 1 iLJ478_flux_$rname default $reward
            sbatch maincleaniterv3.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_002947.4_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_t2_maxsep_df.pkl negative 0 1 iJN1463_flux_$rname default $reward
            sbatch maincleaniterv3.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AM180355.1_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_t2_maxsep_df.pkl positive 0 1 iCN900_flux_$rname default $reward
            sbatch maincleaniterv3.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000913.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000913.3_t4_maxsep_df.pkl negative 0 1 iAF1260_flux_$rname default $reward
            sbatch maincleaniterv3.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/CP000148.1_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1_t4_maxsep_df.pkl negative 0 1 iAF987_flux_$rname default $reward
            sbatch maincleaniterv3.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000962.3_t5.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3_t5_maxsep_df.pkl positive 0 1 iNJ661_flux_$rname default $reward
            sbatch maincleaniterv3.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AL009126.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AL009126.3_t4_maxsep_df.pkl positive 0 1 iYO844_flux_$rname default $reward
        elif [[ "$typrd" == "fluxblock" ]]; then
            sbatch maincleaniterv3.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000853.1_t1.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000853.1_t1_maxsep_df.pkl negative 1 1 iLJ478_fluxblock_$rname  default $reward
            sbatch maincleaniterv3.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_002947.4_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_002947.4_t2_maxsep_df.pkl negative 1 1 iJN1463_fluxblock_$rname  default $reward
            sbatch maincleaniterv3.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AM180355.1_t2.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AM180355.1_t2_maxsep_df.pkl positive 1 1 iCN900_fluxblock_$rname  default $reward
            sbatch maincleaniterv3.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000913.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000913.3_t4_maxsep_df.pkl negative 1 1 iAF1260_fluxblock_$rname  default $reward
            sbatch maincleaniterv3.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/CP000148.1_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/CP000148.1_t4_maxsep_df.pkl negative 1 1 iAF987_fluxblock_$rname  default $reward
            sbatch maincleaniterv3.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/NC_000962.3_t5.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/NC_000962.3_t5_maxsep_df.pkl positive 1 1 iNJ661_fluxblock_$rname  default $reward
            sbatch maincleaniterv3.sh /ibex/user/niuk0a/CLEAN/app/data/inputs/AL009126.3_t4.fasta /ibex/user/niuk0a/CLEAN/app/results/inputs/AL009126.3_t4_maxsep_df.pkl positive 1 1 iYO844_fluxblock_$rname  default $reward

        fi
    done
done

# python mainclea