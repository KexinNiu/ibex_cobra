python test.py --input_file ../data/input/CP000148.1_t4_maxsep_df.pkl --file_type 1 --reward 0.2 \
    --iter 4 --block_flage 0 --flux_flage 0 --media default --tasks 0 --min_frac 0.01 --max_frac 0.5 \
    --threshold 5 --upper 15 --lower 5 --maxweight 100 --minweight 0.0 \
    --gram negative --name test_cp148 --cpu 1 --gapfill yes --exchange 1 --test no 

python funcarve_main_gmm.py --input_file ../data/input/CP000148.1_t4_maxsep_df.pkl --file_type 1 --reward 0.2 \
    --iter 4 --block_flage 0 --flux_flage 0 --media default --tasks 0 --min_frac 0.01 --max_frac 0.5 \
    --threshold 5 --upper 15 --lower 5 --maxweight 100 --minweight 0.0 \
    --gram negative --name test_cp148 --cpu 1 --gapfill yes --exchange 1 --test no