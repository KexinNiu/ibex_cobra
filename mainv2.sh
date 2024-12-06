#!/bin/bash
#SBATCH --job-name=v2t2df2
#SBATCH -o v2t2df2_%A_%a.out
#SBATCH -e v2t2df2_%A_%a.err
#SBATCH --mem=150G
#SBATCH --cpus-per-task=8
#SBATCH --time=10:00:00
#SBATCH --array=1-45

# Set the directory containing input files
# fastadr=/ibex/user/niuk0a/CLEAN/app/data/inputs/
fastadr=/ibex/user/niuk0a/CLEAN/app/data/inputs/testgenomes/t2/
cleandf=/ibex/user/niuk0a/CLEAN/app/results/inputs/t2/
# cleandf=/ibex/user/niuk0a/CLEAN/app/results/inputs/
gramf=/ibex/user/niuk0a/CLEAN/app/data/inputs/testgenomes/taxo_gram_t2.csv
# gramf=/ibex/user/niuk0a/CLEAN/app/data/inputs/testgenomes/taxo_gram_t2_addon.csv

# Create a list of files matching "UP*.fasta"
files=($(ls $fastadr | grep "UP.*\.fasta$"))
echo "Found ${#files[@]} files"
# Calculate the number of files
num_files=${#files[@]}
echo "Number of files: $num_files"§
# Get the file corresponding to the current SLURM_ARRAY_TASK_ID
# Subtract 1 from SLURM_ARRAY_TASK_ID to match zero-based indexing
index=$((SLURM_ARRAY_TASK_ID - 1))

# Ensure the index is within bounds
if [[ $index -ge 0 && $index -lt $num_files ]]; then
    file=${files[$index]}
    name="${file%.fasta}"  # Remove the .fasta extension for naming
    echo "Running $name (file: $file, index: $SLURM_ARRAY_TASK_ID)" 
    
    taxo=$(echo $name | cut -d_ -f2)
    gram=$(grep $taxo $gramf | cut -d, -f3)
    # if gram ==none then continue
    if [ "$gram" != "none" ]; then
        exit 1
    elif [ "$gram" == "none" ]; then
        echo "Running $name (file: $file, index: $SLURM_ARRAY_TASK_ID,taxo:$taxo gram: $gram)"
        python mainclean_interation_v2.py --input_file $fastadr$file --cleanfile $cleandf${name}_maxsep_df.pkl --block_flage 0 --flux_flage 0 --file_type 7 --reward 0.1 --iter 10 --cpu 8 --gram negative --name $name --media default 

    fi
    echo "Running $name (file: $file, index: $SLURM_ARRAY_TASK_ID,taxo:$taxo gram: $gram)"
    #Running UP000000238_349521 (file: UP000000238_349521.fasta, index: 1)
    echo "python mainclean_interation_v2.py --input_file $fastadr$file --cleanfile $cleandf${name}_maxsep_df.pkl --block_flage 0 --flux_flage 0 --file_type 7 --reward 0.1  --iter 10 --cpu 8 --gram $gram --name $name"
    
    # python mainclean_interation_v2.py --input_file $fastadr$file --cleanfile $cleandf${name}_maxsep_df.pkl --block_flage 0 --flux_flage 0 --file_type 7 --reward 0.1 --iter 10 --cpu 8 --gram $gram --name $name --media default 

    
    # Run the Python script with the specified arguments
    # echo "$dr/$file"
    # python CLEAN_infer_fasta_train.py --fasta_data "$name" --traindata split70_bactaxo_train2
    echo "Done with $name"
else
    echo "Index $index is out of bounds for the files array (num_files: $num_files)"
    exit 1
fi


# dr=/ibex/user/niuk0a/CLEAN/app/data/inputs
# files=$(ls $dr)
# SLURM_ARRAY_TASK_ID=0
# for file in $files; do
#     if [[ $file == *"UP"* ]]; then
#         basename=$(basename $file)
#         name=$(basename $basename .fasta)
#         echo running $name
#         $SLURM_ARRAY_TASK_ID=$((SLURM_ARRAY_TASK_ID+1))
#         python CLEAN_infer_fasta_train.py --fasta_data "$name" --traindata split70_bactaxo_train2
#         echo "done with $name"
#         continue
#     fi
   
# done