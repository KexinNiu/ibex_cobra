## grep last line of files in current directory

# for file in *; do
#     ## end with .out
#     if [[ $file == *.out ]]; then
#          tail -n 1 $file >> lastline.txt
#     fi
# done

for file in *; do
    ## end with .err AND start with test_meth1
    if [[ $file == test_meth3*.out ]]; then
    
        grep 'result:' $file >> testmeth.txt
        #  tail -n 1 $file >> lastline.txt
    fi
done

sbatch 01run.sh 7777 1
sbatch 01run.sh 7777 2
sbatch 01run.sh 7777 3
sbatch 01run.sh 7777 4
sbatch 01run.sh 7777 5
sbatch 01run.sh 7777 6
sbatch 01run.sh 7777 7
sbatch 01run.sh 7777 8
sbatch 01run.sh 7777 9
sbatch 01run.sh 6666 1
sbatch 01run.sh 6666 2
sbatch 01run.sh 6666 3
sbatch 01run.sh 6666 4
sbatch 01run.sh 6666 5
sbatch 01run.sh 6666 6
sbatch 01run.sh 6666 7
sbatch 01run.sh 6666 8
sbatch 01run.sh 6666 9
sbatch 01run.sh 5555 1
sbatch 01run.sh 5555 2
sbatch 01run.sh 5555 3
sbatch 01run.sh 5555 4
sbatch 01run.sh 5555 5
sbatch 01run.sh 5555 6
sbatch 01run.sh 5555 7
sbatch 01run.sh 5555 8
sbatch 01run.sh 5555 9
sbatch 01run.sh 4444 1
sbatch 01run.sh 4444 2
sbatch 01run.sh 4444 3
sbatch 01run.sh 4444 4
sbatch 01run.sh 4444 5
sbatch 01run.sh 4444 6
sbatch 01run.sh 4444 7
sbatch 01run.sh 4444 8
sbatch 01run.sh 4444 9

sbatch 01run.sh 7777 10
sbatch 01run.sh 6666 10
sbatch 01run.sh 5555 10
sbatch 01run.sh 4444 10

sbatch 01run.sh 7777 6
sbatch 01run.sh 6666 6
sbatch 01run.sh 5555 6
sbatch 01run.sh 4444 6