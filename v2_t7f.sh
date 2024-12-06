files=$(ls /ibex/user/niuk0a/funcarve/cobra/ | grep "v2_t7rewards.*\.err$")
echo "Found ${#files[@]} files"
for file in $files; do
    # echo $file
# if file ends with .err and start with 'v2_t7rewards", then check the file if is emp
    if [[ -s $file ]]; then
        # echo "$file has data"
      ## get the outfile of the file
        outfile=$(echo $file | cut -d. -f2)
        outfile=v2_t7rewards.${outfile}.out
        cat $outfile
    else
        rm $file
        
    fi
    
   
done