# for xml file in current folder run memote
# for file in *.xml; do
#     echo "Running $file"
#     # memote run $file --filename "report_$file.json"
#     sbatch ../memote_sb.sh $file "report_$file.json"
# done

modelnames="iLJ478 iJN1463 iCN900 iAF1260 iAF987 iNJ661 iYO844"
for modelname in $modelnames; do
    wget https://www.metanetx.org/cgi-bin/mnxget/user_model/bigg_$modelname.chemicals.tsv
    wget https://www.metanetx.org/cgi-bin/mnxget/user_model/bigg_$modelname.reactions.tsv

done
