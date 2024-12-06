

# download from uniprot
datalink=https://ftp.uniprot.org/pub/databases/uniprot/previous_major_releases/release-20{year}_{month}/relnotes.txt
for year in {10..24}; do
    if year > 15; then
        month=04
    else    
        month=01
    fi
    datalink=https://ftp.uniprot.org/pub/databases/uniprot/previous_major_releases/release-20${year}_${month}/relnotes.txt
    echo $datalink
#   for month in {01..12}; do
    wget -O - $datalink -O uniprot_${year}_${month}.txt
done