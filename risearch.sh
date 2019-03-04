#Run RIsearch2 but in an automated way

#./risearch2.x -c Downloads/interactionstest.fa -o interaction.suf
grep -v 'Start' $1 | cut -d ',' -f 1,8 | tr ',' ' ' | perl -lane '{print ">$F[0]\n$F[1]"}' > data.fa
./risearch2.x -q data.fa -i 181217interaction.suf
gunzip *.gz

FILES=/home/helenacooper/risearch_*.out
data=$(ls $FILES | sort -V)
IFS=$'\n'
echo RIsearchMIN, RIsearchAVE > 181217risearch.txt

for file in $data
do
    python3 average.py $file RIsearch >> 181217risearch.txt
done

paste -d ',' $1 181217risearch.txt > $2
rm -rf risearch_*.out
rm -rf data.fa

#./RIblast db -i Downloads/interactionstest.fa -o test_db
#./RIblast ris -i DATA/181211snrna.fa -o output.txt -d test_db
# python3 average.py output.txt RIblast >> 181214riblast.txt
