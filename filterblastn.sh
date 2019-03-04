#A more accurate way of getting blast to work

echo RepeatsTotal,RepeatsChr,Repeats100Query > repeats.csv
grep -v 'Start' $1 | cut -d ',' -f 1,8 > sequences
#grep -v 'Start' $1 | cut -d ',' -f 1,69 > sequences

IFS=$'\n'
for line in $(cat sequences)
do

    echo $line | tr ',' ' ' | perl -lane '{print ">$F[0]\n$F[1]"}' > rna.fa
    blastn -query rna.fa -db human_genome -outfmt "10 qaccver saccver pident" > output.csv
    total=$(cat output.csv | wc -l)
    chromo=$(cut -d ',' -f 2 output.csv | grep 'NC' | uniq | wc -l)
    hundred=$(cut -d ',' -f 3 output.csv | grep '100.000' | wc -l)
    echo $total, $chromo, $hundred >> repeats.csv
done

paste -d ',' $1 repeats.csv > $2.csv
rm -rf rna.fa
rm -rf output.csv
rm -rf sequences

