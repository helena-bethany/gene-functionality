# Calculate number of times ncRNA appears within the human genome

#####################################################################

# $1 is the input ncRNA dataset file 

IFS=$'\n'
d=$(date +%y%m%d) # date identifier
id=$( echo $1 | cut -d '.' -f 1 | cut -d '-' -f 4 )
count=$(( $id + 1 ))

#####################################################################

# Create local human genome blast database

read -p "Make local blastn database? (y/n)" response

if [[ $response == "y" ]]
then
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_genomic.fna.gz
    gunzip GCF_000001405.38_GRCh38.p12_genomic.fna.gz
    makeblastdb -in GCF_000001405.38_GRCh38.p12_genomic.fna -dbtype nucl -parse_seqids -out human_genome
else
    :
fi

#####################################################################

# Run blastn to identify number of matches within human genome

echo RepeatsTotal,RepeatsChr,Repeats100Query > $d-repeats.csv
grep -v 'Start' $1 | cut -d ',' -f 1,8 > sequences

for line in $(cat sequences)
do

    echo $line | tr ',' ' ' | perl -lane '{print ">$F[0]\n$F[1]"}' > rna.fa
    blastn -query rna.fa -db human_genome -outfmt "10 qaccver saccver pident" > output.csv
    total=$(cat output.csv | wc -l)
    chromo=$(cut -d ',' -f 2 output.csv | grep 'NC' | uniq | wc -l)
    hundred=$(cut -d ',' -f 3 output.csv | grep '100.000' | wc -l)
    echo $total, $chromo, $hundred >> $d-repeats.csv
done

paste -d ',' $1 $d-repeats.csv > $d-ncrna-dataset-$count.csv

#####################################################################

rm -rf rna.fa
rm -rf output.csv
rm -rf sequences

