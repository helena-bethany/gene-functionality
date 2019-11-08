# Use RIsearch2 to determine RNA:RNA interactions

#####################################################################

# $1 is the input ncRNA dataset file
# $2 is the test interaction database, if generated

IFS=$'\n'
d=$(date +%y%m%d) # date identifier
id=$( echo $1 | cut -d '.' -f 1 | cut -d '-' -f 4 )
id_count=$(( $id + 1 ))
calc() { awk "BEGIN{print $*}"; }

#####################################################################

# Generate interaction database, if needed

# Input is downloaded file from NCBI (make sure to download human data only)

# Promoters: https://www.ncbi.nlm.nih.gov/gene/?term=(Promoter+region)+AND+%22Homo+sapiens%22%5Bporgn%3A__txid9606%5D
# mRNA: https://www.ncbi.nlm.nih.gov/gene/?term=(mRNA)+AND+%22Homo+sapiens%22%5Bporgn%3A__txid9606%5D
# ncRNA: https://www.ncbi.nlm.nih.gov/gene/?term=(ncRNA)+AND+%22Homo+sapiens%22%5Bporgn%3A__txid9606%5D

get_ncbi_sequence() {
    
    top=$( grep "NC_" $1.txt | head -250 )
    rm -rf $1.fa

    for line in $top
    do
        chr=$( echo $line | cut -f 11 )
        start=$( echo $line | cut -f 13)
        end=$( echo $line | cut -f 14)
        sequence=$( grep -w "chromosome $chr" $2 | cut -f 2 | cut -c$start-$end )
        echo ">$chr-$start-$end" >> $1.fa
        echo $sequence >> $1.fa
    done
}

if [[ $2 == "interaction.suf" ]]
then
    
    grep -v 'Start' $1 | cut -d ',' -f 1,8 | tr ',' ' ' | perl -lane '{print ">$F[0]\n$F[1]"}' | tr -d '"' > data.fa
    interaction_database=$2

else

    read -p "Enter first NCBI text file: " file1
    read -p "Enter second NCBI text file: " file2
    read -p "Enter third NCBI text file: " file3
    
    file1_name=$( echo $file1 | cut -d '.' -f 1 )
    file2_name=$( echo $file2 | cut -d '.' -f 1 )
    file3_name=$( echo $file3 | cut -d '.' -f 1 )
    
    read -p "Enter CSV file for human genome: " genome
    
    get_ncbi_sequence $file1_name $genome
    get_ncbi_sequence $file2_name $genome
    get_ncbi_sequence $file3_name $genome
    
    grep -v 'Start' $1 | cut -d ',' -f 1,8 | tr ',' ' ' | perl -lane '{print ">$F[0]\n$F[1]"}' | tr -d '"' > data.fa
    
    grep -v ">" data.fa | shuf -n 250 > 250-seqs.fa

    for line in $( cat 250-seqs.fa )
    do
	    grep -B 1 $line data.fa >> ncrna-250.fa
    done
    
    cat $file1_name.fa $file2_name.fa $file3_name.fa ncrna-250.fa > interactionstest.fa 
    ./risearch2.x -c interactionstest.fa -o interaction.suf
    
    interaction_database=interaction.suf
    
fi

#####################################################################

# Use RIsearch2 to identify RNA:RNA interactions

risearch=$( locate risearch2.x | head -1 )

./risearch2.x -q data.fa -i $interaction_database
gunzip *.gz

FILES=./risearch_*.out
data=$(ls $FILES | sort -V)

echo RIsearchMIN, RIsearchAVE > risearch_intermediate.txt

for file in $data
do
    # Sums up every number in the 8th column
    count=$(awk '{s+=$8} END {print s}' $file)
    # Sorts the 8th column and takes the most negative number
    minimum=$(cut -f 8 $file | sort | head -1)
    # The number of lines in the file
    total=$(cat $file | wc -l)
    average=$(calc $count/$total)
    echo $minimum,$average >> risearch_intermediate.txt
done

paste -d ',' $1 risearch_intermediate.txt > $d-ncrna-dataset-$id_count.csv

#####################################################################

rm -rf data.fa
rm -rf 250-seqs.fa
rm -rf ncrna-250.fa

mkdir $d-risearch
mv risearch_RNA*.out $d-risearch/ 

