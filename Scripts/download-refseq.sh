# Obtain sequence from NCBI using RefSeq IDs

#set -xv

rm -rf numbers
rm -rf genes.fa
rm -rf seq
rm -rf summary

# $1 is the text file with the RefSeq IDs
# $2 is the human genome GFF file

####################################################

# Script set up

IFS=$'\n'
d=$(date +%y%m%d)
calc() { awk "BEGIN{print $*}"; }

rm -rf $d-protein-coding-genes.fa

echo ID,Name,Functional,Chromosome,Start,End > summary

####################################################

# Select 1000 random RefSeq IDs

max=$( cat $1 | wc -l )

shuf -i 1-$max -n 1000 > numbers

####################################################

# Obtain gene sequences and locations

count=1

for line in $( cat numbers )
do

    id=$( sed -n "$line"p $1 )

    coords=$( grep -E $id $2 | head -1 | cut -f 1,4,5 )

    chr=$( echo $coords | cut -d ' ' -f 1 | tr -d "NC_" | cut -d '.' -f 1 | cut -c5,6 )
    test=$( echo $chr | cut -c1 )
    other=$( echo $chr | cut -c2 )

    if [ "$test" -eq "0" ]    # Remove zero before each chromosome number (ie: 01 becomes 1)
    then
        chr=$other
    else
        :
    fi

    start=$( echo $coords | cut -f 2 )
    end=$( echo $coords | cut -f 3 )

    echo $id,RNA$count,Yes,$chr,$start,$end >> summary

    esearch -db nucleotide -query "$id" | efetch -format fasta >> genes.fa
    count=$(( $count + 1 ))
done

#####################################################

# Reformat so sequence on one line

n=$( grep -c ">" genes.fa )
num=$(( n*2 ))

cat genes.fa | awk '/^>/ {printf("\n%s\n",$1);next; } { printf("%s",$1);}  END {printf("\n");}' | tail -"$num" > formatted-fasta.fa

#####################################################

# Add sequence to summary csv file

echo Sequence >> seq

for line in $( cat summary )
do

    name=$( echo $line | cut -d ',' -f 1 )
    sequence=$( grep -A 1 $name formatted-fasta.fa | tail -1 )
    # Only echo sequence if it's present
    [ ! -z "$sequence" ] && echo $sequence >> seq

done

summary_len=$( cat summary | wc -l )
seq_len=$( cat summary | wc -l)
cat summary | cut -d ',' -f 2- > interim

[[ $summary_len -eq $seq_len ]] && paste -d ',' interim seq > $d-protein-coding-dataset.csv

grep -v "Start" $d-protein-coding-dataset.csv | cut -d ',' -f 1,6 | tr ',' ' ' | perl -lane '{print ">$F[0]\n$F[1]"}' > $d-protein-coding-seq.fa

######################################################

rm -rf numbers
rm -rf genes.fa
rm -rf seq
rm -rf summary
rm -rf interim
rm -rf formatted-fasta.fa
