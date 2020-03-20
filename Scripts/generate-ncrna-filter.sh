# Creation of the training/test dataset

set -xv

#####################################################################

# $1 is the downloaded fasta file from RNAcentral
# $2 is the downloaded chromosome coordinates file
# $3 is the csv file of the human genome
# $4 is the uniprot bed file from UCSC

IFS=$'\n'
d=$(date +%y%m%d)  # date identifier

#####################################################################

rm -rf $d-converted.csv
rm -rf ncrna-data
rm -rf $d-left.csv
rm -rf $d-right.csv
rm -rf $d-left-coordinates.bed
rm -rf $d-left-intersect.bed
rm -rf $d-left-intersect
rm -rf $d-right-coordinates.bed
rm -rf $d-right-intersect.bed
rm -rf $d-right-intersect
rm -rf numbers

#####################################################################

# Convert FASTA file to CSV

fasta_formatter -i $1 -o $d-converted.csv -t

#####################################################################

# Obtain chromosome coordinates for each functional ncRNA

count=0
max=$( cat $d-converted.csv | wc -l )
shuf -i 1-$max -n 1000 > numbers

echo ID,Function-type,Chromosome,Start,End,Sequence > ncrna-data

if [ "$max" -gt "1000" ]
then

    for line in $( cat numbers )
    do

        info=$( sed -n "$line"p $d-converted.csv )
        ID=$( echo $info | cut -f 1 | cut -c1-18 )
        sequence=$( echo $info | cut -f 2 )

        if grep -q $ID $2
        then
            chr=$( grep -m 1 $ID $2 | cut -f 1 )
            start=$( grep -m 1 $ID $2 | cut -f 2 )
            end=$( grep -m 1 $ID $2 | cut -f 3 )
            count=$(( $count + 1 ))

            echo $ID,Functional-ncRNA,$chr,$start,$end,$sequence >> ncrna-data

        else
            :

        fi
    done

else

    for line in $(cat $d-converted.csv)
    do

        ID=$( echo $line | cut -f 1 | cut -c1-18 )
        sequence=$( echo $line | cut -f 2 )

        if grep -q $ID $2
        then
            chr=$( grep -m 1 $ID $2 | cut -f 1 )
            start=$( grep -m 1 $ID $2 | cut -f 2 )
            end=$( grep -m 1 $ID $2 | cut -f 3 )
            count=$(( $count + 1 ))

            echo $ID,Functional-ncRNA,$chr,$start,$end,$sequence >> ncrna-data

        else

        fi

    done
fi

grep -v "chrM\|chrY" ncrna-data > $d-functional-ncrna-dataset.csv

grep -v "Start" $d-functional-ncrna-dataset.csv | cut -d ',' -f 1,6 | tr ',' ' ' | perl -lane '{print ">$F[0]\n$F[1]"}' > $d-functional-ncrna-seq.fa

echo
echo There were $count functional ncRNA found.
echo

#####################################################################

# Create negative control ncRNA from human genome

text=$( grep -v "Start" $d-functional-ncrna-dataset.csv )

for line in $text
do
    #Extract chromosome coordinates and sequence length for functional ncRNA
    chr=$( echo $line | cut -d ',' -f 3 )
    chromo=$( echo $chr | tr -d "chr" )
    start=$( echo $line | cut -d ',' -f 4)
    end=$( echo $line | cut -d ',' -f 5 )
    length=$(( $end - $start ))

    #Generate null sequence 20,000 upstream of ncRNA
    left_end=$(( $start - 20000 ))
    left_start=$(( $left_end - $length ))
    left_sequence=$( grep -w "chromosome $chromo" $3 | cut -f 2 | cut -c$left_start-$left_end )
    echo $chr,$left_start,$left_end,$left_sequence >> left-data

    #Generate null sequence 20,000 downstream of ncRNA
    right_start=$(( $end + 20000 ))
    right_end=$(( $right_start + $length ))
    right_sequence=$( grep -w "chromosome $chromo" $3 | cut -f 2 | cut -c$right_start-$right_end )
    echo $chr,$right_start,$right_end,$right_sequence >> right-data
done

grep -v "chrM\|chrY" left-data > $d-left.csv
grep -v "chrM\|chrY" right-data > $d-right.csv

echo Negative control ncRNA sequences generated
echo

#####################################################################

# Filter left negative control sequences using UniProt

grep -v "Start" $d-left.csv | cut -d ',' -f 1,2,3 | tr ',' ' ' | perl -lane '{print "$F[0] $F[1] $F[2]"}' | tr ' ' '\t' > $d-left-coordinates.bed
bedtools intersect -a $4 -b $d-left-coordinates.bed > $d-left-intersect.bed
cat $d-left-intersect.bed | cut -f 1,2,3 | tr -d "chr" | tr '\t' ',' > $d-left-intersect
left_count=$( cat $d-left-intersect | wc -l )

echo ID,Functional-type,Chromosome,Start,End,Sequence > $d-negative-control-dataset.csv

count=1

for line in $( cat $d-left.csv )
do

    coordinates=$( echo $line | cut -d ',' -f 1,2,3 | tr -d "chr" )

    if grep -q $coordinates $d-left-intersect
    then
        :
    else
        echo NC$count,Negative-control,$line >> $d-negative-control-dataset.csv
        count=$(($count+1))
    fi

done

#####################################################################

# Filter right negative control sequences using UniProt

grep -v "Start" $d-right.csv | cut -d ',' -f 1,2,3 | tr ',' ' ' | perl -lane '{print "$F[0] $F[1] $F[2]"}' | tr ' ' '\t' > $d-right-coordinates.bed
bedtools intersect -a $4 -b $d-right-coordinates.bed > $d-right-intersect.bed
cat $d-right-intersect.bed | cut -f 1,2,3 | tr -d "chr" | tr '\t' ',' > $d-right-intersect
right_count=$( cat $d-right-intersect | wc -l )

for line in $( cat $d-right.csv )
do

    coordinates=$( echo $line | cut -d ',' -f 1,2,3 | tr -d "chr" )

    if grep -q $coordinates $d-right-intersect
    then
        :
    else
        echo NC$count,Negative-control,$line >> $d-negative-control-dataset.csv
        count=$(($count+1))
    fi

done

grep -v "Start" $d-negative-control-dataset.csv | cut -d ',' -f 1,6 | tr ',' ' ' | perl -lane '{print ">$F[0]\n$F[1]"}' > $d-negative-control-seq.fa

total=$(( $left_count + $right_count ))
echo There were $total negative control sequences that overlapped with annotated coding regions.
echo
count=$(( $count - 1 ))
echo There were $count non-overlapping negative control sequences generated.
echo
