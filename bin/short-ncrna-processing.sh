#!/bin/bash
#
# Script Name: short-ncrna-processing.sh
#
# Author: Helena Cooper
# Last edited: 20/09/2020
#
# Description: This script filters the functional short ncRNAs and generates negative control sequences. 
#
# Input: $1 is the folder containing the additional files and dependencies not exported to $PATH (do not include '/' at the end)
#        $2 is the downloaded fasta file from RNAcentral containing the short ncRNAs. 
#        $3 is the downloaded chromosome coordinates file from RNAcentral for all ncRNAs in the database. Eg: Homo_sapiens.GRCh38.bed
#        $4 is the converted CSV file of the human genome. Eg: GRCh38.p13_genome.csv
#        $5 is the uniprot bed file from UCSC for filtering against known protein-coding genes. Eg: uni-prot-ucsc.bed
#        $6 is the ncrna bed file from GENCODE for filtering against known ncRNAs. Eg: gencode-ncrna-annotation.bed
#

###########################################################################################################################

######## General setup

# set -xv
IFS=$'\n'
d=$(date +%y%m%d) 

######## Delete files if they currently exist, as script continuously appends to the same file
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
rm -rf left-intersect
rm -rf right-intersect
rm -rf left-intersect.bed
rm -rf right-intersect.bed

###########################################################################################################################

# Reformat functional short ncRNAs

###########################################################################################################################

######## Convert input from FASTA to CSV file for easier parsability
[ -f $1/fasta_formatter ] && $1/fasta_formatter -i $1/$2 -o $d-converted.csv -t || fasta_formatter -i $1/$2 -o $d-converted.csv -t

######## Obtain chromosome coordinates for each functional ncRNA
count=0

for line in $(cat $d-converted.csv)
do
    ID=$( echo $line | cut -f 1 | cut -c1-18 )   # ncRNA ID
    sequence=$( echo $line | cut -f 2 )          # ncRNA Sequence

    if grep -q $ID $1/$3   # If ncRNA ID has chromosome coordinates available
    then
        meta=$( grep -m 1 $ID $1/$3 )
        chr=$( echo $meta | cut -f 1 )
        start=$( echo $meta | cut -f 2 )
        end=$( echo $meta | cut -f 3 )
        sequence_length=$(( $end - $start + 1 ))  # Add 1 to match up the sequence length to its true value

        if [ $chr == 'chrM' ] || [ $chr == 'chrY' ] # Removes ncRNA from mitochondria and Y chromosome
        then
            :
        elif [ "$sequence_length" -gt "3000" ]  # Removes sequences greater than 3000 bp
        then
            :
        else
            count=$(( $count + 1 ))
            echo RNA$count,Yes,$chr,$start,$end,$sequence >> ncrna-dataset  # Append data to temporary dataset
        fi

    else
        :

    fi

done

######## Generate functional dataset
echo ID,Functional,Chromosome,Start,End,Sequence > $d-functional-ncrna-dataset.csv

if [ "$count" -gt "1000" ]  # Filter for processing precursor miRNAs
then
    max=$( cat ncrna-dataset | wc -l )
    shuf -i 1-$max -n 89 > numbers   # Change this value depending on how many precursor miRNAs should be processed
    id_count=1
    for line in $( cat numbers )
    do
        info=$( sed -n "$line"p ncrna-dataset | cut -d ',' -f 2- )
        echo RNA$id_count,$info >> $d-functional-ncrna-dataset.csv
        id_count=$(( id_count + 1 ))
    done
else
    cat ncrna-dataset >> $d-functional-ncrna-dataset.csv
fi

######## Generate functional FASTA file
grep -v "Start" $d-functional-ncrna-dataset.csv | cut -d ',' -f 1,6 | tr ',' ' ' | perl -lane '{print ">$F[0]\n$F[1]"}' > $d-functional-ncrna-seq.fa

######## Report final number of functional sequences
other_count=$( grep -v "Start" $d-functional-ncrna-dataset.csv | wc -l )

echo
echo There were $other_count functional ncRNA found.
echo

###########################################################################################################################

# Create negative control ncRNA sequences from human genome

###########################################################################################################################

######## Generate sequences upstream and downstream of functional gene
text=$( grep -v "Start" $d-functional-ncrna-dataset.csv )

for line in $text
do
    ######## Extract chromosome coordinates and sequence length for functional ncRNA
    chr=$( echo $line | cut -d ',' -f 3 )
    chromo=$( echo $chr | tr -d "chr" )
    start=$( echo $line | cut -d ',' -f 4)
    end=$( echo $line | cut -d ',' -f 5 )
    length=$(( $end - $start ))

    ######## Generate null sequence 20,000 upstream of ncRNA
    left_end=$(( $start - 20000 ))
    left_start=$(( $left_end - $length ))
    left_sequence=$( grep -w "chromosome $chromo" $1/$4 | cut -f 2 | cut -c$left_start-$left_end )
    if [ -z $left_sequence ]  # If no sequence extracted, remove
    then
        :
    elif [[ $left_sequence == *"N"* ]]  # If null region contains undefined nuclotide(s), remove  
    then
        :
    else
        echo $chr,$left_start,$left_end,$left_sequence >> $d-left.csv
    fi

    ######## Generate null sequence 20,000 downstream of ncRNA
    right_start=$(( $end + 20000 ))
    right_end=$(( $right_start + $length ))
    right_sequence=$( grep -w "chromosome $chromo" $1/$4 | cut -f 2 | cut -c$right_start-$right_end )
    if [ -z $right_sequence ] # If no sequence extracted, remove
    then
        :
    elif [[ $right_sequence == *"N"* ]] # If null region contains undefined nuclotide(s), remove  
    then
        :
    else
        echo $chr,$right_start,$right_end,$right_sequence >> $d-right.csv
    fi
done

echo Negative control ncRNA sequences generated
echo

###########################################################################################################################

# Filter upstream negative control sequences using UniProt and GENCODE

###########################################################################################################################

######## Reformat data for bedtools
grep -v "Start" $d-left.csv | cut -d ',' -f 1,2,3 | tr ',' ' ' | perl -lane '{print "$F[0] $F[1] $F[2]"}' | tr ' ' '\t' > $d-left-coordinates.bed

######## Uniprot filtering
[ -f $1/bedtools ] && $1/bedtools intersect -a $1/$5 -b $d-left-coordinates.bed > $d-left-intersect.bed || bedtools intersect -a $1/$5 -b $d-left-coordinates.bed > $d-left-intersect.bed

######## Extract chromosome coordinates of overlapping negative control sequences
cat $d-left-intersect.bed | cut -f 1,2,3 | tr -d "chr" | tr '\t' ',' > $d-left-intersect
left_count=$( cat $d-left-intersect | wc -l )

######## GENCODE filtering
[ -f $1/bedtools ] && $1/bedtools intersect -a $1/$6 -b $d-left-coordinates.bed > left-intersect.bed || bedtools intersect -a $1/$6 -b $d-left-coordinates.bed > left-intersect.bed

######## Extract chromosome coordinates of overlapping negative control sequences
cat left-intersect.bed | cut -f 1,2,3 | tr -d "chr" | tr '\t' ',' > left-intersect

######## Remove negative control sequences that overlap with known functional sequences
echo ID,Functional-type,Chromosome,Start,End,Sequence > $d-negative-control-dataset.csv

for line in $( cat $d-left.csv )
do
    coordinates=$( echo $line | cut -d ',' -f 1,2,3 | tr -d "chr" )

    if grep -q $coordinates $d-left-intersect  # If sequence overlaps a known protein-coding gene, remove
    then
        :
    elif grep -q $coordinates left-intersect   # If sequenec overlaps a known ncRNA gene, remove
    then
        :
    else
        echo RNA$count,No,$line >> $d-negative-control-dataset.csv
        count=$(($count+1))  
    fi

done

###########################################################################################################################

# Filter downstream negative control sequences using UniProt and GENCODE

###########################################################################################################################

######## Reformat data for bedtools
grep -v "Start" $d-right.csv | cut -d ',' -f 1,2,3 | tr ',' ' ' | perl -lane '{print "$F[0] $F[1] $F[2]"}' | tr ' ' '\t' > $d-right-coordinates.bed

######## Uniprot filtering
[ -f $1/bedtools ] && $1/bedtools intersect -a $1/$5 -b $d-right-coordinates.bed > $d-right-intersect.bed || bedtools intersect -a $1/$5 -b $d-right-coordinates.bed > $d-right-intersect.bed

######## Extract chromosome coordinates of overlapping negative control sequences
cat $d-right-intersect.bed | cut -f 1,2,3 | tr -d "chr" | tr '\t' ',' > $d-right-intersect
right_count=$( cat $d-right-intersect | wc -l )

######## GENCODE filtering
[ -f $1/bedtools ] && $1/bedtools intersect -a $1/$6 -b $d-right-coordinates.bed > right-intersect.bed || bedtools intersect -a $1/$6 -b $d-right-coordinates.bed > right-intersect.bed

######## Extract chromosome coordinates of overlapping negative control sequences
cat right-intersect.bed | cut -f 1,2,3 | tr -d "chr" | tr '\t' ',' > right-intersect

######## Remove negative control sequences that overlap with known functional sequences

for line in $( cat $d-right.csv )
do
    coordinates=$( echo $line | cut -d ',' -f 1,2,3 | tr -d "chr" )

    if grep -q $coordinates $d-right-intersect  # If sequence overlaps a known protein-coding gene, remove
    then
        :
    elif grep -q $coordinates right-intersect   # If sequenec overlaps a known ncRNA gene, remove
    then
        :
    else
        echo NC$count,Negative-control,$line >> $d-negative-control-dataset.csv
        count=$(($count+1))  # Same count as with the upstream negative control sequences
    fi

done

######## Generate negative control FASTA file
grep -v "Start" $d-negative-control-dataset.csv | cut -d ',' -f 1,6 | tr ',' ' ' | perl -lane '{print ">$F[0]\n$F[1]"}' > $d-negative-control-seq.fa

######## Report number of generated negative control sequences, post filtering
total=$(( $left_count + $right_count ))
echo There were $total negative control sequences that overlapped with annotated coding regions.
echo
file_count=$( cat $d-negative-control-dataset.csv | wc -l )
count=$(( $file_count - 1 ))
echo There were $count non-overlapping negative control sequences generated.
echo

###########################################################################################################################

######## Delete excess files
rm -rf left-data
rm -rf right-data
rm -rf ncrna-data
rm -rf numbers
rm -rf $d-converted.csv
rm -rf $d-left.csv
rm -rf $d-right.csv
rm -rf $d-left-coordinates.bed
rm -rf $d-left-intersect.bed
rm -rf $d-left-intersect
rm -rf $d-right-coordinates.bed
rm -rf $d-right-intersect.bed
rm -rf $d-right-intersect
rm -rf left-intersect
rm -rf right-intersect
rm -rf left-intersect.bed
rm -rf right-intersect.bed

