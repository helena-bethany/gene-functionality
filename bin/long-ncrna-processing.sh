#!/bin/bash
#
# Script Name: long-ncrna-processing.sh
#
# Author: Helena Cooper
# Last edited: 30/07/2020
#
# Description: This script filters the functional long ncRNAs and generates negative control sequences. 
#
# Input: $1 is the folder containing the additional files and dependencies not exported to $PATH (do not include '/' at the end)
#        $2 is the downloaded fasta file from RNAcentral containing the lncRNAs. 
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
rm -rf ncrna-dataset-*
rm -rf $d-left.csv
rm -rf $d-right.csv
rm -rf $d-left-coordinates.bed
rm -rf $d-left-intersect.bed
rm -rf $d-left-intersect
rm -rf $d-right-coordinates.bed
rm -rf $d-right-intersect.bed
rm -rf $d-right-intersect
rm -rf true-coordinates
rm -rf $d-all-exon-coords
rm -rf left-intersect.bed
rm -rf left-intersect
rm -rf right-intersect.bed
rm -rf right-intersect

###########################################################################################################################

# Reformat functional short ncRNAs

###########################################################################################################################

######## Convert input from FASTA to CSV file for easier parsability
[ -f $1/fasta_formatter ] && $1/fasta_formatter -i $1/$2 -o $d-converted.csv -t || fasta_formatter -i $1/$2 -o $d-converted.csv -t

######## Obtain chromosome coordinates for each functional ncRNA
count=1
max=$( cat $d-converted.csv | wc -l )
shuf -i 1-$max -n $max > numbers   # Very easy for sequences to be filtered out, so have 1000 cap later.

echo ID,Functional,Chromosome,Start,End,Sequence > $d-functional-ncrna-exon2-dataset.csv
echo ID,Functional,Chromosome,Start,End,Sequence > $d-functional-ncrna-exon3-dataset.csv

for line in $( cat $d-converted.csv )
do

    ID=$( echo $line | cut -f 1 | cut -c1-18 )
    sequence=$( echo $line | cut -f 2 )
    # Parse RNAcentral file
    if grep -q $ID $1/$3
    then
        meta=$( grep -m 1 $ID $1/$3 )
        chr=$( echo $meta | cut -f 1 )
        start=$( echo $meta | cut -f 2 )
        true_end=$( echo $meta | cut -f 3 )
        len=$( echo $meta | cut -f 11 | cut -d ',' -f 1 )
        end=$(( $start + $len ))
        exon_count=$( echo $meta | cut -f 10 )

        if [ -z $sequence ]   # If no sequence available for lncRNA, remove
        then
            :
	elif [ $chr == 'chrM' ] || [ $chr == 'chrY' ]  # Removes ncRNA from mitochondria and Y chromosome
        then
            :
	elif [ "$count" -gt "1000" ]   # Only need 1000 functional sequences at most
        then
            :
        elif [ "$exon_count" -lt "3" ]   # Need at least three exons
        then
            :
	else
            len_one=$( echo $meta | cut -f 11 | cut -d ',' -f 1 )
            len_two=$( echo $meta | cut -f 11 | cut -d ',' -f 2 )
            len_three=$( echo $meta | cut -f 11 | cut -d ',' -f 3 )
            two=$( echo $meta | cut -f 12 | cut -d ',' -f 2 ) # Start position with introns, but RNAcentral seq without
            three=$( echo $meta | cut -f 12 | cut -d ',' -f 3 ) # These are relative to StartPos, not genome coords
            start_two=$(( $true_start + $two + 1 ))  # +1 to account for the first relative start pos being 0
            start_three=$(( $true_start + $three + 1 ))   # Start position relative to chromosome coordinates
            end_two=$(( $start_two + $len_two ))
            end_three=$(( $start_three + $len_three ))
             
            intronless_two=$(( $len_one + 1 ))  # add one to previous end position (total length of first exon) to calculate 
                                                # intronless start position of exon two
            total_two=$(( $intronless_two + $len_two ))
            seq_two=$( echo $sequence | cut -c$intronless_two-$total_two )

            intronless_three=$(( $total_two + 1 ))
            total_three=$(( $intronless_three + $len_three ))
            seq_three=$( echo $sequence | cut -c$intronless_three-$total_three )

            if [ "$len_two" -lt "3000" ] || [ "$len_three" -lt "3000" ]   # Sequence length cannot be greater than 3000 bp
            then
                echo RNA$count,Yes,$chr,$start_two,$end_two,$seq_two >> $d-functional-lncrna-exon2-dataset.csv
                echo RNA$count,Yes,$chr,$start_three,$end_three,$seq_three >> $d-functional-lncrna-exon3-dataset.csv
                if [ "$true_start" -gt "$true_end" ]         # Take reverse strand into account
                then
                    echo $chr,$true_end,$true_start,$len_two,$len_three >> true-coordinates   # Coordinates of whole gene, rather than exons
                else
                    echo $chr,$true_start,$true_end,$len_two,$len_three >> true-coordinates
                fi
                count=$(( $count + 1 ))
            else
                :
            fi
        fi
    else
        :
    fi
done

######## Generate functional FASTA file
grep -v "Start" $d-functional-lncrna-exon2-dataset.csv | cut -d ',' -f 1,6 | tr ',' ' ' | perl -lane '{print ">$F[0]\n$F[1]"}' > $d-functional-ncrna-exon2-seq.fa
grep -v "Start" $d-functional-;ncrna-exon3-dataset.csv | cut -d ',' -f 1,6 | tr ',' ' ' | perl -lane '{print ">$F[0]\n$F[1]"}' > $d-functional-ncrna-exon3-seq.fa

######## Report final number of functional sequences
other_count=$( grep -v "Start" $d-functional-lncrna-exon2-dataset.csv | wc -l )

echo
echo There were $other_count functional ncRNA found.
echo

###########################################################################################################################

# Create negative control lncRNA sequences from human genome

###########################################################################################################################

######## Generate sequences upstream and downstream of functional gene

for line in $( cat true-coordinates )
do
    ######## Extract chromosome coordinates and sequence length for functional lncRNA
    chr=$( echo $line | cut -d ',' -f 1 )
    chromo=$( echo $chr | tr -d "chr" )
    start=$( echo $line | cut -d ',' -f 2)
    end=$( echo $line | cut -d ',' -f 3 )
    left_length=$( echo $line | cut -d ',' -f 4 )   # null sequence still generated up or downstream of original coordinates
    right_length=$( echo $line | cut -d ',' -f 5 )  # length is based on exons to match short ncRNA

    ######## Generate null sequence 20,000 upstream of lncRNA
    left_end=$(( $start - 20000 ))
    left_start=$(( $left_end - $left_length ))
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
    
    ######## Generate null sequence 20,000 downstream of lncRNA
    right_start=$(( $end + 20000 ))
    right_end=$(( $right_start + $right_length ))
    right_sequence=$( grep -w "chromosome $chromo" $1/$4 | cut -f 2 | cut -c$right_start-$right_end )
    if [ -z $right_sequence ]  # If no sequence extracted, remove
    then
        :
    elif [[ $right_sequence == *"N"* ]]  # If null region contains undefined nuclotide(s), remove 
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
echo ID,Functional,Chromosome,Start,End,Sequence > $d-negative-control-dataset.csv
count=$( tail -1 $d-functional-lncrna-exon2-dataset.csv | cut -d ',' -f 1 | tr -d "RNA" )

for line in $( cat $d-left.csv )
do
    coordinates=$( echo $line | cut -d ',' -f 1,2,3 | tr -d "chr" )

    if grep -q $coordinates $d-left-intersect   # If sequence overlaps a known protein-coding gene, remove
    then
        :
    elif grep -q $coordinates left-intersect    # If sequenec overlaps a known ncRNA gene, remove
    then
        :
    else
        count=$(($count+1))
        echo RNA$count,No,$line >> $d-negative-control-dataset.csv
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

for line in $( cat $d-right.csv )
do
    coordinates=$( echo $line | cut -d ',' -f 1,2,3 | tr -d "chr" )

    if grep -q $coordinates $d-right-intersect   # If sequence overlaps a known protein-coding gene, remove
    then
        :
    elif grep -q $coordinates right-intersect    # If sequenec overlaps a known ncRNA gene, remove
    then
        :
    else
        count=$(($count+1))
        echo RNA$count,No,$line >> $d-negative-control-dataset.csv
    fi
done

######## Generate negative control FASTA file
grep -v "Start" $d-negative-control-dataset.csv | cut -d ',' -f 1,6 | tr ',' ' ' | perl -lane '{print ">$F[0]\n$F[1]"}' > $d-negative-control-seq.fa

######## Report number of generated negative control sequences, post filtering
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
rm -rf true-coordinates
rm -rf $d-all-exon-coords
