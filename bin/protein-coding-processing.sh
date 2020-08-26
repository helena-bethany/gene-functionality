#!/bin/bash
#
# Script Name: protein-coding-processing.sh
#
# Author: Helena Cooper
# Last edited: 12/07/2020
#
# Description: This script filters the functional protein-coding genes and generates negative control sequences. 
#
# Input: $1 is the folder containing the additional files and dependencies not exported to $PATH (do not include '/' at the end)
#        $2 is the text file containing RefSeq IDs from HGNC
#        $3 is the GFF file for the downloaded human genome Eg: GCF_000001405.39_GRCh38.p13_genomic.gff
#        $4 is the converted CSV file of the human genome. Eg: GRCh38.p13_genome.csv
#        $5 is the uniprot bed file from UCSC for filtering against known protein-coding genes. Eg: uni-prot-ucsc.bed
#        $6 is the ncrna bed file from GENCODE for filtering against known ncRNAs. Eg: gencode-ncrna-annotation.bed
#

###########################################################################################################################

######## General setup

# set -xv
IFS=$'\n'
d=$(date +%y%m%d) 
calc() { awk "BEGIN{print $*}"; }

######## Delete files if they currently exist, as script continuously appends to the same file
rm -rf $d-left.csv
rm -rf $d-right.csv
rm -rf $d-left-coordinates.bed
rm -rf $d-left-intersect.bed
rm -rf $d-left-intersect
rm -rf $d-right-coordinates.bed
rm -rf $d-right-intersect.bed
rm -rf $d-right-intersect
rm -rf true-coordinates
rm -rf left-intersect.bed
rm -rf left-intersect
rm -rf right-intersect.bed
rm -rf right-intersect
rm -rf true-coordinates

###########################################################################################################################

# Obtain gene sequences and locations

###########################################################################################################################

######## Create list of randomly selected numbers, so sequences can be selected at random. 
max=$( cat $1/$2 | wc -l )
shuf -i 1-$max -n $max > numbers
count=1

######## Obtain chromosome coordinates for protein-coding exons
echo ID,Functional,Chromosome,Start,End,Sequence > $d-protein-exon2-dataset.csv
echo ID,Functional,Chromosome,Start,End,Sequence > $d-protein-exon3-dataset.csv

# Add sequence to summary csv file

for line in $( cat numbers )
do
    id=$( sed -n "$line"p $1/$2 )  # Choose a randomly selected ID

    grep "exon-$id" $1/$3 > info                  # Grep sequence information from GFF based on RefSeq ID
    coords_one=$( head -1 info | cut -f 1,4,5 )   # Required to generate the upstream negative control sequences
    coords_two=$( head -2 info | tail -1 | cut -f 1,4,5 )   # Exon two coordinates
    coords_three=$( head -3 info | tail -1 | cut -f 1,4,5 ) # Exon three coordinates
    final_end=$( tail -1 info | cut -f 1,4,5 )    # Required to generate the downstream negative control sequences.

    chr=$( echo $coords_one | cut -d ' ' -f 1 | tr -d "NC_" | cut -d '.' -f 1 | cut -c5,6 )   # Chromosome variable
    test=$( echo $chr | cut -c1 )   # Records any zeros in the chromosome variable
    other=$( echo $chr | cut -c2 )  # If zero is in chromosome variable, only record the single digit (ie: 01 becomes 1).
    mt_test=$( echo $coords_one | cut -d ' ' -f 1 )   # Variable to check if gene is located on the mitochondrial genome.
    
    # Reformat chr variable
    if [ -z "$chr" ]  # If chromosome variable empty (genes/mRNA that have been removed), then rename to allow it to be filtered out.
    then
        chr=26   
    elif [[ "$mt_test" == "NC_012920.1" ]]  # If gene is encoded on the mitochondrial genome, then rename to allow it to be filtered out.
    then
        chr=25     
    elif [[ "$test" == "0" ]] # If chromosome variable begins with zero, then rename as a single digit (ie: 01 becomes 1).
    then
        chr=$other
    elif [[ "$chr" == "23" ]]   # Chromosome X is NC_000023, but should be recorded as X in the final dataset for readability.
    then
        chr=X
    else
        :
    fi

    # Process exon data to create a dataset of 1000 sequences (max)
    if [ "$count" -gt "1000" ] 
    then
        :
    # Process if gene belongs to autosomal chromosome or chrX (ie: no unassembled scaffolds, Y or MT). 
    elif [ "$chr" -le "22" ] || [[ "$chr" == "X" ]]   
    then
        start_one=$( echo $coords_one | cut -f 2 )           # Start position of exon one
        end_one=$( echo $coords_one | cut -f 3 )             # End position of exon one
        start_two=$( echo $coords_two | cut -f 2 )           # Start position of exon two
        end_two=$( echo $coords_two | cut -f 3 )             # End position of exon two
        start_three=$( echo $coords_three | cut -f 2 )       # Start position of exon three
        end_three=$( echo $coords_three | cut -f 3 )         # End position of exon three
        seq_two=$( grep -w "chromosome $chr" $1/$4 | cut -f 2 | cut -c$start_two-$end_two )             # Exon two sequence
        seq_three=$( grep -w "chromosome $chr" $1/$4 | cut -f 2 | cut -c$start_three-$end_three )       # Exon three sequence
        end_final=$( echo $final_end | cut -f 3 )            # End position of final exon
	
        len_two=$(( $end_two - $start_two ))                 # Length of exon two
        len_three=$(( $end_three - $start_three ))           # Length of exon three

        # Included sequence must have a least two exons and not contained any unknown nucleotides (N)
        if [ ! -z "$seq_two" ] && [ ! -z "$seq_three" ] && [[ "$seq_two" != *"N"* ]] && [[ "$seq_three" != *"N"* ]]
        then
	    # Exclude sequences longer than 3000 nt
	    if [ "$len_two" -lt "3000" ] && [ "$len_three" -lt "3000" ] 
	    then
                echo RNA$count,Yes,chr$chr,$start_two,$end_two,$seq_two >> $d-protein-exon2-dataset.csv
                echo RNA$count,Yes,chr$chr,$start_three,$end_three,$seq_three >> $d-protein-exon3-dataset.csv
                if [ "$start_one" -gt "$end_final" ]   # Reverse transcripts can alter order of start/end positions
                then
		    # True coordinates file is used to generate negative control sequences that are the same length as exons two and three
                    echo chr$chr,$end_final,$start_one,$len_two,$len_three >> true-coordinates
                else
                    echo chr$chr,$start_one,$end_final,$len_two,$len_three >> true-coordinates
                fi
                count=$(( $count + 1 ))   # Counter for Sequence ID
            else
                :
	    fi
	else
	    :
	fi
     else
         :
     fi

done

rm -rf numbers

######## Generate functional FASTA file
for line in $( grep -v "Start" $d-protein-exon2-dataset.csv ) ; do echo $line | cut -d ',' -f 1,6 | tr ',' ' ' | perl -lane '{print ">$F[0]\n$F[1]"}' >> $d-protein-exon2-seq.fa ; done
for line in $( grep -v "Start" $d-protein-exon3-dataset.csv ) ; do echo $line | cut -d ',' -f 1,6 | tr ',' ' ' | perl -lane '{print ">$F[0]\n$F[1]"}' >> $d-protein-exon3-seq.fa ; done

###########################################################################################################################

# Create negative control protein-coding sequences from human genome

###########################################################################################################################

######## Generate sequences upstream and downstream of functional gene

for line in $( cat true-coordinates )
do
    #Extract chromosome coordinates and sequence length for functional sequence
    chr=$( echo $line | cut -d ',' -f 1 )
    chromo=$( echo $chr | tr -d "chr" )
    start=$( echo $line | cut -d ',' -f 2)
    end=$( echo $line | cut -d ',' -f 3 )
    left_length=$( echo $line | cut -d ',' -f 4 )   # null sequence still generated up or downstream of original coordinates.
    right_length=$( echo $line | cut -d ',' -f 5 )  # length is based on exons to match protein-coding exons.

    #Generate null sequence 20,000 upstream of sequence
    left_end=$(( $start - 20000 ))
    left_start=$(( $left_end - $length ))
    if [[ $left_end -lt 0 ]] || [[ $left_start -lt 0 ]] ; then left_sequence= ; else left_sequence=$( grep -w "chromosome $chromo" $1/$4 | cut -f 2 | cut -c$left_start-$left_end ) ; fi
   
    if [ -z $left_sequence ]   # If no sequence extracted, then remove.
    then
        :
    elif [[ $left_sequence == *"N"* ]]   # If sequence extracted contains unknown nucleotides (N), then remove.
    then
        :
    else
        echo $chr,$left_start,$left_end,$left_sequence >> $d-left.csv
    fi

    #Generate null sequence 20,000 downstream of sequence
    right_start=$(( $end + 20000 ))
    right_end=$(( $right_start + $length ))
    if [[ $right_start -lt 0 ]] || [[ $right_end -lt 0 ]] ; then right_sequence= ; else right_sequence=$( grep -w "chromosome $chromo" $1/$4 | cut -f 2 | cut -c$right_start-$right_end ) ; fi
    
    if [ -z $right_sequence ]     # If no sequence extracted, then remove.
    then
        :
    elif [[ $right_sequence == *"N"* ]]   # If sequence extracted contains unknown nucleotides (N), then remove.
    then
        :
    else
        echo $chr,$right_start,$right_end,$right_sequence >> $d-right.csv
    fi
done

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
count=$( tail -1 $d-protein-exon2-dataset.csv | cut -d ',' -f 1 | tr -d "RNA" )

for line in $( cat $d-left.csv )
do
    coordinates=$( echo $line | cut -d ',' -f 1,2,3 | tr -d "chr" )

    if grep -q $coordinates $d-left-intersect           # If overlap with protein-coding gene
    then
        :
    elif grep -q $coordinates left-intersect            # If overlap with ncRNA gene
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

    if grep -q $coordinates $d-right-intersect      # If overlap with protein-coding gene
    then
        :
    elif grep -q $coordinates right-intersect       # If overlap with ncRNA gene
    then
        :
    else
        count=$(($count+1))
        echo RNA$count,No,$line >> $d-negative-control-dataset.csv
    fi
done

######## Generate negative control FASTA file
grep -v "Start" $d-negative-control-dataset.csv | cut -d ',' -f 1,6 | tr ',' ' ' | perl -lane '{print ">$F[0]\n$F[1]"}' > $d-negative-control-seq.fa

###########################################################################################################################

######## Delete excess files
rm -rf numbers
rm -rf left-intersect.bed
rm -rf left-intersect
rm -rf $d-left-intersect.bed
rm -rf $d-left-intersect
rm -rf right-intersect.bed
rm -rf right-intersect
rm -rf $d-right-intersect.bed
rm -rf $d-right-intersect
rm -rf $d-left-coordinates.bed
rm -rf $d-right-coordinates.bed
rm -rf $d-left.csv
rm -rf $d-right.csv
