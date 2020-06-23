# Obtain sequence from NCBI using RefSeq IDs

#set -xv

rm -rf numbers
rm -rf genes.fa
rm -rf seq
rm -rf summary
rm -rf interim
rm -rf formatted-fasta.fa
rm -rf true-coordinates

# $1 is the text file with the RefSeq IDs
# $2 is the human genome GFF file
# $3 is the dependencies location
# $4 is the CSV human genome
# $5 is the uniprot file
# $6 is the GENCODE file

####################################################

# Script set up

IFS=$'\n'
d=$(date +%y%m%d)
calc() { awk "BEGIN{print $*}"; }

rm -rf $d-protein-coding-genes.fa

####################################################

# Obtain gene sequences and locations

max=$( cat $1 | wc -l )
shuf -i 1-$max -n 1000 > numbers
count=1

echo ID,Functional,Chromosome,Start,End,Sequence > $d-protein-exon1-dataset.csv
echo ID,Functional,Chromosome,Start,End,Sequence > $d-protein-exon2-dataset.csv

# Add sequence to summary csv file

for line in $( cat numbers )
do

    id=$( sed -n "$line"p $1 )

    grep "exon-$id" $2 > info
    coords_one=$( head -1 info | cut -f 1,4,5 )
    coords_two=$( head -2 info | tail -1 | cut -f 1,4,5 )
    final_end=$( tail -1 info | cut -f 1,4,5 )

    chr=$( echo $coords | cut -d ' ' -f 1 | tr -d "NC_" | cut -d '.' -f 1 | cut -c5,6 )
    test=$( echo $chr | cut -c1 )
    other=$( echo $chr | cut -c2 )

    if [ "$test" -eq "0" ]    # Remove zero before each chromosome number (ie: 01 becomes 1)
    then
        chr=$other
    else
        :
    fi

    # Only process sequence if from Chr1-22 or ChrX
    if [ "$chr" -le "22" ] || [[ "$chr" == "X" ]]
    then
        start_one=$( echo $coords_one | cut -f 2 )
        end_one=$( echo $coords_one | cut -f 3 )
        start_two=$( echo $coords_two | cut -f 2 )
        end_two=$( echo $coords_two | cut -f 3 )
        seq_one=$( grep -w "chromosome $chr" $4 | cut -f 2 | cut -c$start_one-$end_one )
        seq_two=$( grep -w "chromosome $chr" $4 | cut -f 2 | cut -c$start_two-$end_two )
        end_final=$( echo $final_end | cut -f 3 )

        if [ ! -z "$seq_one" ] && [ ! -z "$seq_two" ] # Included sequence must have a least two exons
        then
            echo RNA$count,Yes,chr$chr,$start_one,$end_one,$seq_one >> $d-protein-exon1-dataset.csv
            echo RNA$count,Yes,chr$chr,$start_two,$end_two,$seq_two >> $d-protein-exon2-dataset.csv
            if [ "$start_one" -gt "$end_final" ]   # Reverse transcripts can alter order of start/end positions
            then
                echo chr$chr,$end_final,$start_one >> true-coordinates
            else
                echo chr$chr,$start_one,$end_final >> true-coordinates
            fi
            count=$(( $count + 1 ))
        else
            :
	    fi
     else
         :
     fi

done

# Conversion of fasta with '\n' to without
# cat genes.fa | awk '/^>/ {printf("\n%s\n",$1);next; } { printf("%s",$1);}  END {printf("\n");}' | tail -"$num" > formatted-fasta.fa

#####################################################

# Generate FASTA files
grep -v "Start" $d-protein-exon1-dataset.csv | cut -d ',' -f 1,6 | tr ',' ' ' | perl -lane '{print ">$F[0]\n$F[1]"}' > $d-protein-exon1-seq.fa
grep -v "Start" $d-protein-exon2-dataset.csv | cut -d ',' -f 1,6 | tr ',' ' ' | perl -lane '{print ">$F[0]\n$F[1]"}' > $d-protein-exon2-seq.fa

######################################################

rm -rf numbers
rm -rf genes.fa
rm -rf seq
rm -rf summary
rm -rf interim
rm -rf formatted-fasta.fa
rm -rf rna_id

#####################################################

# Generate negative control sequences from human genome

for line in $( cat true-coordinates )
do
    #Extract chromosome coordinates and sequence length for functional sequence
    chr=$( echo $line | cut -d ',' -f 1 )
    chromo=$( echo $chr | tr -d "chr" )
    start=$( echo $line | cut -d ',' -f 2)
    end=$( echo $line | cut -d ',' -f 3 )
    length=$( shuf -i 100-1000 -n 1 )  # null sequence still generated up or downstream of original coordinates

    #Generate null sequence 20,000 upstream of sequence
    left_end=$(( $start - 20000 ))
    left_start=$(( $left_end - $length ))
    left_sequence=$( grep -w "chromosome $chromo" $4 | cut -f 2 | cut -c$left_start-$left_end )
    if [ -z $left_sequence ]
    then
        :
    elif [[ $left_sequence == *"N"* ]]
    then
        :
    else
        echo $chr,$left_start,$left_end,$left_sequence >> $d-left.csv
    fi

    #Generate null sequence 20,000 downstream of sequence
    right_start=$(( $end + 20000 ))
    right_end=$(( $right_start + $length ))
    right_sequence=$( grep -w "chromosome $chromo" $4 | cut -f 2 | cut -c$right_start-$right_end )
    if [ -z $right_sequence ]
    then
        :
    elif [[ $right_sequence == *"N"* ]]
    then
        :
    else
        echo $chr,$right_start,$right_end,$right_sequence >> $d-right.csv
    fi
done

#####################################################################

# Filter left negative control sequences using UniProt (protein) and GENCODE (ncRNA)

grep -v "Start" $d-left.csv | cut -d ',' -f 1,2,3 | tr ',' ' ' | perl -lane '{print "$F[0] $F[1] $F[2]"}' | tr ' ' '\t' > $d-left-coordinates.bed
[ -f $3/bedtools ] && $3/bedtools intersect -a $5 -b $d-left-coordinates.bed > $d-left-intersect.bed || bedtools intersect -a $5 -b $d-left-coordinates.bed > $d-left-intersect.bed
cat $d-left-intersect.bed | cut -f 1,2,3 | tr -d "chr" | tr '\t' ',' > $d-left-intersect
left_count=$( cat $d-left-intersect | wc -l )

[ -f $3/bedtools ] && $3/bedtools intersect -a $6 -b $d-left-coordinates.bed > left-intersect.bed || bedtools intersect -a $6 -b $d-left-coordinates.bed > left-intersect.bed
cat left-intersect.bed | cut -f 1,2,3 | tr -d "chr" | tr '\t' ',' > left-intersect

echo ID,Functional,Chromosome,Start,End,Sequence > $d-negative-control-dataset.csv
count=$( tail -1 $d-protein-exon1-dataset.csv | cut -d ',' -f 1 | tr -d "RNA" )

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

#####################################################################

# Filter right negative control sequences using UniProt

grep -v "Start" $d-right.csv | cut -d ',' -f 1,2,3 | tr ',' ' ' | perl -lane '{print "$F[0] $F[1] $F[2]"}' | tr ' ' '\t' > $d-right-coordinates.bed
[ -f $3/bedtools ] && $3/bedtools intersect -a $5 -b $d-right-coordinates.bed > $d-right-intersect.bed || bedtools intersect -a $5 -b $d-right-coordinates.bed > $d-right-intersect.bed
cat $d-right-intersect.bed | cut -f 1,2,3 | tr -d "chr" | tr '\t' ',' > $d-right-intersect
right_count=$( cat $d-right-intersect | wc -l )

[ -f $3/bedtools ] && $3/bedtools intersect -a $6 -b $d-right-coordinates.bed > right-intersect.bed || bedtools intersect -a $6 -b $d-right-coordinates.bed > right-intersect.bed
cat right-intersect.bed | cut -f 1,2,3 | tr -d "chr" | tr '\t' ',' > right-intersect

for line in $( cat $d-right.csv )
do

    coordinates=$( echo $line | cut -d ',' -f 1,2,3 | tr -d "chr" )

    if grep -q $coordinates $d-right-intersect
    then
        :
    elif grep -q $coordinates right-intersect
    then
        :
    else
        count=$(($count+1))
        echo RNA$count,No,$line >> $d-negative-control-dataset.csv
    fi

done

grep -v "Start" $d-negative-control-dataset.csv | cut -d ',' -f 1,6 | tr ',' ' ' | perl -lane '{print ">$F[0]\n$F[1]"}' > $d-negative-control-seq.fa

#############################################################################

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
