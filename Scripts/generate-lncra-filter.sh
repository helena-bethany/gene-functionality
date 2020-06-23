# Creation of the training/test dataset for the lncRNA (split by exons)

# set -xv

#####################################################################

# $1 is the downloaded fasta file from RNAcentral
# $2 is the downloaded chromosome coordinates file
# $3 is the csv file of the human genome
# $4 is the uniprot bed file from UCSC
# $5 is the ncRNA bed file from GENCODE
# $6 is the additional dependencies folder (optional)

IFS=$'\n'
d=$(date +%y%m%d)  # date identifier

#####################################################################

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

#####################################################################

# Convert FASTA file to CSV

fasta_formatter -i $1 -o $d-converted.csv -t

#####################################################################

# Obtain chromosome coordinates for each functional ncRNA

count=0

echo ID,Functional,Chromosome,Start,End,Sequence > $d-functional-ncrna-exon1-dataset.csv
echo ID,Functional,Chromosome,Start,End,Sequence > $d-functional-ncrna-exon2-dataset.csv

for line in $( cat $d-converted.csv )
do

    ID=$( echo $line | cut -f 1 | cut -c1-18 )
    sequence=$( echo $line | cut -f 2 )

    # Parse RNAcentral file
    if grep -q $ID $2
    then
        chr=$( grep -m 1 $ID $2 | cut -f 1 )
        start=$( grep -m 1 $ID $2 | cut -f 2 )
        true_end=$( grep -m 1 $ID $2 | cut -f 3 )
        len=$( grep -m 1 $ID $2 | cut -f 11 | cut -d ',' -f 1 )
        end=$(( $start + $len ))
        exon_count=$( grep -m 1 $ID $2 | cut -f 10 )

        if [ -z $sequence ]
        then
            :
	      elif [ $chr == 'chrM' ]
        then
            :
        elif [ $chr == 'chrY' ]
        then
            :
	      elif [ "$exon_count" -eq "1" ]
        then
            :
	      else
            count=$(( $count + 1 ))
            exon_one=$( echo $sequence | cut -c1-$len )
            echo $true_end,$chr,$start,$end,$exon_one >> ncrna-dataset-one

            start_exon=$( grep -m 1 $ID $2 | cut -f 12 | cut -d ',' -f 2 )
            start_two=$(( $start + $start_exon ))
            len_exon=$( grep -m 1 $ID $2 | cut -f 11 | cut -d ',' -f 2 )
            len=$(( $len + 1 ))
            end_two=$(( $start_two + $len_exon ))
            total=$(( $len + $len_exon ))
            exon_two=$( echo $sequence | cut -c$len-$total )
            echo $chr,$start_two,$end_two,$exon_two >> ncrna-dataset-two
        fi
    else
        :

    fi

done

if [ "$count" -gt "1000" ]
then

    max=$( cat ncrna-dataset-one | wc -l )
    shuf -i 1-$max -n 1000 > numbers
    id_count=1

    for line in $( cat numbers )
    do
      	info_one=$( sed -n "$line"p ncrna-dataset-one | cut -d ',' -f 2- )
        echo RNA$id_count,Yes,$info_one >> $d-functional-ncrna-exon1-dataset.csv
        info_two=$( sed -n "$line"p ncrna-dataset-two )
        echo RNA$id_count,Yes,$info_two >> $d-functional-ncrna-exon2-dataset.csv
        data=$( sed -n "$line"p ncrna-dataset-one | cut -d ',' -f 1 )
        other=$( echo $info_one | cut -d ',' -f 1,2 )
        echo $other,$data >> true-coordinates
        id_count=$(( $id_count + 1 ))
    done
else
    id_count=1
    for line in $( cat ncrna-dataset-one | cut -d ',' -f 2- ) ; do echo RNA$id_count,Yes,$line >> $d-functional-ncrna-exon1-dataset.csv ; id_count=$(( $id_count + 1 )) ; done

    id_count=1
    for line in $( cat ncrna-dataset-two) ; do echo RNA$id_count,Yes,$line >> $d-functional-ncrna-exon2-dataset.csv ; id_count=$(( $id_count + 1 )) ; done

    cat ncrna-dataset-one | cut -d ',' -f 2,3,4 >> true-coordinates
fi

grep -v "Start" $d-functional-ncrna-exon1-dataset.csv | cut -d ',' -f 1,6 | tr ',' ' ' | perl -lane '{print ">$F[0]\n$F[1]"}' > $d-functional-ncrna-exon1-seq.fa
grep -v "Start" $d-functional-ncrna-exon2-dataset.csv | cut -d ',' -f 1,6 | tr ',' ' ' | perl -lane '{print ">$F[0]\n$F[1]"}' > $d-functional-ncrna-exon2-seq.fa

other_count=$( grep -v "Start" $d-functional-ncrna-exon1-dataset.csv | wc -l )

echo
echo There were $other_count functional ncRNA found.
echo

#####################################################################

# Create negative control ncRNA from human genome

#text=$( grep -v "Start" $d-functional-ncrna-dataset.csv )

for line in $( cat true-coordinates )
do
    #Extract chromosome coordinates and sequence length for functional ncRNA
    chr=$( echo $line | cut -d ',' -f 1 )
    chromo=$( echo $chr | tr -d "chr" )
    start=$( echo $line | cut -d ',' -f 2)
    end=$( echo $line | cut -d ',' -f 3 )
    length=$( shuf -i 100-500 -n 1 )  # null sequence still generated up or downstream of original coordinates

    #Generate null sequence 20,000 upstream of ncRNA
    left_end=$(( $start - 20000 ))
    left_start=$(( $left_end - $length ))
    left_sequence=$( grep -w "chromosome $chromo" $3 | cut -f 2 | cut -c$left_start-$left_end )
    if [ -z $left_sequence ]
    then
        :
    elif [[ $left_sequence == *"N"* ]]
    then
        :
    else
        echo $chr,$left_start,$left_end,$left_sequence >> $d-left.csv
    fi
    #Generate null sequence 20,000 downstream of ncRNA
    right_start=$(( $end + 20000 ))
    right_end=$(( $right_start + $length ))
    right_sequence=$( grep -w "chromosome $chromo" $3 | cut -f 2 | cut -c$right_start-$right_end )
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

echo Negative control ncRNA sequences generated
echo

#####################################################################

# Filter left negative control sequences using UniProt and GENCODE

grep -v "Start" $d-left.csv | cut -d ',' -f 1,2,3 | tr ',' ' ' | perl -lane '{print "$F[0] $F[1] $F[2]"}' | tr ' ' '\t' > $d-left-coordinates.bed
[[ ! -z "$6" && -f $6/bedtools ]] && $6/bedtools intersect -a $4 -b $d-left-coordinates.bed > $d-left-intersect.bed || bedtools intersect -a $4 -b $d-left-coordinates.bed > $d-left-intersect.bed
cat $d-left-intersect.bed | cut -f 1,2,3 | tr -d "chr" | tr '\t' ',' > $d-left-intersect
left_count=$( cat $d-left-intersect | wc -l )

[[ ! -z "$6" && -f $6/bedtools ]] && $6/bedtools intersect -a $5 -b $d-left-coordinates.bed > left-intersect.bed || bedtools intersect -a $5 -b $d-left-coordinates.bed > left-intersect.bed
cat left-intersect.bed | cut -f 1,2,3 | tr -d "chr" | tr '\t' ',' > left-intersect

echo ID,Functional,Chromosome,Start,End,Sequence > $d-negative-control-dataset.csv
count=$( tail -1 $d-functional-ncrna-exon1-dataset.csv | cut -d ',' -f 1 | tr -d "RNA" )


for line in $( cat $d-left.csv )
do

    coordinates=$( echo $line | cut -d ',' -f 1,2,3 | tr -d "chr" )

    if grep -q $coordinates $d-left-intersect
    then
        :
    elif grep -q $coordinates left-intersect
    then
        :
    else
        count=$(($count+1))
        echo RNA$count,No,$line >> $d-negative-control-dataset.csv
    fi

done

#####################################################################

# Filter right negative control sequences using UniProt and GENCODE

grep -v "Start" $d-right.csv | cut -d ',' -f 1,2,3 | tr ',' ' ' | perl -lane '{print "$F[0] $F[1] $F[2]"}' | tr ' ' '\t' > $d-right-coordinates.bed
[[ ! -z "$6" && -f $6/bedtools ]] && $6/bedtools intersect -a $4 -b $d-right-coordinates.bed > $d-right-intersect.bed || bedtools intersect -a $4 -b $d-right-coordinates.bed > $d-right-intersect.bed
cat $d-right-intersect.bed | cut -f 1,2,3 | tr -d "chr" | tr '\t' ',' > $d-right-intersect
right_count=$( cat $d-right-intersect | wc -l )

[[ ! -z "$6" && -f $6/bedtools ]] && $6/bedtools intersect -a $5 -b $d-right-coordinates.bed > right-intersect.bed || bedtools intersect -a $5 -b $d-right-coordinates.bed > right-intersect.bed
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

total=$(( $left_count + $right_count ))
echo There were $total negative control sequences that overlapped with annotated coding regions.
echo
file_count=$( cat $d-negative-control-dataset.csv | wc -l )
count=$(( $file_count - 1 ))
echo There were $count non-overlapping negative control sequences generated.
echo

#############################################################################

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
