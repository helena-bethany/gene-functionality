# Creation of the training/test dataset

#####################################################################

# $1 is the downloaded fasta file from RNAcentral
# $2 is the downloaded chromosome coordinates file
# $3 is the csv file of the human genome

IFS=$'\n'
d=$(date +%y%m%d)  # date identifier
mkdir $d-ncrna-excess # location of additional files generated

#####################################################################

# Convert FASTA file to CSV

fasta_formatter -i $1 -o $d-converted.csv -t

#####################################################################

# Obtain chromosome coordinates for each functional ncRNA

count=0

echo ncRNA,RNAcentral ID,Chromosome,Start,End,Sequence > $d-coords.csv

for line in $(cat $d-converted.csv)
do

    ID=$( echo $line | cut -f 1 | cut -c1-18 )
    sequence=$( echo $line | cut -f 2 )
    info=$( echo $line | cut -f 1 | cut -d "(" -f 2- | cut -d "(" -f 2 | tr -d ")" | cut -d ',' -f 1 )

    if grep -q $ID $2
    then
	chr=$( grep -m 1 $ID $2 | cut -f 1 | cut -c4- )
	start=$( grep -m 1 $ID $2 | cut -f 2 )
	end=$( grep -m 1 $ID $2 | cut -f 3 )
	count=$(( $count + 1 ))

	echo $info,$ID,$chr,$start,$end,$sequence >> $d-coords.csv
    
	echo ">$ID" >> $d-ncrna.fa
	echo $sequence >> $d-ncrna.fa

    else
	:

    fi

done

echo There were $count coordinates found.
echo

#####################################################################

# Run cmscan to calculate bit score from pre-defined Rfam CMs

cmscan --cpu 20 --rfam --cut_ga --nohmmonly --tblout $d.tblout --fmt 2 --clanin Rfam.12.1.clanin Rfam.cm $d-ncrna.fa %> /dev/null
echo cmscan complete
echo

#####################################################################

# Extract top cmscan hit and add metrics to one file

file=$( grep -v "#" $d.tblout )
echo Chromosome,Start,End,CM-scan,eValue-scan,Sequence > $d-final.csv

for line in $file
do

    rank=$(echo $line | tr -s ' ' | cut -d ' ' -f 1 )
    
    if (( $(echo "$rank == 1" | bc -l) ))
    then
	#Extract relevant information from cmscan output
	rna=$(echo $line | tr -s ' ' | cut -d ' ' -f 2 )
	id=$(echo $line | tr -s ' ' | cut -d ' ' -f 4 )
	cm=$(echo $line | tr -s ' ' | cut -d ' ' -f 17 )
	evalue=$(echo $line | tr -s ' ' | cut -d ' ' -f 18 )

	#Extract corresponding information from coords output
	chr=$( grep $id $d-coords.csv | cut -d ',' -f 3 )
	start=$( grep $id $d-coords.csv | cut -d ',' -f 4 )
	end=$( grep $id $d-coords.csv | cut -d ',' -f 5 )
	sequence=$( grep $id $d-coords.csv | cut -d ',' -f 6 )

	#Add information to _final.csv output
	echo $chr,$start,$end,$cm,$evalue,$sequence >> $d-final.csv
	
    else
	:
    fi

done

echo Functional ncRNA data generated    
echo

#####################################################################

# Create negative control ncRNA from human genome

text=$( grep -v "Start" $d-final.csv )

for line in $text
do
    #Extract chromosome coordinates and sequence length for functional ncRNA
    chr=$( echo $line | cut -d ',' -f 1 )
    start=$( echo $line | cut -d ',' -f 2)
    end=$( echo $line | cut -d ',' -f 3 )
    length=$(( $end - $start ))
    
    #Generate null sequence 20,000 upstream of ncRNA
    left_end=$(( $start - 20000 ))
    left_start=$(( $left_end - $length ))
    left_sequence=$( grep -w "chromosome $chr" $3 | cut -f 2 | cut -c$left_start-$left_end )
    echo $chr,$left_start,$left_end,$left_sequence >> $d-left.csv
    
    #Generate null sequence 20,000 downstream of ncRNA
    right_start=$(( $end + 20000 ))
    right_end=$(( $right_start + $length ))
    right_sequence=$( grep -w "chromosome $chr" $3 | cut -f 2 | cut -c$right_start-$right_end )
    echo $chr,$right_start,$right_end,$right_sequence >> $d-right.csv
done

echo Negative control ncRNA sequences generated

#####################################################################

# Run cmscan to calculate bit score from pre-defined Rfam CMs for upstream sequences

grep -v 'Sequence' $d-left.csv | cut -d ',' -f 2,4 | tr ',' ' ' | perl -lane '{print ">$F[0]\n$F[1]"}' > $d-left.fa
cmscan --cpu 20 --rfam -E 10000 --nohmmonly --tblout $d-left.tblout --fmt 2 --clanin Rfam.12.1.clanin Rfam.cm $d-left.fa &> /dev/null
echo cmscan complete
echo

file=$( grep -v "#" $d-left.tblout )

for line in $file
do

    rank=$(echo $line | tr -s ' ' | cut -d ' ' -f 1 )
    
    if (( $(echo "$rank == 1" | bc -l) ))
    then
	#Extract relevant information from cmscan output
	rna_left=$(echo $line | tr -s ' ' | cut -d ' ' -f 2 )
	id_left=$(echo $line | tr -s ' ' | cut -d ' ' -f 4 )
	cm_left=$(echo $line | tr -s ' ' | cut -d ' ' -f 17 )
	evalue_left=$(echo $line | tr -s ' ' | cut -d ' ' -f 18 )
    
	#Extract corresponding information from coords output
	chr_left=$( grep $id $d-coords.csv | cut -d ',' -f 1 | cut -c4-)
	start_left=$( grep $id $d-coords.csv | cut -d ',' -f 2 )
	end_left=$( grep $id $d-coords.csv | cut -d ',' -f 3 )
	sequence_left=$( grep $id $d-coords.csv | cut -d ',' -f 4 )

	#Add information to final.csv output
	echo $chr_left,$start_left,$end_left,$cm_left,$evalue_left,$sequence_left >> $d-final-left.csv 
	
    else
	:
    fi

done

#####################################################################

# Run cmscan to calculate bit score from pre-defined Rfam CMs for downstream sequences

grep -v 'Sequence' $d-right.csv | cut -d ',' -f 2,4 | tr ',' ' ' | perl -lane '{print ">$F[0]\n$F[1]"}' > $d-right.fa
cmscan --cpu 20 --rfam -E 10000 --nohmmonly --tblout $d-right.tblout --fmt 2 --clanin Rfam.12.1.clanin Rfam.cm $d-right.fa &> /dev/null
echo cmscan complete
echo

file=$( grep -v "#" $d-right.tblout )

for line in $file
do

    rank=$(echo $line | tr -s ' ' | cut -d ' ' -f 1 )
    
    if (( $(echo "$rank == 1" | bc -l) ))
    then
	#Extract relevant information from cmscan output
	rna_right=$(echo $line | tr -s ' ' | cut -d ' ' -f 2 )
	id_right=$(echo $line | tr -s ' ' | cut -d ' ' -f 4 )
	cm_right=$(echo $line | tr -s ' ' | cut -d ' ' -f 17 )
	evalue_right=$(echo $line | tr -s ' ' | cut -d ' ' -f 18 )
    
	#Extract corresponding information from coords output
	chr_right=$( grep $id $d-coords.csv | cut -d ',' -f 1 | cut -c4-)
	start_right=$( grep $id $d-coords.csv | cut -d ',' -f 2 )
	end_right=$( grep $id $d-coords.csv | cut -d ',' -f 3 )
	sequence_right=$( grep $id $d-coords.csv | cut -d ',' -f 4 )

	#Add information to final.csv output
	echo $chr_right,$start_right,$end_right,$cm_right,$evalue_right,$sequence_right >> $d-final-right.csv 
	
    else
	:
    fi

done

#####################################################################

# Combine ncRNA together and add IDs/Functionality status

cat $d-final.csv $d-final-left.csv $d-final-right.csv > $d-combined.csv

echo Functionality > status.txt
samples=$(grep -v "Start" $d-final.csv | wc -l)
left=$(grep -v "Start" $d-final-left.csv | wc -l)
right=$(grep -v "Start" $d-final-right.csv | wc -l)
total=$(($left+$right))
y=0
n=0
while (($y != $(($samples))))
do
    echo Yes >> status.txt
    y=$(($y+1))
done
while (($n != $total))
do
    echo No >> status.txt
    n=$(($n+1))
done

paste -d ',' status.txt $d-combined.csv > $d-combined-2.csv

grep -v [XYM], $d-combined-2.csv > $d-intermediate
grep -v ' ' $d-intermediate > $d-to-add.csv

echo Name > labels.txt
count=1
all=$(grep -v "Start" $d-to-add.csv | wc -l)
while (($count <= $all))
do
    echo RNA$count >> labels.txt
    count=$(($count+1))
done

paste -d ',' labels.txt $d-to-add.csv > $d-ncrna-dataset-1.csv
echo $d-ncrna-dataset-1.csv completed
echo

#####################################################################

# Move excess files to dated folder

mv $d-converted.csv $d-ncrna-excess/
mv $d-coords.csv $d-ncrna-excess/
mv $d-ncrna.fa $d-ncrna-excess/
mv $d.tblout $d-ncrna-excess/
mv $d-final.csv $d-ncrna-excess/
mv $d-left.csv $d-ncrna-excess/
mv $d-right.csv $d-ncrna-excess/
mv $d-left.fa $d-ncrna-excess/
mv $d-left.tblout $d-ncrna-excess/
mv $d-right.fa $d-ncrna-excess/
mv $d-right.tblout $d-ncrna-excess/
mv $d-final-left.csv $d-ncrna-excess/
mv $d-final-right.csv $d-ncrna-excess/
mv $d-combined.csv $d-ncrna-excess/

rm -rf labels.txt
rm -rf status.txt
rm -rf $d-combined-2.csv
rm -rf $d-intermediate
rm -rf $d-to-add.csv

#####################################################################

# Obtain correctly formatted chromosome positions for UCSC Table Browser

grep -v 'Chromosome' $d-ncrna-dataset.csv | cut -d ',' -f 3,4,5 | tr ',' ' ' | perl -lane '{print "$F[0]:$F[1]-$F[2]"}' > $d-ncrna-dataset-ucsc.txt

echo RNAcentral.sh has finished running

