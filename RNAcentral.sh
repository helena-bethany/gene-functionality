#Converting the data downloaded from RNAcentral all the way to the end product

# $1 is the downloaded fasta file from RNAcentral
# $2 is the date the dataset was generated, or any other identifiers

#Format Dataset

fasta_formatter -i $1 -o converted.csv -t
python3 coordinates2.py converted.csv $2coords.csv $2.fa
echo
#cmscan --rfam -E 10000 --nohmmonly --tblout $2.tblout --fmt 2 --clanin /home/helenacooper/infernal-1.1.2/testsuite/Rfam.12.1.clanin /home/helenacooper/infernal-1.1.2/Rfam.cm $2.fa &> /dev/null
cmscan --rfam --cut_ga --nohmmonly --tblout $2.tblout --fmt 2 --clanin /home/helenacooper/infernal-1.1.2/testsuite/Rfam.12.1.clanin /home/helenacooper/infernal-1.1.2/Rfam.cm $2.fa &> /dev/null
echo cmscan complete
echo
python3 RfamCM.py $2.tblout $2coords.csv $2_final.csv Test
echo

#Create Null Datasets

python3 random_sequences.py $2_final.csv $2_left.csv $2_right.csv
echo

#Obtain CM scores for Null Datasets

grep -v 'Sequence' $2_left.csv | cut -d ',' -f 2,4 | tr ',' ' ' | perl -lane '{print ">$F[0]\n$F[1]"}' > $2_left.fa
cmscan --rfam -E 10000 --nohmmonly --tblout $2_left.tblout --fmt 2 --clanin /home/helenacooper/infernal-1.1.2/testsuite/Rfam.12.1.clanin /home/helenacooper/infernal-1.1.2/Rfam.cm $2_left.fa &> /dev/null
python3 RfamCM.py $2_left.tblout $2_left.csv $2_Left.csv Null
echo
#(paste -d',' <(cut -d',' -f-3 $2_left.csv) text.file <(cut -d',' -f4- $2_left.csv)) > $2_Left.csv
#rm -rf $2_left.csv
#rm -rf text.file
grep -v 'Sequence' $2_right.csv | cut -d ',' -f 2,4 | tr ',' ' ' | perl -lane '{print ">$F[0]\n$F[1]"}' > $2_right.fa
cmscan --rfam -E 10000 --nohmmonly --tblout $2_right.tblout --fmt 2 --clanin /home/helenacooper/infernal-1.1.2/testsuite/Rfam.12.1.clanin /home/helenacooper/infernal-1.1.2/Rfam.cm $2_right.fa &> /dev/null
python3 RfamCM.py $2_right.tblout $2_right.csv $2_Right.csv Null
echo
#(paste -d',' <(cut -d',' -f-3 $2_right.csv) text.file <(cut -d',' -f4- $2_right.csv)) > $2_Right.csv
#rm -rf $2_right.csv
#rm -rf text.file
echo cmscan complete
echo

#Combine the Datasets together and add IDs/Functionality status

cat $2_final.csv $2_Left.csv $2_Right.csv > $2combined.csv

echo Functionality > status.txt
samples=$(grep -v "Start" $2_final.csv | wc -l)
IFS=$'\n'
left=$(grep -v "Start" $2_Left.csv | wc -l)
right=$(grep -v "Start" $2_Right.csv | wc -l)
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

#echo $samples $y
#echo $total $n
#echo

echo Name > labels.txt
count=1
all=$(($total+$samples+1))
while (($count != $all))
do
    echo RNA$count >> labels.txt
    count=$(($count+1))
done

#echo $all $count
#echo

(paste -d ',' labels.txt status.txt $2combined.csv) > $2.csv
echo Final dataset completed
echo

#Obtain Chromosome positions for UCSC Table Browser

grep -v 'Chromosome' $2.csv | cut -d ',' -f 3,4,5 | tr ',' ' ' | perl -lane '{print "$F[0]:$F[1]-$F[2]"}' > $2finalnames.txt

mv $2.csv ..
mv $2finalnames.txt ..

echo RNAcentral.sh has finished running

