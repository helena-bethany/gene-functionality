# Calculating ncRNA secondary structure metrics

#####################################################################

# $1 is the input ncRNA dataset file

IFS=$'\n'
date=$(date +%y%m%d) # date identifier
id=$( echo $1 | cut -d '.' -f 1 | cut -d '-' -f 4 )
id_count=$(( $id + 1 ))

rm -rf $date-rscapedata
mkdir $date-rscapedata

var=1

#####################################################################

# Define biopython function to convert MAF to Stockholm

maf_to_stockholm() {

PYTHON_ARG="$1" python3 - <<END

from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys

a = []
RNA = str(os.environ['PYTHON_ARG'])

alignments = AlignIO.parse("mafOut", 'maf')
for alignment in alignments:
	a.append(alignment)	# Extracts each alignment from MAF file

my_records = []
file = open("RNA.stk", 'w')
file.write('# STOCKHOLM 1.0\n')	# Setting up beginning of stk file
file.write('\n')
file.write('#=GF ID '+RNA)
file.write('\n')
file.write('\n')

alignment = a[0]		
for record in alignment:
	n = 1
	ID = record.id		# Extracts ID from alignment information
	seq = str(record.seq)	# Extracts sequence from alignment information
	copies = 0
	while n < len(a):
		other = a[n]
		for record in other:
			d = record.id
			if d == ID:
				seq += str(record.seq)	# MAF file kept the sequence split on numerous lines
				copies += 1		# so this connects them together into one sequence.
			else:
				pass
		n += 1
	ID = ID.replace('.','/')
	if copies != (len(a)-1):
		pass
	else:
		file.write(ID)
		file.write(' '*(40-len(ID)))	# Making sure than spacing is the same for all sequences
		file.write(seq)
		file.write('\n')

file.close()

count = SeqIO.convert("RNA.stk", 'stockholm', RNA+".fa", "fasta")

END

}


#####################################################################

# Obtain multiple alignment files and calculate covariance

grep -v 'Chromosome' $1 | cut -d ',' -f 3,4,5 | tr ',' ' ' | perl -lane '{print "Chr$F[0] $F[1] $F[2]"}' > coordinates

for line in  $(cat coordinates)
do
	echo $line > overBed
	./mafFetch hg38 multiz100way overBed mafOut
	name='RNA'$var
	maf_to_stockholm $name
	var=$((var+1))
	./RNAalifold -f S --aln-stk=$name RNA.stk &> /dev/null
	./RNAalifold -f S --aln-stk=final RNA.stk &> /dev/null
	rscape_v1.0.4/bin/R-scape -E 100 -s $name.stk &> /dev/null
	echo $name
	rm -rf $name.sorted.out
	mv $name.out \$date-rscapedata
	rm -rf overBed
	rm -rf blank.txt
	rm -rf *.pdf
	rm -rf *ss.ps
	rm -rf *.svg
	rm -rf *.surv
	rm -rf *.sum
	rm -rf *.sto
	rm -rf *.fa
done

#####################################################################

# Build covariance model and calculate CM-HMM

cmbuild RNAfinal.cm final.stk > $date-cmbuild.txt
echo nseq, effnseq, CM-build, HMM-build > $date-cmbuildfinal.txt
grep -v '#' $date-cmbuild.txt | tr -s " " | cut -d ' ' -f 3,4,5,10,11 | tr ' ' ',' >> $date-cmbuild1.txt
var=1
for line in $(cat $date-cmbuild1.txt)
do

    echo $line > snp
    r=$(cat snp | cut -d ',' -f 1)
    data=$(cat snp | cut -d ',' -f 2,3,4,5)
    rna=RNA$var
    if (( $r == $rna ))
    then
	echo $data >> $date-cmbuildfinal.txt
    else
	echo 0,0,0,0 >> $date-cmbuildfinal.txt
    fi
    var=$(($var+1))
done

echo CM-HMM > $date-cmhmm.txt
calc() { awk "BEGIN{print $*}"; }

for line in $(cat $date-cmbuildfinal.txt | grep -v 'CM')
do
	echo $line > file.txt
	CM=$(cut -d ',' -f 3 file.txt)
	HMM=$(cut -d ',' -f 4 file.txt)
	diff=$(calc $CM-$HMM)
	echo $diff >> $date-cmhmm.txt
done

rm -rf file.txt

#####################################################################

# Process R-scape results

echo covMin10,covMax10,averageCov,#sigpair,#compatiable,#incompatiable,totalpairs > $date-rscape.csv
FILES=./$date-rscapedata/*.out
data=$(ls $FILES | sort -V)

for f in $data
do
	min=$(grep -r 'GTp' $f | cut -d '[' -f 2 | tr ']' ',') #>> $1cov.csv
   	a=$(grep -r '*' $f | wc -l)   # number of *
    	b=$(grep -r '~' $f | wc -l)   # number of ~
    	c=$(grep -v '#' $f | wc -l)   # total hits
    	d=$(($c-$a-$b))               # number of ' '
   	grep -v '#' $f | cut -f 4 > score
    	total=0
    	count=0
    	for number in $(cat score)
    	do
		total=$(calc $total+$number)
		count=$(($count+1))
    	done

    	if (( $(echo "$count == 0" | bc -l) ))
    	then
		count=1
	else
		:
	fi

    	average=$(calc $total/$count)

	if (( $(echo "$average == 0" | bc -l) ))
	then
	    echo 0, 0, 0, 0, 0, 0, 0 > info
	else
	    echo $min $average, $a, $b, $d, $c > info
	fi
	
    	cat info >> $date-rscape.csv
done

paste -d ',' $1 $date-rscape.csv $date-cmbuildfinal.txt $date-cmhmm.txt > $date-ncrna-dataset-$id_count.csv

#####################################################################

rm -rf *.stk
rm -rf coordinates
rm -rf overBed
rm -rf final.stk
rm -rf RNAfinal.cm
