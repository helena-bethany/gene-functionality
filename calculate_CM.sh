#Calculating CM scores for a tonne of sequences

#grep -v 'Chromosome' $1 | cut -d ',' -f 2,3,4 | tr ',' ' ' | perl -lane '{print "Chr$F[0] $F[1] $F[2]"}' > coordinates
grep -v 'Chromosome' $1 | cut -d ',' -f 3,4,5 | tr ',' ' ' | perl -lane '{print "Chr$F[0] $F[1] $F[2]"}' > coordinates

rm -rf rscapedata
mkdir rscapedata

IFS=$'\n'
rm -rf overBed
rm -rf *.stk
rm -rf final.stk
rm -rf RNAfinal.cm
var=1

for line in  $(cat coordinates)
do
	echo $line > overBed
	./mafFetch hg38 multiz100way overBed mafOut
	#grep -r "score" mafOut | wc -l
	name='RNA'$var
	python3 conversion3.py mafOut $name
	#tr '.' '/' RNA.stk
	#mv $name.fa \FASTA2
	var=$((var+1))
	./RNAalifold -f S --aln-stk=$name RNA.stk &> /dev/null
	./RNAalifold -f S --aln-stk=final RNA.stk &> /dev/null
	rscape_v1.0.4/bin/R-scape -E 100 -s $name.stk &> /dev/null
	echo $name
	rm -rf $name.sorted.out
	mv $name.out \rscapedata
	#mv $name.txt \rscapedata
	rm -rf overBed
	#rm -rf RNA.stk
	rm -rf blank.txt
	rm -rf *.pdf
	rm -rf *ss.ps
	rm -rf *.svg
	rm -rf *.surv
	rm -rf *.sum
	rm -rf *.sto
	rm -rf *.fa
done

cmbuild RNAfinal.cm final.stk > $2cmbuild.txt
echo nseq, effnseq, CM-build, HMM-build > $2cmbuildfinal.txt
grep -v '#' $2cmbuild.txt | tr -s " " | cut -d ' ' -f 3,4,5,10,11 | tr ' ' ',' >> $2cmbuild1.txt
var=1
for line in $(cat $2cmbuild1.txt)
do

    echo $line > snp
    r=$(cat snp | cut -d ',' -f 1)
    data=$(cat snp | cut -d ',' -f 2,3,4,5)
    rna=RNA$var
    if (( $r == $rna ))
    then
	echo $data >> $2cmbuildfinal.txt
    else
	echo 0,0,0,0 >> $2cmbuildfinal.txt
    fi
    var=$(($var+1))
done

echo CM-HMM > $2cmhmm.txt
calc() { awk "BEGIN{print $*}"; }

for line in $(cat $2cmbuildfinal.txt | grep -v 'CM')
do
	echo $line > file.txt
	CM=$(cut -d ',' -f 3 file.txt)
	HMM=$(cut -d ',' -f 4 file.txt)
	diff=$(calc $CM-$HMM)
	echo $diff >> $2cmhmm.txt
done

rm -rf file.txt

echo covMin10,covMax10,#sigpair,#compatiable,#incompatiable,totalpairs,averageCov > $2rscape.csv
FILES=/home/helenacooper/rscapedata/*.out
data=$(ls $FILES | sort -V)

for f in $data
do
	min=$(grep -r 'GTp' $f | cut -d '[' -f 2 | tr ']' ',') #>> $1cov.csv
   	a=$(grep -r '*' $f | wc -l)   #number of *
    	b=$(grep -r '~' $f | wc -l)   #number of ~
    	c=$(grep -v '#' $f | wc -l)   #total hits
    	d=$(($c-$a-$b))               #number of ' '
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
	    echo $min $a, $b, $d, $c, $average > info
	fi
	
    	cat info >> $2rscape.csv
done

paste -d ',' $1 $2rscape.csv $2cmbuildfinal.txt $2cmhmm.txt > $2.csv
rm -rf *.stk
