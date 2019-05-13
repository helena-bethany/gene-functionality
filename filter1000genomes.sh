IFS=$'\n'
FILES=/home/helenacooper/bigDATA/1000genomes/RNA*
#FILES=/home/helenacooper/FASTA/
data=$(ls $FILES | sort -V)

A='A'
G='G'
T='T'
C='C'
calc() { awk "BEGIN{print $*}"; }
echo transition,transversion,TTratio,minMAF,maxMAF,aveMAF > freqsummary.csv
#min=0

for f in $data
do
    vcftools --vcf $f --freq --out freq-afr &> /dev/null
    grep -v "CHROM" freq-afr.frq | cut -f 5,6 | perl -lane '{print "$F[0]:$F[1]"}' > snps.file
    #grep -v "CHROM" freq-afr.frq | cut -f 6 | cut -d ":" -f 2 | tr ":" " " > maf
    a=0
    b=0
    total=0
    count=$(grep -v "CHROM" freq-afr.frq | wc -l)
    max=0
    min=0.9

    for line in $(cat snps.file)
    do
	snpa=$(echo $line | cut -d ':' -f 1,3 | tr ":" " " | perl -lane '{print "$F[0]"}')
	snpb=$(echo $line | cut -d ':' -f 1,3 | tr ":" " " | perl -lane '{print "$F[1]"}')
	maf=$(echo $line | cut -d ":" -f 4)
	#echo $maf
	total=$(calc $total+$maf)
	#echo $maf $min $max
	if (( $(echo "$max == 0" | bc -l) ))
	then
	    min=$maf
	    max=$maf
	elif (( $(echo "$maf > $max" | bc -l) ))
	then
	    max=$maf
	elif (( $(echo "$maf < $min" | bc -l) ))
	then
	    min=$maf
	else
	    :
	fi

	if [ "$snpa" == "$G" ]
	then
	    if [ "$snpb" == "$A" ]
	    then
		a=$((a+1))
	    elif [ "$snpb" == "$T" ]
	    then
		b=$((b+1))
	    elif [ "$snpb" == "$C" ]
	    then
		b=$((b+1))
	    else
		:
	    fi
	elif [ "$snpa" == "$A" ]
	then
	    if [ "$snpb" == "$G" ]
	    then
		a=$((a+1))
	    elif [ "$snpb" == "$T" ]
	    then
		b=$((b+1))
	    elif [ "$snpb" == "$C" ]
	    then
		b=$((b+1))
	    else
		:
	    fi
	elif [ "$snpa" == "$C" ]
	then
	    if [ "$snpb" == "$T" ]
	    then
		a=$((a+1))
	    elif [ "$snpb" == "$A" ]
	    then
		b=$((b+1))
	    elif [ "$snpb" == "$G" ]
	    then
		b=$((b+1))
	    else
		:
	    fi
	elif [ "$snpa" = "$T" ]
	then
	    if [ "$snpb" == "$C" ]
	    then
		a=$((a+1))
	    elif [ "$snpb" == "$A" ]
	    then
		b=$((b+1))
	    elif [ "$snpb" == "$G" ]
	    then
		b=$((b+1))
	    else
		:
	    fi
	else
	    :
	fi
    done

    if (( $(echo "$count == 0" | bc -l) ))
    then
	count=1
    else
	:
    fi
    
    average=$(calc $total/$count)
    top=$(calc $a+1)
    bottom=$(calc $a+$b+2)
    ratio=$(calc $top/$bottom)
    #echo $average
    if (( $(echo "$ratio == 0.5" | bc -l) ))
    then
	echo $a,$b,0,0,$max,$average > line
    else
	echo $a,$b,$ratio,$min,$max,$average > line
    fi
    cat line >> freqsummary.csv

done

paste -d ',' $1 freqsummary.csv > $2.csv
