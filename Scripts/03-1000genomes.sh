# Obtain VCF Files from 1000 Genomes Project and calculate transitions, transversions and minor allele frequency

#####################################################################

# $1 is the input ncRNA dataset file

IFS=$'\n'
d=$(date +%y%m%d) # date identifier
id=$( echo $1 | cut -d '.' -f 1 | cut -d '-' -f 4 )
count=$(( $id + 1 ))

#####################################################################

# Obtain VCF Files

grep -v 'Chromosome' $1 | cut -d ',' -f 3,4,5 | tr ',' ' ' | perl -lane '{print "$F[0]:$F[1]-$F[2]"}' > coordinates

IFS=$'\n'
var=1
rm -rf \FASTA
mkdir FASTA

for line in $(cat coordinates)
do
        name='RNA'$var
	echo $line | cut -d ':' -f 1 > text.file
	chr=$(cat text.file)
	chr="${chr%\"}"
	chr="${chr#\"}"
	
	# There	are two	different tags for the 1000genomes ftp files
	if [ $chr == 'X' ] 
	then
		tabix -f -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz $line > $name.vcf
	else
		tabix -f -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz $line > $name.vcf
	fi

	mv $name.vcf \FASTA
	var=$((var+1))
	rm -rf text.file
done

#####################################################################

# Calculate transitions, transversions and minor allele frequency

FILES=./FASTA/
data=$(ls $FILES | sort -V)

cd $FILES
cp ../$1 ./

A='A'
G='G'
T='T'
C='C'
calc() { awk "BEGIN{print $*}"; }
echo transition,transversion,TTratio,minMAF,maxMAF,aveMAF > freqsummary.csv

for f in $data
do
    vcftools --vcf $f --freq --out freq-afr &> /dev/null
    grep -v "CHROM" freq-afr.frq | cut -f 5,6 | perl -lane '{print "$F[0]:$F[1]"}' > snps.file
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
	total=$(calc $total+$maf)
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
    if (( $(echo "$ratio == 0.5" | bc -l) ))
    then
	echo $a,$b,0,0,$max,$average > line
    else
	echo $a,$b,$ratio,$min,$max,$average > line
    fi
    cat line >> freqsummary.csv

done

paste -d ',' $1 freqsummary.csv > $d-ncrna-dataset-$count.csv

#####################################################################

rm -rf freqsummary.csv
rm -rf line
rm -rf snps.file
