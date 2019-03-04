#Please let this work

grep -v 'Chromosome' $1 | cut -d ',' -f 3,4,5 > locations
#grep -v '"' $2 | grep -v 'value' > phyloP
#grep -v '"' $3 | grep -v 'value' > phastCons

IFS=$'\n'
echo UCSCsnp,UCSCsnpaverage > $2.csv
calc() { awk "BEGIN{print $*}"; }
var=1

for line in $(cat locations)
do

    chrom=$(echo $line | cut -d ',' -f 1)
    #echo $chrom
    chromo=$(echo "'Chr$chrom'")
    start=$(echo $line | cut -d ',' -f 2)
    s=$(echo "'$start'")
    #echo $start
    end=$(echo $line | cut -d ',' -f 3)
    e=$(echo "'$end'")
    #echo $end
    snp=$(sudo mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -D hg38 -N -e "SELECT COUNT(chrom) FROM snp151 WHERE chrom = $chromo AND chromStart BETWEEN $s AND $e AND chromEnd BETWEEN $s AND $e;")
    #grep $chromo phyloP | grep $start | grep $end | cut -d ',' -f 8,10 > P
    #cat P
    #grep $chromo phastCons | grep $start | grep $end | cut -d ',' -f 8,10 > C
    #cat C
    seq=$(($end-$start))
    ave=$(calc $snp/$seq)
    #C=$(cat C)
    #P=$(cat P)
    echo $snp, $ave>> $2.csv
    echo RNA$var
    var=$(($var+1))

done

#paste -d ',' $1 cons.csv > $4csv

