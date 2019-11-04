# Obtain dbSNP151 from UCSC and calculate number of SNPs in ncRNA

#####################################################################

# $1 is the input ncRNA dataset file
# $2 is either the gunzip or bed snp151 file

d=$(date +%y%m%d) # date identifier
id=$( echo $1 | cut -d '.' -f 1 | cut -d '-' -f 4 )
count=$(( $id + 1 ))
calc() { awk "BEGIN{print $*}"; }
IFS=$'\n'

#####################################################################

# Create local copy of dbSNP151 database

if [[ $2 == "snp151.txt.gz" ]]
then
    grep -v chr[A-Z]_* I think > snp151_no_alts
    cut -f 2,3,4 snp151_no_alts > snp151_trimmed.bed
    sort -k1,1 -k2,2n snp151_trimmed.bed > snp151_trimmed_sorted.bed
    dbsnp=snp151_trimmed_sorted.bed
else
    dbsnp=$2
fi
 
#####################################################################

# Create bedfile of ncRNA data

grep -v Chromosome $1 |  cut -f 3,4,5 -d "," | tr ',' '\t' > rnacentraloutput_FINAL.bed
sort -k1,1 -k2,2n rnacentraloutput_FINAL.bed > rnacentraloutput_FINAL_sorted.bed 
awk '$1="chr"$1' rnacentraloutput_FINAL_sorted.bed | tr ' ' '\t' > rnacentraloutput_FINAL_sorted1.bed

#####################################################################

# Calculate number of SNPs per ncRNA

bedtools intersect -c -a rnacentraloutput_FINAL_sorted1.bed -b $dbsnp -wa -sorted > $d-snp-intersection.bed

echo Chromosome,Start,End,SNPs.UCSC,SNPsAverage.UCSC > $d-snp-intersection-average.csv

for line in $( cat $d-snp-intersection.bed )
do
    chr=$( echo $line | cut -f 1 | cut -c 4- )
    start=$( echo $line | cut -f 2 )
    end=$( echo $line | cut -f 3 )
    snp=$( echo $line | cut -f 4 )
    length=$(( $end - $start ))
    average=$( calc $snp/$length )
    echo $chr,$start,$end,$snp,$average >> $d-snp-intersection-average.csv
done

#####################################################################

# Match up snps to ncRNA in R (as bedtools is unordered)

cat $1 > file1.csv
cat $d-snp-intersection-average.csv > file2.csv

R --save << RSCRIPT

df1 <- read.csv("file1.csv", stringsAsFactors=F)
df2 <- read.csv("file2.csv", stringsAsFactors=F)

for(i in 1:nrow(df2)){
	index <- grep(df2[i, 2], df1[,4])
	df1[index,'snp_num'] <- df2[i,4]
	df1[index,'snp_ave'] <- df2[i,5]
}

write.csv(df1, "file3.csv", quote=F, row.names=F)

RSCRIPT

cat file3.csv > $d-ncrna-dataset-$count.csv

#####################################################################

rm -rf file1.csv
rm -rf file2.csv
rm -rf file3.csv

