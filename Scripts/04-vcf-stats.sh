# Obtain population statistics

#####################################################################

# $1 is the input ncRNA dataset file

d=$(date +%y%m%d) # date identifier
id=$( echo $1 | cut -d '.' -f 1 | cut -d '-' -f 4 )
count=$(( $id + 1 ))
calc() { awk "BEGIN{print $*}"; }
IFS=$'\n'

#####################################################################

# Calculate population statistics

R --save << RSCRIPT

library(PopGenome)

# Make sure the FASTA folder only contains the VCF files

genome.class <- readData("FASTA", format="VCF",gffpath=FALSE)
genome.class <- neutrality.stats(genome.class)
dataset <- get.neutrality(genome.class)[[1]]
  
subset <- dataset[,c(1,2,4,5)]

write.csv(subset,"neutrality-stats.csv")     

RSCRIPT

#####################################################################

# Calculate average number of SNPs per ncRNA

cat neutrality-stats.csv | cut -d ',' -f 2,3,4,5 > stats
paste -d ',' $1 stats > vcfstats.csv
sed -i 's/n.segregating.sites/SNPs.1000/g' vcfstats.csv 

echo SNPsAverage.1000 > snps-average

for line in $( grep -v "Start" vcfstats.csv )
do
    chr=$( echo $line | cut -d ',' -f 3 )
    start=$( echo $line | cut -d ',' -f 4 )
    end=$( echo $line | cut -d ',' -f 5 )
    snp=$( echo $line | rev | cut -d ',' -f 3 | rev )
    length=$(( $end - $start ))
    average=$( calc $snp/$length )
    echo $average >> snps-average
done

paste -d ',' vcfstats.csv snps-average > $d-ncrna-dataset-$count.csv

#####################################################################

rm -rf stats
rm -rf snps-average
