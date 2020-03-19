# TEST Combining all code related to generating ncrna functionality traits

#################################################################

#set -xv
#sudo updatedb

IFS=$'\n'
date=$(date +%y%m%d)
calc() { awk "BEGIN{print $*}"; }
[ -d additional-output/ ] && rm -rf additional-output/* || mkdir additional-output

#################################################################

# Beginning needs to be loading in all the correct variables and a way of keeping track of them.

echo
echo Initial parameter set-up:
echo
read -p "Enter CSV file for initial dataset: " initial_data
echo
read -p "Enter FASTA file for initial dataset: " initial_fasta
echo
echo "________________________________________________________________________________________"
echo
read -p "Enter folder location containing additional input: " additional_folder
echo
echo "Type in the names of the following input within $additional_folder:"
echo
read -p "Rfam database (cm): " rfam_name
echo
read -p "Pfam database (hmm): " pfam_name
echo
read -p "Dfam database (hmm): " dfam_name
echo
read -p "RNA:RNA interaction database (suf): " interaction_database
echo
read -p "Human genome blastn database: " human_genome
echo
read -p "Processed snp151 database (bed): " snp_database
echo
echo "________________________________________________________________________________________"
echo
read -p "Enter folder location containing additional dependencies: " additional_dependencies
echo

#################################################################

# Calculate bit scores from Rfam, Pfam and Dfam

##################################################################

# Run Rfam

rfam_model=$additional_folder/$rfam_name

if [ -f $additional_dependencies/cmscan ]
then
    $additional_dependencies/cmscan --rfam --cut_ga --nohmmonly --tblout cm_rfam.tblout --fmt 2 $rfam_model $initial_fasta &> /dev/null
    $additional_dependencies/cmscan --hmmonly --tblout hmm_rfam.tblout --fmt 2 $rfam_model $initial_fasta &> /dev/null
else
    cmscan --rfam --cut_ga --nohmmonly --tblout cm_rfam.tblout --fmt 2 $rfam_model $initial_fasta &> /dev/null
    cmscan --hmmonly --tblout hmm_rfam.tblout --fmt 2 $rfam_model $initial_fasta &> /dev/null
fi

[ -f cm_rfam.tblout ] && cm_file=$( grep -v "#" cm_rfam.tblout ) || cm_file=' '
[ -f hmm_rfam.tblout ] && hmm_file=$( grep -v "#" hmm_rfam.tblout ) || hmm_file=' '

echo ID,Rfam-BS > $date-cm-rfam

for line in $cm_file
do

    rank=$(echo $line | tr -s ' ' | cut -d ' ' -f 1 )

    if (( $(echo "$rank == 1" | bc -l) ))
    then
        # Extract relevant information from cmscan output
        id=$(echo $line | tr -s ' ' | cut -d ' ' -f 4 )
        cm=$(echo $line | tr -s ' ' | cut -d ' ' -f 17 )

        echo $id,$cm >> $date-cm-rfam

    else
        :
    fi
    
done

echo ID,Rfam-BS > $date-hmm-rfam

for line in $hmm_file
do

    rank=$(echo $line | tr -s ' ' | cut -d ' ' -f 1 )

    if (( $(echo "$rank == 1" | bc -l) ))
    then
        #Extract relevant information from cmscan output
        id=$(echo $line | tr -s ' ' | cut -d ' ' -f 4 )
        hmm=$(echo $line | tr -s ' ' | cut -d ' ' -f 17 )

        echo $id,$hmm >> $date-hmm-rfam

    else
        :
    fi

done

###################################################################

# Run Pfam (note more than one hit per sequence)

pfam_model=$additional_folder/$pfam_name

if [ -f $additional_dependencies/transeq ]
then
    $additional_dependencies/transeq $initial_fasta $date-pfam.pep -frame=6 &> /dev/null
else
    transeq $initial_fasta $date-pfam.pep -frame=6 &> /dev/null
fi

if [ -f $additional_dependencies/hmmscan ]
then
    $additional_dependencies/hmmscan --tblout pfam.tbl -T 0 $pfam_model $date-pfam.pep &> /dev/null
else
    hmmscan --tblout pfam.tbl -T 0 $pfam_model $date-pfam.pep &> /dev/null
fi

echo ID,Pfam-BS > $date-pfam

[ -f pfam.tbl ] && pfam_file=$( grep -v "#" pfam.tbl ) || pfam_file=' '

for line in $pfam_file
do

    id=$( echo $line | tr -s ' ' | cut -d ' ' -f 3 | rev | cut -d '_' -f 2- | rev )
    score=$( echo $line | tr -s ' ' | cut -d ' ' -f 6 )
    echo $id,$score >> $date-pfam

done

##################################################################

# Run Dfam

dfam_model=$additional_folder/$dfam_name
dfam_script=$additional_dependencies/dfamscan.pl

perl $dfam_script -fastafile $initial_fasta -hmmfile $dfam_model -dfam_outfile dfam.out &> /dev/null

[ -f dfam.out ] && dfam_file=$( grep -v "#" dfam.out ) || dfam_file=' '

echo ID,Dfam-BS > $date-dfam

for line in $dfam_file
do

    id=$( echo $line | tr -s ' ' | cut -d ' ' -f 3 )
    score=$( echo $line | tr -s ' ' | cut -d ' ' -f 4 )
    echo $id,$score >> $date-dfam
done

##################################################################

# Combine scores into one file and enter in theoretical zeros

echo Rfam-CM,Rfam-HMM,Rfam-CM-HMM,Pfam-BS,Dfam-BS > bitscore-dataset.csv

for line in $( grep -v "Start" $initial_data )
do

    name=$( echo $line | cut -d ',' -f 1 )
    rfam_cm=$( grep -w "$name" $date-cm-rfam | cut -d ',' -f 2 )
    rfam_hmm=$( grep -w "$name" $date-hmm-rfam | cut -d ',' -f 2 )
    dfam_bs=$( grep -w "$name" $date-dfam | cut -d ',' -f 2 )
    pfam_bs=$( grep -w "$name" $date-pfam | cut -d ',' -f 2 | sort -V -r | head -1 )

    if [ -z "$rfam_cm" ]
    then
        rfam_cm=0
        rfam_diff=0
    else
        :
    fi

    if [ -z "$rfam_hmm" ]
    then
        rfam_hmm=0
    else
        :
    fi

    if [ -z "$dfam_bs" ]
    then
        dfam_bs=0
    else
        :
    fi

    if [ -z "$pfam_bs" ]
    then
        pfam_bs=0
    else
        :
    fi

    rfam_diff=$( calc $rfam_cm-$rfam_hmm )

    echo $rfam_cm,$rfam_hmm,$rfam_diff,$pfam_bs,$dfam_bs >> bitscore-dataset.csv

done

##################################################################

rm -rf cm_rfam.tblout
rm -rf hmm_rfam.tblout
mv $date-cm-rfam additional-output/
mv $date-hmm-rfam additional-output/
rm -rf pfam.tbl
mv $date-pfam additional-output/
rm -rf dfam.out
mv $date-dfam additional-output/
mv $date-pfam.pep additional-output/

##################################################################

# Obtain VCF Files from 1000 Genomes Project and calculate transitions, transversions and minor allele frequency

##################################################################

# Obtain VCF Files

grep -v 'Chromosome' $initial_data | cut -d ',' -f 3,4,5 | tr ',' ' ' | perl -lane '{print "$F[0]:$F[1]-$F[2]"}' | tr -d "chr" > coordinates

var=1
rm -rf FASTA.zip
mkdir FASTA

for line in $(cat coordinates)
do
        name='RNA'$var
        chr=$( echo $line | cut -d ':' -f 1 )

        if [ -f $additional_dependencies/tabix ]
        then
            if [ $chr == 'X' ]
            then
                $additional_dependencies/tabix -f -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz $line > $name.vcf
            else
                $additional_dependencies/tabix -f -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz $line > $name.vcf
            fi
        else
            if [ $chr == 'X' ]
            then
                tabix -f -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz $line > $name.vcf
            else
                tabix -f -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz $line > $name.vcf
            fi
        fi

        mv $name.vcf \FASTA
        var=$((var+1))
        rm -rf text.file

done

#####################################################################

# Calculate transitions, transversions and minor allele frequency

FILES=./FASTA/
data=$(ls $FILES | sort -V)

A='A'
G='G'
T='T'
C='C'

echo transition,transversion,TTratio,minMAF,aveMAF > 1000g-freqsummary.csv

for f in $data
do

    if [ -f $additional_dependencies/vcftools ]
    then
        $additional_dependencies/vcftools --vcf FASTA/$f --freq --out freq-afr &> /dev/null
    else
        vcftools --vcf FASTA/$f --freq --out freq-afr &> /dev/null
    fi

    grep -v "CHROM" freq-afr.frq | cut -f 5,6 | perl -lane '{print "$F[0]:$F[1]"}' > snps.file
    a=0
    b=0
    total=0
    count=$(grep -v "CHROM" freq-afr.frq | wc -l)
    max=0
    min=0.9

    for line in $(cat snps.file)
    do

        #echo $line
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
        echo $a,$b,0,0,$average > line
    else
        echo $a,$b,$ratio,$min,$average > line
    fi
    cat line >> 1000g-freqsummary.csv

done

#####################################################################

rm -rf line
rm -rf snps.file

#####################################################################

# Obtain population statistics

#####################################################################

# Calculate population statistics

run_popgenome() {

R --save << RSCRIPT
library(PopGenome)
# Make sure the FASTA folder only contains the VCF files
genome.class <- readData("FASTA", format="VCF",gffpath=FALSE)
genome.class <- neutrality.stats(genome.class)
dataset <- get.neutrality(genome.class)[[1]]

subset <- dataset[,c(1,2,4,5)]
write.csv(subset,"neutrality-stats.csv")
RSCRIPT

}

run_popgenome &> /dev/null

#####################################################################

# Calculate average number of SNPs per ncRNA

cat neutrality-stats.csv | cut -d ',' -f 2,3,4,5 > stats
sed -i 's/n.segregating.sites/SNPs.1000/g' stats
paste -d ',' $initial_data stats > vcfstats.csv

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

paste -d ',' stats snps-average > 1000g-popstats.csv

#####################################################################

rm -rf stats
mv snps-average additional-output/
mv neutrality-stats.csv additional-output/
mv vcfstats.csv additional-output/
rm -rf ALL.*.vcf.gz.tbi
rm -rf coordinates
rm -rf freq-afr.*

zip -r FASTA FASTA &> /dev/null
rm -rf FASTA/

#####################################################################

# Calculating ncRNA secondary structure metrics

#####################################################################

rm -rf rscapedata
mkdir rscapedata

var=1

#####################################################################

# Define biopython function to convert MAF to Stockholm

maf_to_stockholm() {

grep -v "#\|score" mafOut > mafOut.txt
awk -v RS= '{print > ("mafOut-" NR ".txt")}' mafOut.txt

IFS=$'\n'
max=$( ls mafOut-*.txt | wc -l )
count=1

while [ $count -le $max ]
do

    for line in $( cat mafOut-$count.txt )
    do
        start=$( echo $line | tr -s ' ' | cut -d ' ' -f 1 )
        if (( $(echo "$start == s" | bc -l ) ))
        then
            if (( $(echo "$count == 1" | bc -l) ))
            then
                name=$( echo $line | tr -s ' ' | cut -d ' ' -f 2 | tr '.' '/' )
                seq=$( echo $line | tr -s ' ' | cut -d ' ' -f 7 )
                if [ -z "$seq" ]
                then
                    :
                else
                    printf "%-40s\n" $name >> mafOut-0-seq
                    echo $seq >> mafOut-$count-seq
                fi      
            else
                seq=$( echo $line | tr -s ' ' | cut -d ' ' -f 7 )
                if [ -z "$seq" ]
                then
                    :
                else
                    echo $seq >> mafOut-$count-seq
                fi
            fi
        else
            :
        fi

    done
    count=$(( $count + 1 ))  

done

echo "# STOCKHOLM 1.0" > RNA.stk
echo >> RNA.stk
echo "#=GF ID $1" >> RNA.stk
echo >> RNA.stk

list=$( ls mafOut-*-seq | sort -V )
paste -d'\0' $list >> alignment

human_len=$( grep "hg38" alignment | tr -s ' ' | cut -d ' ' -f 2 | wc -m )

for line in $( cat alignment )
do

    seq_len=$( echo $line | tr -s ' ' | cut -d ' ' -f 2 | wc -m )
    if [[ $seq_len -eq $human_len ]]
    then
        echo $line >> RNA.stk
    else
        :
    fi
done

rm -rf mafOut-*.txt
rm -rf mafOut-*-seq
rm -rf alignment

}

#####################################################################

# Obtain multiple alignment files and calculate covariance

var=1
rm -rf *.stk

grep -v 'Chromosome' $initial_data | cut -d ',' -f 3,4,5 | tr ',' ' ' | perl -lane '{print "$F[0] $F[1] $F[2]"}' > coordinates

for line in  $(cat coordinates)
do
        echo $line > overBed

        if [ -f $additional_dependencies/mafFetch ]
        then
            $additional_dependencies/mafFetch hg38 multiz100way overBed mafOut
        else
            mafFetch hg38 multiz100way overBed mafOut
        fi

        rna_id='RNA'$var
        maf_to_stockholm $rna_id
        var=$((var+1))

        if [ -f $additional_dependencies/RNAalifold ]
        then
            $additional_dependencies/RNAalifold -f S --aln-stk=$rna_id RNA.stk &> /dev/null
        else
            RNAalifold -f S --aln-stk=$rna_id RNA.stk &> /dev/null
        fi

        if [ -f $additional_dependencies/R-scape ]
        then
            $additional_dependencies/R-scape -E 100 -s $rna_id.stk &> /dev/null
        else
            R-scape -E 100 -s $rna_id.stk &> /dev/null
        fi

        rm -rf *.sorted.cov
        mv *.cov rscapedata
        rm -rf overBed
        rm -rf blank.txt
        rm -rf *.pdf
        rm -rf *ss.ps
        rm -rf *.svg
        rm -rf *.surv
        rm -rf *.sum
        rm -rf *.sto
        rm -rf *.power

done

#####################################################################

# Process R-scape results

echo covMin10,covMax10,averageCov,#sigpair,#compatiable,#incompatiable,totalpairs > rscape-dataset.csv
FILES=./rscapedata/*.cov
data=$(ls $FILES | sort -V)

for f in $data
do
        min=$(grep -r 'GTp' $f | cut -d '[' -f 2 | tr ']' ',' | cut -d ',' -f 1,2)
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
            echo $min,$average,$a,$b,$d,$c > info
        fi

        cat info >> rscape-dataset.csv

done

#####################################################################

rm -rf *.stk
rm -rf coordinates
rm -rf overBedx
rm -rf final.stk
rm -rf RNAfinal.cm
rm -rf RNA*.fa
rm -rf info
rm -rf mafOut
rm -rf score
rm -rf mafOut.txt
zip -r rscapedata rscapedata &> /dev/null
rm -rf rscapedata/

#####################################################################

# Use RIsearch2 to determine RNA:RNA interactions

#####################################################################

# Use RIsearch2 to identify RNA:RNA interactions

database=$additional_folder/$interaction_database

if [ -f $additional_dependencies/risearch2.x ]
then
    $additional_dependencies/risearch2.x -q $initial_fasta -i $database
else
    risearch2.x -q $initial_fasta -i $database
fi

gunzip -f *.gz

FILES=./risearch_*.out
data=$(ls $FILES | sort -V)

echo RIsearchMIN,RIsearchAVE > risearch-intermediate.csv

for file in $data
do
    # Sums up every number in the 8th column
    count=$(awk '{s+=$8} END {print s}' $file)
    # Sorts the 8th column and takes the most negative number
    minimum=$(cut -f 8 $file | sort | head -1)
    # The number of lines in the file
    total=$(cat $file | wc -l)
    average=$(calc $count/$total)
    echo $minimum,$average >> risearch-intermediate.csv
done

#####################################################################

rm -rf data.fa

mkdir risearch
mv risearch_*.out risearch
zip -r risearch risearch &> /dev/null
rm -rf risearch/

#####################################################################

# Calculate number of times ncRNA appears within the human genome

#####################################################################

# Run blastn to identify number of matches within human genome

cp $additional_folder/$human_genome.n* .

echo Genome_copy_number,Genome_complete_match > copy-number.csv
grep -v 'Start' $initial_data | cut -d ',' -f 1,6 > sequences

for line in $(cat sequences)
do

    echo $line | tr ',' ' ' | perl -lane '{print ">$F[0]\n$F[1]"}' > rna.fa

    if [ -f $additional_dependencies/blastn ]
    then
        $additional_dependencies/blastn -query rna.fa -db $human_genome -outfmt "10 qaccver saccver pident" > output.csv
    else
        blastn -query rna.fa -db $human_genome -outfmt "10 qaccver saccver pident" > output.csv
    fi

    total=$(cat output.csv | wc -l)
    hundred=$(cut -d ',' -f 3 output.csv | grep '100.000' | wc -l)
    echo $total,$hundred >> copy-number.csv

done

#####################################################################

rm -rf rna.fa
rm -rf output.csv
rm -rf sequences
rm -rf $human_genome.n*

#####################################################################

# dbsnp151 SNP count

grep -v Chromosome $initial_data |  cut -f 3,4,5 -d "," | tr ',' '\t' > rnacentraloutput_FINAL.bed
sort -k1,1 -k2,2n rnacentraloutput_FINAL.bed | tr ' ' '\t' > rnacentraloutput_FINAL_sorted.bed
#awk '$1="chr"$1' rnacentraloutput_FINAL_sorted.bed | tr ' ' '\t' > rnacentraloutput_FINAL_sorted1.bed

if [ -f $additional_dependencies/bedtools ]
then
    $additional_dependencies/bedtools intersect -c -a rnacentraloutput_FINAL_sorted.bed -b $additional_folder/$snp_database -wa -sorted > $date-snp-intersection.bed
else
    bedtools intersect -c -a rnacentraloutput_FINAL_sorted.bed -b $additional_folder/$snp_database -wa -sorted > $date-snp-intersection.bed
fi

echo Chromosome,Start,End,SNPs.UCSC,SNPsAverage.UCSC > $date-snp-intersection-average.csv

for line in $( cat $date-snp-intersection.bed )
do
    chr=$( echo $line | cut -f 1 | cut -c 4- )
    start=$( echo $line | cut -f 2 )
    end=$( echo $line | cut -f 3 )
    snp=$( echo $line | cut -f 4 )
    length=$(( $end - $start ))
    average=$( calc $snp/$length )
    echo $chr,$start,$end,$snp,$average >> $date-snp-intersection-average.csv
done

#####################################################################

# Match up snps to ncRNA in R (as bedtools is unordered)

cat $initial_data > file1.csv
cat $date-snp-intersection-average.csv > file2.csv

reformat_snp() {

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

}

reformat_snp &> /dev/null

cat file3.csv > $date-ncrna-dataset-snp.csv

####################################################################

rm -rf rnacentraloutput_FINAL_sorted1.bed
rm -rf rnacentraloutput_FINAL.bed
rm -rf rnacentraloutput_FINAL_sorted.bed
rm -rf file1.csv
rm -rf file2.csv
rm -rf file3.csv

mv $date-snp-intersection.bed additional-output/

#####################################################################

# Extracting expression for a region across 12 ENCODE datasets

#####################################################################

# Cell line datasets

grep -v "Start" $initial_data | cut -d ',' -f 3,4,5 | tr ',' ' ' | perl -lane '{print "$F[0]:$F[1]-$F[2]"}' > locations

echo HepG2,hESC,K562,GM12878,AverageLINE,AveragebpLINE > $date-line.csv

if [ -f $additional_dependencies/samtools ]
then
        samtools_file=$additional_dependencies/samtools
else
        samtools_file=samtools
fi

for line in $(cat locations)
do
        HepG2=$($samtools_file view $additional_folder/ENCFF067CVP.bam $line | wc -l)
        hESC=$($samtools_file view $additional_folder/ENCFF089EWC.bam $line | wc -l)
        K562=$($samtools_file view $additional_folder/ENCFF796BVP.bam $line | wc -l)
        GM12878=$($samtools_file view $additional_folder/ENCFF893HSY.bam $line | wc -l)
        total=$(calc $HepG2+$hESC+$K562+$GM12878)
        average=$(calc $total/4)
        echo $line > line.txt
        start=$(cat line.txt | tr ":" "," | cut -d ',' -f 2 | tr '-' ' ' | perl -lane '{print "$F[0]"}')
        end=$(cat line.txt | tr ":" "," | cut -d ',' -f 2 | tr '-' ' ' | perl -lane '{print "$F[1]"}')
        length=$(calc $end-$start)
        bp=$(calc $average/$length)
        echo $HepG2,$hESC,$K562,$GM12878,$average,$bp > info
        cat info >> $date-line.csv
done

#####################################################################

# Differentiated cell type datasets

echo SmoothMuscleCell,Hepatocyte,NeuralProgenitorCell,Myocyte,BipolarNeuron,Myotube,HematopoieticMultipotent,CardiacMuscle,AverageTYPE,AveragebpTYPE > $date-type.csv

for line in $(cat locations)
do
        smoothmuscle=$($samtools_file view $additional_folder/ENCFF369DYD.bam $line | wc -l)
        hepatocyte=$($samtools_file view $additional_folder/ENCFF907AIK.bam $line | wc -l)
        progenitor=$($samtools_file view $additional_folder/ENCFF766WSM.bam $line | wc -l)
        myocyte=$($samtools_file view $additional_folder/ENCFF722EAR.bam $line | wc -l)
        neuron=$($samtools_file view $additional_folder/ENCFF713UNS.bam $line | wc -l)
        myotube=$($samtools_file view $additional_folder/ENCFF178TTA.bam $line | wc -l)
        multipotent=$($samtools_file view $additional_folder/ENCFF065MVD.bam $line | wc -l)
        cardiac=$($samtools_file view $additional_folder/ENCFF475WLJ.bam $line | wc -l)
        total=$(calc $smoothmuscle+$hepatocyte+$progenitor+$myocyte+$neuron+$myotube+$multipotent+$cardiac)
        average=$(calc $total/8)
        echo $line > line.txt
        start=$(cat line.txt | tr ":" "," | cut -d ',' -f 2 | tr '-' ' ' | perl -lane '{print "$F[0]"}')
        end=$(cat line.txt | tr ":" "," | cut -d ',' -f 2 | tr '-' ' ' | perl -lane '{print "$F[1]"}')
        length=$(calc $end-$start)
        bp=$(calc $average/$length)
        echo $smoothmuscle,$hepatocyte,$progenitor,$myocyte,$neuron,$myotube,$multipotent,$cardiac,$average,$bp > info
        cat info >> $date-type.csv
done

paste -d ',' $date-line.csv $date-type.csv > encode-rnaseq.csv

#####################################################################

rm -rf line.txt
rm -rf locations
rm -rf info
mv $date-line.csv additional-output/
mv $date-type.csv additional-output/

#####################################################################

# Now combine all the generated data together

date=$(date +%y%m%d)
name=$( echo $initial_data | cut -d '-' -f 2,3 )

paste -d ',' $date-ncrna-dataset-snp.csv bitscore-dataset.csv 1000g-freqsummary.csv 1000g-popstats.csv rscape-dataset.csv risearch-intermediate.csv copy-number.csv encode-rnaseq.csv > $date-$name-final-dataset.csv

mv bitscore-dataset.csv additional-output/
mv $date-snp-intersection-average.csv additional-output/
mv 1000g-freqsummary.csv additional-output/
mv 1000g-popstats.csv additional-output/
mv rscape-dataset.csv additional-output/
mv risearch-intermediate.csv additional-output/
mv copy-number.csv additional-output/
mv encode-rnaseq.csv additional-output/
mv $date-ncrna-dataset-snp.csv additional-output/

#####################################################################

echo Finished calculating functionality traits.
echo
