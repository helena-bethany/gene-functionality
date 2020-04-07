# Combining all code related to generating ncrna functionality traits

#################################################################

#exec 19>logfile
#BASH_XTRACEFD=19

#set -xv

IFS=$'\n'
date=$(date +%y%m%d)
calc() { awk "BEGIN{print $*}"; }
[ -d additional-output/ ] && rm -rf additional-output/* || mkdir additional-output

#################################################################

# Loading in all the correct user-made variables required.

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
read -p "RNA:RNA interaction database (fa): " interaction_database
echo
read -p "Human genome blastn database: " human_genome
echo
read -p "Processed snp151 database (bed): " snp_database
echo
echo "________________________________________________________________________________________"
echo
read -p "Enter folder location containing additional dependencies: " additional_dependencies
echo
# Need to specify boost lib folder for IntaRNA if you don't have pseudo permissions when installing.
if [ -f $additional_dependencies/IntaRNA ]
then
    read -p "Enter path for Boost C++ lib folder: " lib_directory
    echo
else
    :
fi

#################################################################

# Calculate CM-HMM bit scores from Rfam

##################################################################

# Run Rfam

rfam_model=$additional_folder/$rfam_name      # Define location of CM model

# Define location of required Infernal tools
if [ -f $additional_dependencies/cmscan ]
then
    cmscan_path=$additional_dependencies/cmscan
    cmfetch_path=$additional_dependencies/cmfetch
    cmpress_path=$additional_dependencies/cmpress
else
    cmscan_path=cmscan
    cmfetch_path=cmfetch
    cmpress_path=cmpress
fi

# Determine whether RNA variable ID starts at 1 (functional/positive control) or higher (negative control)
var=$( head -2 $initial_data | tail -1 | cut -d ',' -f 1 | tr -d RNA )
max=$( tail -1 $initial_data | cut -d ',' -f 1 | tr -d RNA )

echo ID,Rfam-BS > $date-cm-rfam
echo ID,Rfam-BS > $date-hmm-rfam

while [ $var -le $max ]
do
    grep -A 1 -m 1 "RNA$var" $initial_fasta > rna.fa
    # Run sequence on CM model with lower sensitivity to find top scoring model
    $cmscan_path --tblout cm_all.tblout --fmt 2 $rfam_model rna.fa &> /dev/null
    # Run sequence with only HMM and higher sensitivity
    $cmscan_path --hmmonly --hmmmax --tblout hmm_rfam.tblout --fmt 2 $rfam_model rna.fa &> /dev/null
    
    # Determine best scoring CM model from either CM or HMM-only search (if no hit to CM).
    if [ -f cm_all.tblout ]
    then
        model_name=$( grep "RNA$var" cm_all.tblout | head -1 | tr -s ' ' | cut -d ' ' -f 3 )
    elif [ -f hmm_rfam.tblout ]
    then
        model_name=$( grep "RNA$var" hmm_rfam.tblout | head -1 | tr -s ' ' | cut -d ' ' -f 3 )
    else
        model_name=
    fi

    if [ -z "$model_name" ]
    then
        :
    else
        # Extract the CM model of interest
        $cmfetch_path $rfam_model $model_name > $model_name.cm
        $cmpress_path $model_name.cm &> /dev/null
        # Run cmscan with high sensitivity parameters: runs on the inside and doesn't score to HMM at all.
        $cmscan_path --nohmmonly --max --tblout cm_rfam.tblout --fmt 2 $model_name.cm rna.fa &> /dev/null
    fi

    # Record CM bitscore
    if [ -f cm_rfam.tblout ]
    then
        cm=$( grep -m 1 "RNA$var" cm_rfam.tblout | tr -s ' ' | cut -d ' ' -f 17 )
        echo RNA$var,$cm >> $date-cm-rfam
    else
        :
    fi

    # Record HMM bitscore
    if [ -f hmm_rfam.tblout ]
    then
        hmm=$( grep -m 1 "RNA$var" hmm_rfam.tblout | tr -s ' ' | cut -d ' ' -f 17 )
        echo RNA$var,$hmm >> $date-hmm-rfam
    else
        :
    fi

    var=$(( $var+1 ))
    rm -rf *.cm.*
    rm -rf $model_name.cm

done

##################################################################

# Combine scores into one file and enter in theoretical zeros

echo Rfam-CM-HMM > bitscore-dataset.csv

for line in $( grep -v "Start" $initial_data )
do

    name=$( echo $line | cut -d ',' -f 1 )
    rfam_cm=$( grep -w "$name" $date-cm-rfam | cut -d ',' -f 2 )
    rfam_hmm=$( grep -w "$name" $date-hmm-rfam | cut -d ',' -f 2 )

    if [ -z "$rfam_cm" ] # If no match to CM model, then CM-HMM is recorded as zero
    then
        rfam_cm=0
        rfam_hmm=0
    elif [ -z "$rfam_hmm" ]  # If no match to HMM-only, then CM-HMM = CM
    then
        rfam_hmm=0
    else
        :
    fi

    rfam_diff=$( calc $rfam_cm-$rfam_hmm )

    echo $rfam_diff >> bitscore-dataset.csv

done

##################################################################

rm -rf cm_rfam.tblout
rm -rf hmm_rfam.tblout
rm -rf rna.fa
rm -rf cm_all.tblout
mv $date-cm-rfam additional-output/
mv $date-hmm-rfam additional-output/

##################################################################

# Obtain VCF Files from 1000 Genomes Project and calculate transitions, transversions and minor allele frequency

##################################################################

# Obtain VCF Files

# Alter count according to what number the RNA IDs start at
max=$( tail -1 $initial_data | cut -d ',' -f 1 | tr -d RNA )
var=$( head -2 $initial_data | tail -1 | cut -d ',' -f 1 | tr -d RNA )
count=1

# Reformat chromosome coordinates for tabix
grep -v 'Chromosome' $initial_data | cut -d ',' -f 3,4,5 | tr ',' ' ' | perl -lane '{print "$F[0]:$F[1]-$F[2]"}' | tr -d "chr" > coordinates

# Delete any previous VCF files
rm -rf FASTA.zip
mkdir FASTA

while [ $var -le $max ]
do
        # Grab the coordinates for each RNA as the counter increases
        line=$( sed -n "$count"p coordinates )
        name='RNA'$var
        chr=$( echo $line | cut -d ':' -f 1 )

        if [ -f $additional_dependencies/tabix ]
        then
            # The X chromosome is downloaded from a different folder to the autosomal chromosomes.
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
        count=$(( $count + 1 ))
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
    # Determine SNP frequencies for each ncRNA sequence
    if [ -f $additional_dependencies/vcftools ]
    then
        $additional_dependencies/vcftools --vcf FASTA/$f --freq --out freq-afr &> /dev/null
    else
        vcftools --vcf FASTA/$f --freq --out freq-afr &> /dev/null
    fi

    # Extract each SNP at frequency observed at that position into a parsable format
    grep -v "CHROM" freq-afr.frq | cut -f 5,6 | perl -lane '{print "$F[0]:$F[1]"}' > snps.file
    a=0
    b=0
    total=0
    count=$(grep -v "CHROM" freq-afr.frq | wc -l)
    max=0
    min=0.9

    for line in $(cat snps.file)
    do
        # Record frequencies for each SNP
        snpa=$(echo $line | cut -d ':' -f 1,3 | tr ":" " " | perl -lane '{print "$F[0]"}')
        snpb=$(echo $line | cut -d ':' -f 1,3 | tr ":" " " | perl -lane '{print "$F[1]"}')
        maf=$(echo $line | cut -d ":" -f 4)
        total=$(calc $total+$maf)
        # If the maximum is recorded as zero, set the new values automatically as the max
        if (( $(echo "$max == 0" | bc -l) ))
        then
            min=$maf
            max=$maf
        # If the new SNP analysed has a frequency greater than the recorded max, update the max
        elif (( $(echo "$maf > $max" | bc -l) ))
        then
            max=$maf
        # If the new SNP analysed has a frequency smaller than the recorded min, update the min
        elif (( $(echo "$maf < $min" | bc -l) ))
        then
            min=$maf
        else
            :
        fi
        # Counts the number of transitions and transversions
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

    # Prevents a division by zero if no SNPs were present
    if (( $(echo "$count == 0" | bc -l) ))
    then
        count=1
    else
        :
    fi

    # Calculate avergae allele frequency
    average=$(calc $total/$count)
    top=$(calc $a+1) 
    bottom=$(calc $a+$b+2)
    # Calculate ratio of transitions:transversion
    ratio=$(calc $top/$bottom)
    # Default ratio is 0.5 if no SNPs are present, should be changed to zero instead.
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
# Change title to make it more readable
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

#####################################################################

# Define function to convert MAF to Stockholm

maf_to_stockholm() {

# Extract all the multiz100way alignment information 
grep -v "#\|score" mafOut > mafOut.txt
# Split each file into individual alignments
awk -v RS= '{print > ("mafOut-" NR ".txt")}' mafOut.txt

IFS=$'\n'
max=$( ls mafOut-*.txt | wc -l )
count=1

while [ $count -le $max ]
do
    # Parse through each file, in order 
    for line in $( cat mafOut-$count.txt )
    do
        start=$( echo $line | tr -s ' ' | cut -d ' ' -f 1 )
        # If the line contains sequence information
        if (( $(echo "$start == s" | bc -l ) ))
        then
            # If this is the first file, then also extract sequence name
            if (( $(echo "$count == 1" | bc -l) ))
            then
                name=$( echo $line | tr -s ' ' | cut -d ' ' -f 2 | tr '.' '/' )
                seq=$( echo $line | tr -s ' ' | cut -d ' ' -f 7 )
                if [ -z "$seq" ]
                then
                    :
                else
                    # Ensures all the sequences start in the same position in the STK file
                    printf "%-40s\n" $name >> mafOut-0-seq
                    echo $seq >> mafOut-$count-seq
                fi      
            else
                # Only extract sequence information in subsequent files
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

# Format STK file
echo "# STOCKHOLM 1.0" > RNA.stk
echo >> RNA.stk
echo "#=GF ID $1" >> RNA.stk
echo >> RNA.stk

list=$( ls mafOut-*-seq | sort -V )
paste -d'\0' $list >> alignment

# Find the length of the sequence of the full sequence
human_len=$( grep "hg38" alignment | tr -s ' ' | cut -d ' ' -f 2 | wc -m )

for line in $( cat alignment )
do

    seq_len=$( echo $line | tr -s ' ' | cut -d ' ' -f 2 | wc -m )
    if [[ $seq_len -eq $human_len ]]   # Only include full sequences and removes partial alignments
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

var=$( head -2 $initial_data | tail -1 | cut -d ',' -f 1 | tr -d RNA )
rm -rf *.stk

# Reformats chromosome coordinates for mafFetch
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

        # Generate secondary structure consensus sequence
        if [ -f $additional_dependencies/RNAalifold ]
        then
            $additional_dependencies/RNAalifold -f S --aln-stk=$rna_id RNA.stk &> /dev/null
        else
            RNAalifold -f S --aln-stk=$rna_id RNA.stk &> /dev/null
        fi

        # Run R-scape wit high E-value threshold to better score negative control sequences
        if [ -f $additional_dependencies/R-scape ]
        then
            $additional_dependencies/R-scape -E 100 -s $rna_id.stk &> /dev/null
        else
            R-scape -E 100 -s $rna_id.stk &> /dev/null
        fi

        rm -rf *.sorted.cov
        mv *.cov rscapedata &> /dev/null
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

echo Min_covariance,Max_covariance,significant_bp,covariance_sum > rscape-dataset.csv

# Counter to process each rscape file individually and in order
total_rna=$( tail -1 $initial_data | tail -1 | cut -d ',' -f 1 | tr -d RNA )
count_file=$( head -2 $initial_data | tail -1 | cut -d ',' -f 1 | tr -d RNA )

while [ $count_file -le $total_rna ]
do
        if [ -f rscapedata/RNA$count_file.cov ]
        then
            f=./rscapedata/RNA$count_file.cov
            min=$(grep -r 'GTp' $f | cut -d '[' -f 2 | tr ']' ',' | cut -d ',' -f 1,2)   # Min and max covariance
            a=$(grep -r '*' $f | wc -l)   # number of signifcant basepairings
            grep -v "#" $f | cut -f 4 > score
            sum=0
            # Sum over all covariance
            for number in $( cat score )
            do
                sum=$(calc $sum+$number)   
            done

            echo $min,$a,$sum > info

        else
            # If no co-variance calculated/present
            echo 0,0,0,0 > info
        fi

        cat info >> rscape-dataset.csv
        count_file=$(( $count_file + 1 ))

done

#####################################################################

rm -rf *.ps
rm -rf *.stk
rm -rf coordinates
rm -rf overBed
rm -rf final.stk
rm -rf RNA*.fa
rm -rf info
rm -rf mafOut
rm -rf score
rm -rf mafOut.txt
zip -r rscapedata rscapedata &> /dev/null
rm -rf rscapedata/

#####################################################################

# Use IntaRNA to determine RNA:RNA interactions

#####################################################################

database=$additional_folder/$interaction_database

echo InteractionMIN,InteractionAVE > interaction-intermediate.csv
# Need to clear any previously set lib path, as otherwise the defined lib path will be appended onto the previous
unset LD_LIBRARY_PATH

for seq in $( grep ">" $initial_fasta )
do
    grep -A 1 $seq $initial_fasta > target.fa
    if [ -f $additional_dependencies/IntaRNA ]
    then
        LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$lib_directory $additional_dependencies/IntaRNA -q target.fa -t $database > intaRNA-results
    else
        IntaRNA -q target.fa -t $database > intaRNA-results
    fi

    # Grab all recorded interactions
    grep "interaction energy" intaRNA-results | cut -d '=' -f 2 | tr -d "kcal/mol" | tr -d ' ' > kcal_mol

    min=0
    count=$( grep -c "interaction energy" intaRNA-results )
    sum=0

    for number in $( cat kcal_mol )
    do
        # If the interaction energy is smaller than the recorded min, update the min
        if (( $( echo "$number < $min" | bc -l ) ))
        then
            min=$number
            sum=$( calc $sum+$number )
        else
            sum=$( calc $sum+$number )
        fi
    done

    # If no interactions were recorded, set average as 0 to prevent a division by zero error
    if (( $( echo "$count == 0" | bc -l) ))
    then
        ave=0
    else
        ave=$( calc $sum/$count )
    fi

    echo $min,$ave >> interaction-intermediate.csv

done

#####################################################################

rm -rf target.fa
rm -rf intaRNA-results
rm -rf kcal_mol

#####################################################################

# Calculate genomic copy number

#####################################################################

# Run blastn to identify number of matches within human genome

cp $additional_folder/$human_genome.n* .

echo Genome_copy_number,Genome_complete_match > copy-number.csv

if [ -f $additional_dependencies/blastn ]
then
    $additional_dependencies/blastn -query rna.fa -db $human_genome -outfmt "10 qaccver saccver pident" > output.csv
else
    blastn -query rna.fa -db $human_genome -outfmt "10 qaccver saccver pident" > output.csv
fi

count_file=$( head -2 $initial_data | tail -1 | cut -d ',' -f 1 | tr -d RNA )

while [ $count_file -le $total_rna ]
do
    # Process each RNA in the file (by ID) as the counter increases
    grep "RNA$count_file" output.csv > results
    total=$(cat results | wc -l)
    hundred=$(cut -d ',' -f 3 results | grep '100.000' | wc -l)
    echo $total,$hundred >> copy-number.csv
    count_file=$(( $count_file + 1 ))
done

#####################################################################

rm -rf results
rm -rf rna.fa
rm -rf output.csv
rm -rf sequences
rm -rf $human_genome.n*

#####################################################################

# dbsnp151 SNP count

# Format initial data so that it is in bed format
grep -v Chromosome $initial_data |  cut -f 3,4,5 -d "," | tr ',' '\t' > rnacentraloutput_FINAL.bed
sort -k1,1 -k2,2n rnacentraloutput_FINAL.bed | tr ' ' '\t' > rnacentraloutput_FINAL_sorted.bed

if [ -f $additional_dependencies/bedtools ]
then
    $additional_dependencies/bedtools intersect -c -a rnacentraloutput_FINAL_sorted.bed -b $additional_folder/$snp_database -wa -sorted > snp-intersection.bed
else
    bedtools intersect -c -a rnacentraloutput_FINAL_sorted.bed -b $additional_folder/$snp_database -wa -sorted > snp-intersection.bed
fi

echo Chromosome,Start,End,SNPs.UCSC,SNPsAverage.UCSC > snp-intersection-average.csv

for line in $( cat snp-intersection.bed )
do
    chr=$( echo $line | cut -f 1 | cut -c 4- )
    start=$( echo $line | cut -f 2 )
    end=$( echo $line | cut -f 3 )
    snp=$( echo $line | cut -f 4 )
    length=$(( $end - $start ))
    average=$( calc $snp/$length )
    echo $chr,$start,$end,$snp,$average >> snp-intersection-average.csv
done

#####################################################################

# Match up snps to ncRNA in R (as bedtools is unordered)

cat $initial_data > file1.csv
cat snp-intersection-average.csv > file2.csv

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

cat file3.csv > ncrna-dataset-snp.csv

####################################################################

rm -rf rnacentraloutput_FINAL_sorted1.bed
rm -rf rnacentraloutput_FINAL.bed
rm -rf rnacentraloutput_FINAL_sorted.bed
rm -rf file1.csv
rm -rf file2.csv
rm -rf file3.csv

mv snp-intersection.bed additional-output/
mv snp-intersection-average.csv additional-output/

#####################################################################

# Extracting expression for a region across 12 ENCODE datasets

#####################################################################

grep -v "Start" $initial_data | cut -d ',' -f 3,4,5 | tr ',' ' ' | perl -lane '{print "$F[0]:$F[1]-$F[2]"}' > locations

echo HepG2,hESC,K562,GM12878,SmoothMuscleCell,Hepatocyte,NeuralProgenitorCell,Myocyte,BipolarNeuron,Myotube,HematopoieticMultipotent,CardiacMuscle,MaxReadDepth > encode-rnaseq.csv

# Define location of samtools
if [ -f $additional_dependencies/samtools ]
then
        samtools_file=$additional_dependencies/samtools
else
        samtools_file=samtools
fi

for line in $(cat locations)
do
        # Obtain number of counts for each RNAseq dataset
        HepG2=$($samtools_file view -c $additional_folder/ENCFF067CVP.bam $line )
        hESC=$($samtools_file view -c $additional_folder/ENCFF089EWC.bam $line )
        K562=$($samtools_file view -c $additional_folder/ENCFF796BVP.bam $line )
        GM12878=$($samtools_file view -c $additional_folder/ENCFF893HSY.bam $line )
        smoothmuscle=$($samtools_file view -c $additional_folder/ENCFF369DYD.bam $line )
        hepatocyte=$($samtools_file view -c $additional_folder/ENCFF907AIK.bam $line )
        progenitor=$($samtools_file view -c $additional_folder/ENCFF766WSM.bam $line )
        myocyte=$($samtools_file view -c $additional_folder/ENCFF722EAR.bam $line )
        neuron=$($samtools_file view -c $additional_folder/ENCFF713UNS.bam $line )
        myotube=$($samtools_file view -c $additional_folder/ENCFF178TTA.bam $line )
        multipotent=$($samtools_file view -c $additional_folder/ENCFF065MVD.bam $line )
        cardiac=$($samtools_file view -c $additional_folder/ENCFF475WLJ.bam $line )

        # Obtain max read depth for each dataset, in the same order as above
        $samtools_file depth -r $line $additional_folder/ENCFF067CVP.bam | cut -f 3 | sort -V | tail -1 >> max
        $samtools_file depth -r $line $additional_folder/ENCFF089EWC.bam | cut -f 3 | sort -V | tail -1 >> max
        $samtools_file depth -r $line $additional_folder/ENCFF796BVP.bam | cut -f 3 | sort -V | tail -1 >> max
        $samtools_file depth -r $line $additional_folder/ENCFF893HSY.bam | cut -f 3 | sort -V | tail -1 >> max
        $samtools_file depth -r $line $additional_folder/ENCFF369DYD.bam | cut -f 3 | sort -V | tail -1 >> max
        $samtools_file depth -r $line $additional_folder/ENCFF907AIK.bam | cut -f 3 | sort -V | tail -1 >> max
        $samtools_file depth -r $line $additional_folder/ENCFF766WSM.bam | cut -f 3 | sort -V | tail -1 >> max
        $samtools_file depth -r $line $additional_folder/ENCFF722EAR.bam | cut -f 3 | sort -V | tail -1 >> max
        $samtools_file depth -r $line $additional_folder/ENCFF713UNS.bam | cut -f 3 | sort -V | tail -1 >> max
        $samtools_file depth -r $line $additional_folder/ENCFF178TTA.bam | cut -f 3 | sort -V | tail -1 >> max
        $samtools_file depth -r $line $additional_folder/ENCFF065MVD.bam | cut -f 3 | sort -V | tail -1 >> max
        $samtools_file depth -r $line $additional_folder/ENCFF475WLJ.bam | cut -f 3 | sort -V | tail -1 >> max

        max_depth=$( cat max | sort -V | tail -1 )
        
        # Obtain sequence length
        echo $line > line.txt
        start=$(cat line.txt | tr ":" "," | cut -d ',' -f 2 | tr '-' ' ' | perl -lane '{print "$F[0]"}')
        end=$(cat line.txt | tr ":" "," | cut -d ',' -f 2 | tr '-' ' ' | perl -lane '{print "$F[1]"}')
        length=$(calc $end-$start)

        # Calculate average number of counts
        HepG2_ave=$(calc $HepG2/$length)
        hESC_ave=$( calc $hESC/$length)
        K562_ave=$( calc $K562/$length)
        GM12878_ave=$( calc $GM12878/$length)
        smoothmuscle_ave=$( calc $smoothmuscle/$length)
        hepatocyte_ave=$( calc $hepatocyte/$length)
        progenitor_ave=$( calc $progenitor/$length)
        myocyte_ave=$( calc $myocyte/$length)
        neuron_ave=$( calc $neuron/$length)
        myotube_ave=$( calc $myotube/$length)
        multipotent_ave=$( calc $multipotent/$length)
        cardiac_ave=$( calc $cardiac/$length)

        echo $HepG2_ave,$hESC_ave,$K562_ave,$GM12878_ave,$smoothmuscle_ave,$hepatocyte_ave,$progenitor_ave,$myocyte_ave,$neuron_ave,$myotube_ave,$multipotent_ave,$cardiac_ave,$max_depth >> encode-rnaseq.csv
        rm -rf max
done

#####################################################################

rm -rf line.txt
rm -rf locations
rm -rf max

#####################################################################

# Calculate GC%

echo GC_percentage > GC.csv

for line in $( grep -v "Start" $initial_data )
do
    echo $line | cut -d ',' -f 6 > sequence
    G_content=$( grep -o "G\|g" sequence | wc -l )
    C_content=$( grep -o "C\|c" sequence | wc -l )
    A_content=$( grep -o "A\|a" sequence | wc -l )
    T_content=$( grep -o "T\|t" sequence | wc -l )
    GC_count=$(( $G_content + $C_content ))
    total=$(( $GC_count + $A_content + $T_content ))
    if (( $( echo "$GC_count == 0" | bc -l ) ))
    then
        GC=0
        percentage=0
    else
        GC=$( calc $GC_count/$total)
        percentage=$( calc $GC*100)
    fi
    echo $percentage >> GC.csv
done

rm -rf sequence

#####################################################################

# Now combine all the generated data together

# Update data parameter, in case script runs over multiple days
date=$(date +%y%m%d)
name=$( echo $initial_data | cut -d '-' -f 2,3 )

paste -d ',' ncrna-dataset-snp.csv bitscore-dataset.csv 1000g-freqsummary.csv 1000g-popstats.csv rscape-dataset.csv interaction-intermediate.csv copy-number.csv encode-rnaseq.csv GC.csv > $date-$name-final-dataset.csv

mv bitscore-dataset.csv additional-output/
mv 1000g-freqsummary.csv additional-output/
mv 1000g-popstats.csv additional-output/
mv rscape-dataset.csv additional-output/
mv interaction-intermediate.csv additional-output/
mv copy-number.csv additional-output/
mv encode-rnaseq.csv additional-output/
mv ncrna-dataset-snp.csv additional-output/
mv GC.csv additional-output/

#####################################################################

echo Finished calculating functionality traits.
echo
