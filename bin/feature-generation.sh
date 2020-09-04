#!/bin/bash
#
# Script Name: feature-generation.sh
#
# Author: Helena Cooper
# Last edited: 04/09/2020
#
# Description: This script calculates all chosen features of functions for protein-coding exons, lncRNA exons, short ncRNAs 
#              and negative control sequences.
#
# Input: $1 is the file identifier for the dataset and fasta file. Eg: 200702-functional-ncrna
#        $2 is the location of the folder with the local databases/datasets or version specific executables
# 
# Any additional required files, directories or dependencies will be requested when the script is run and require manual
# input.
#
# Log files: tabix.log records any potential VCF files that weren't downloaded correctly.
#

###########################################################################################################################

######## General setup

#set -xv
IFS=$'\n'
date=$(date +%y%m%d)

######## Function for calculation with float numbers

calc() { awk "BEGIN{print $*}"; }

######## Create folder for additional files generated

[ -d additional-output/ ] && rm -rf additional-output/* || mkdir additional-output
echo

######## Delete previous log files

rm -rf tabix.log
rm -rf errors.log

##########################################################################################################################

# Assign all required input files to individual variables

######## CSV and FASTA for initial data

initial_data=$1-dataset.csv     
initial_fasta=$1-seq.fa         
if [ ! -f $initial_data ] || [ ! -f $initial_fasta ]
then
    echo "Initial data files do not exist."
    exit 1
else
    :
fi

######## Folder with local databases and version specific executables

additional_folder=$2          
if [ ! -d $additional_folder/ ]
then
    echo "Folder with required local databases does not exist."
    exit 1
else
    :
fi

######## Interaction database for IntaRNA

if [ -f $additional_folder/curated-interaction-database.fa ] 
then
    interaction_database=$additional_folder/curated-interaction-database.fa
else
    until [ -f $interaction_database ] ; do read -p "Please enter custom name of RNA:RNA interaction database (fa): " i ; interaction_database=$additional_folder/$i ; echo ; done
fi

######## Local database for blastn

blast_total=$( ls $additional_folder/human_genome* | wc -l )
if [[ $blast_total -eq 6 ]]
then
    human_genome=$additional_folder/human_genome
else
    until [ -f $human_genome.nhr ] ; do read -p "Please enter custom name of human genome blastn database: " h ; human_genome=$additional_folder/$h ; echo ; done
fi

######## Local bedfile of Dfam hits for bedtools

if [ -f $additional_folder/dfam-hg38-sorted.bed ]        
then
    dfam_hits=$additional_folder/dfam-hg38-sorted.bed
else
    until [ -f $dfam_hits ] ; do read -p "Please enter custom name of Dfam non-redundant hits (bed): " d ; dfam_hits=$additional_folder/$d ; echo ; done
fi
######## Conversion file for changing coordinates between hg19 and hg38 for liftOver

if [ -f $additional_folder/hg38ToHg19.over.chain ]
then
    chain_file=$additional_folder/hg38ToHg19.over.chain
else
    until [ -f $chain_file ] ; do read -p "Please enter custom name of coordinates file for converting between hg38 and Hg19 genome versions (chain): " chain ; chain_file=$additional_file/$chain ; echo ; done
fi

######## Local bigwig file of UCSC phyloP scores for BigWigSummary

if [ -f $additional_folder/hg38.phyloP100way.bw ] 
then
    phylo_bw=$additional_folder/hg38.phyloP100way.bw
else
    until [ -f $phylo_bw ] ; do read -p "Please enter custom name of hg38 conservation scoring by phyloP (bigWig): " phylo ; phylo_bw=$additional_folder/$phylo ; echo ; done
fi

######## Local bigwig file of UCSC phastCons scores for BigWigSummary

if [ -f $additional_folder/hg38.phastCons100way.bw ] 
then
    phast_bw=$additional_folder/hg38.phastCons100way.bw 
else
    until [ -f $phast_bw ] ; do read -p "Please enter custom name of hg38 conservation scoring by phastCons (bigWig): " phast ; phast_bw=$additional_folder/$phast ; echo ; done
fi

######## ENCODE RNA-seq BAM files for samtools

bam_total=$( ls $additional_folder/*.bam | wc -l )
if [[ $bam_total -eq 12 ]]
then
    :
else
    echo "Please confirm you have all 12 ENCODE Total RNA-Seq bam files downloaded."
    echo
    exit 1 
fi

###########################################################################################################################

# Assign all required executables a variable

######## liftOver

if find $additional_folder/ -executable -type f | grep -q liftOver
then
    liftOver_exe=$additional_folder/liftOver
elif command -v liftOver &> /dev/null
then
    liftOver_exe=liftOver
else
    echo Please check that liftOver is installed.
    exit 1
fi

######## tabix

if find $additional_folder/ -executable -type f | grep -q tabix
then
    tabix_exe=$additional_folder/tabix
elif command -v tabix &> /dev/null
then
    tabix_exe=tabix
else
    echo Please check that tabix is installed.
    exit 1
fi

######## vcftools

if find $additional_folder/ -executable -type f | grep -q vcftools
then
    vcftools_exe=$additional_folder/vcftools
elif command -v vcftools &> /dev/null
then
    vcftools_exe=vcftools
else
    echo Please check that vcftools is installed.
    exit 1
fi

######## mafFetch

if find $additional_folder/ -executable -type f | grep -q mafFetch
then
    mafFetch_exe=$additional_folder/mafFetch
elif command -v mafFetch &> /dev/null
then
    mafFetch_exe=mafFetch
else
    echo Please check that mafFetch is installed.
    exit 1
fi

######## bigWigSummary

if find $additional_folder/ -executable -type f | grep -q bigWigSummary
then
    bigWigSummary_exe=$additional_folder/bigWigSummary
elif command -v bigWigSummary &> /dev/null
then
    bigWigSummary_exe=bigWigSummary
else
    echo Please check that bigWigSummary is installed.
    exit 1
fi

######## RNAcode

if find $additional_folder/ -executable -type f | grep -q RNAcode
then
    RNAcode_exe=$additional_folder/RNAcode
elif command -v RNAcode &> /dev/null
then
    RNAcode_exe=RNAcode
else
    echo Please check that RNAcode is installed.
    exit 1
fi

######## RNAalifold

if find $additional_folder/ -executable -type f | grep -q RNAalifold
then
    RNAalifold_exe=$additional_folder/RNAalifold
elif command -v RNAalifold &> /dev/null
then
    RNAalifold_exe=RNAalifold
else
    echo Please check that RNAalifold is installed.
    exit 1
fi

######## R-scape

if find $additional_folder/ -executable -type f | grep -q 'R-scape'
then
    rscape_exe=$additional_folder/R-scape
elif command -v R-scape &> /dev/null
then
    rscape_exe=R-scape
else
    echo Please check that R-scape is installed.
    exit 1
fi

######## IntaRNA

if find $additional_folder/ -executable -type f | grep -q IntaRNA
then
    IntaRNA_exe=$additional_folder/IntaRNA
elif command -v IntaRNA &> /dev/null
then
    IntaRNA_exe=IntaRNA
else
    echo Please check that IntaRNA is installed.
    exit 1
fi

intaRNA_error=$( $IntaRNA_exe -h 2>&1 )  # Need to find out if Boost library needs to be defined separately for IntaRNA to run

if [ ! -z "$error" ] && [[ "$error" != "IntaRNA predictors RNA-RNA interactions"* ]]
then
    read -p "Enter path for Boost C++ lib folder for IntaRNA: " lib_directory
    echo
    if [ -d $lib_directory/ ] ; then : ; else echo Please check that the Boost C++ lib directory is correct. ; exit 1 ; fi
else
    :
fi

######## blastn

if find $additional_folder/ -executable -type f | grep -q blastn
then
    blastn_exe=$additional_folder/blastn
elif command -v blastn &> /dev/null
then
    blastn_exe=blastn
else
    echo Please check that blastn is installed.
    exit 1
fi

######## bedtools

if find $additional_folder/ -executable -type f | grep -q bedtools
then
    bedtools_exe=$additional_folder/bedtools
elif command -v bedtools &> /dev/null
then
    bedtools_exe=bedtools
else
    echo Please check that bedtools is installed.
    exit 1
fi

######## samtools

if find $additional_folder/ -executable -type f | grep -q samtools
then
    samtools_exe=$additional_folder/samtools
elif command -v samtools &> /dev/null
then
    samtools_exe=samtools
else
    echo Please check that samtools is installed.
    exit 1
fi

######## RNAfold

if find $additional_folder/ -executable -type f | grep -q RNAfold
then
    RNAfold_exe=$additional_folder/RNAfold
elif command -v RNAfold &> /dev/null
then
    RNAfold_exe=RNAfold
else
    echo Please check that RNAfold is installed.
    exit 1
fi

######## RNAplfold

if command -v RNAplfold &> /dev/null
then
    :
else
    echo Please check that RNAplfold is installed and exported to $PATH.
    exit 1
fi

######## access_py.py

if [ -f $additional_folder/access_py.py ]
then
    :
else
    echo Please check that access_py.py is available in $additional_folder.
fi

######## CPC2

if [ -f $additional_folder/CPC2.py ]
then
    cpc2_directory=$additional_folder/CPC2.py
else
    read -p "Enter path to installed CPC2 bin: " cpc2_directory
    echo
    if [ -d $cpc2_directory/ ] ; then : ; else echo Please check that CPC2 is installed. ; exit 1 ; fi
fi

############################################################################################################################

# Obtain VCF Files from 1000 Genomes Project.

############################################################################################################################

######## Obtain VCF Files
rm -rf coordinates

######## Alter count according to what number the RNA IDs start at
max=$( tail -1 $initial_data | cut -d ',' -f 1 | tr -d RNA )
var=$( head -2 $initial_data | tail -1 | cut -d ',' -f 1 | tr -d RNA )
count=1

######## Reformat chromosome coordinates for to obtain hg19 coordinates
grep -v 'Chromosome' $initial_data | cut -d ',' -f 3,4,5 | tr ',' '\t' > coordinates
grep -v 'Chromosome' $initial_data | cut -d ',' -f 1 > id
paste -d '\t' coordinates id > input.bed

######## Convert coordinates to hg19 genome version
$liftOver_exe input.bed $chain_file output.bed unlifted.bed &> /dev/null

######## If available, unzip previously downloaded VCF files
if [ -f FASTA.zip ]
then
    unzip FASTA.zip &> /dev/null
else
    mkdir FASTA
    
    ######## Downloaded VCF files according to hg19 chromosome coordinates
    while [ $var -le $max ]
    do
        ######## Grab the coordinates for each RNA as the counter increases
        line=$( grep -w "RNA$var" output.bed | cut -f 1,2,3 | tr '\t' ' ' | tr -d "chr" | perl -lane '{print "$F[0]:$F[1]-$F[2]"}' )
        
        ######## If previous coordinate was associated with an annoated chromosome, use hg38 coordinates
        if [ -z $line ] || [[ $line == *"Un_"* ]]
        then
            line=$( grep -w "RNA$var" $initial_data | cut -d ',' -f 3,4,5 | tr ',' ' ' | tr -d "chr" | perl -lane '{print "$F[0]:$F[1]-$F[2]"}' )
                
        ######## Reformatting scaffolds of annotated chromosomes
        elif [[ $line == *"_"* ]]
        then
            a=$( echo $line | cut -d "_" -f 1 )
            b=$( echo $line | cut -d ":" -f 2 )
            line=$( echo $a:$b )
        else
            :
        fi
        
        name='RNA'$var
        chr=$( echo $line | cut -d ':' -f 1 )

        ######## The X chromosome is downloaded from a different folder to the autosomal chromosomes.
        if [ $chr == 'X' ]
        then
            until $tabix_exe -f -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz $line > $name.vcf 2>>tabix.log ; do echo Downloading $name.vcf >> tabix.log ; echo >> tabix.log ; done
            
        ######## Tabix is very sensitive to crashing, so repeat the command until a valid VCF file has been downloaded.
        else
            until $tabix_exe -f -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz $line > $name.vcf 2>>tabix.log ; do echo Downloading $name.vcf >> tabix.log ; echo >> tabix.log ; done
        fi

        mv $name.vcf \FASTA
        var=$((var+1))
        rm -rf text.file
        count=$(( $count + 1 ))
    done
fi

######## Remove excess files
rm -rf output.bed
rm -rf input.bed
rm -rf unlifted.bed

############################################################################################################################

# Calculate transitions:transversions ratio and minor allele frequency.

############################################################################################################################

######## Specify location of VCF files
FILES=./FASTA/
data=$(ls $FILES | sort -V)

######## Define nucleotides as variables
A='A'
G='G'
T='T'
C='C'

######## Loop through VCF files to calculate features
echo TTratio,aveMAF > 1000g-freqsummary.csv

for f in $data
do
    ######## Determine SNP frequencies for each ncRNA sequence
    $vcftools_exe --vcf FASTA/$f --freq --out freq-afr &> /dev/null

    ######## Extract each SNP at frequency observed at that position into a parsable format
    grep -v "CHROM" freq-afr.frq | cut -f 5,6 | perl -lane '{print "$F[0]:$F[1]"}' > snps.file
    a=0  # No. of transitions
    b=0  # No. of transversion
    total=0   # No. of SNPs observed
    count=$(grep -v "CHROM" freq-afr.frq | wc -l)   # Maximum data available
    max=0    # Max MAF observed
    min=0.9  # Min MAF observed

    for line in $(cat snps.file)
    do
        echo $line > test
        if grep -q ">" test   # Stop loop crashing if VCF file is empty
        then
            rm -rf test
        else
            ######## Record frequencies for each SNP
            snpa=$(echo $line | cut -d ':' -f 1,3 | tr ":" " " | perl -lane '{print "$F[0]"}')
            snpb=$(echo $line | cut -d ':' -f 1,3 | tr ":" " " | perl -lane '{print "$F[1]"}')
            maf=$(echo $line | cut -d ":" -f 4)
            total=$(calc $total+$maf)
            ######## If the maximum is recorded as zero, set the new values automatically as the max
            if (( $(echo "$max == 0" | bc -l) ))
            then
                min=$maf
                max=$maf
            ######## If the new SNP analysed has a frequency greater than the recorded max, update the max
            elif (( $(echo "$maf > $max" | bc -l) ))
            then
                max=$maf
            ######## If the new SNP analysed has a frequency smaller than the recorded min, update the min
            elif (( $(echo "$maf < $min" | bc -l) ))
            then
                min=$maf
            else
                :
            fi
            ######## Counts the number of transitions and transversions
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
        fi
    done

    ######## Prevents a division by zero if no SNPs were present, include NA instead.
    if (( $(echo "$count == 0" | bc -l) ))
    then
        average=NA
    else
        average=$( calc $a+1 )
    fi

    ######## Calculate transitions:transversion ratio
    top=$(calc $a+1)
    bottom=$(calc $a+$b+2)
    if [[ $a -eq 0 ]] && [[ $b -eq 0 ]] ; then ratio=NA ; else ratio=$( calc $top/$bottom ) ; fi 
    echo $ratio,$average >> 1000g-freqsummary.csv

done

rm -rf line
rm -rf snps.file

############################################################################################################################

# Obtain population statistics from VCF files

############################################################################################################################

######## Function for calculating population statistics in R

run_popgenome() {

R --save << RSCRIPT

######## If PopGenome not install, install package. Otherwise load required library.
packages <- "PopGenome"

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

######## Make sure the FASTA folder only contains the VCF files
genome.class <- readData("FASTA", format="VCF",gffpath=FALSE)
genome.class <- neutrality.stats(genome.class)
dataset <- get.neutrality(genome.class)[[1]]

subset <- dataset[,c(1,2,4,5)]
write.csv(subset,"neutrality-stats.csv")
RSCRIPT

}

############################################################################################################################

######## Calculate Tajima's D and Fu and Li's D using R
run_popgenome &> /dev/null

######## Change title to make it more readable
cat neutrality-stats.csv | cut -d ',' -f 2,3,4,5 > stats
sed -i 's/n.segregating.sites/SNPsCount_1000G/g' stats
sed -i 's/"Tajima.D"/Tajima_D/g' stats
sed -i 's/"Fu.Li.D"/Fu_Li_D/g' stats 
paste -d ',' $initial_data stats > vcfstats.csv

echo SNPsDensity_1000G > snps-average

######## Calculate average number of SNPs per ncRNA
for line in $( grep -v "Start" vcfstats.csv )
do
    chr=$( echo $line | cut -d ',' -f 3 )
    start=$( echo $line | cut -d ',' -f 4 )
    end=$( echo $line | cut -d ',' -f 5 )
    snp=$( echo $line | rev | cut -d ',' -f 3 | rev )
    length=$(( $end - $start ))
    average=$( calc $snp/$length )
    [[ $snp -eq 0 ]] && average=NA
    echo $average >> snps-average
done

paste -d ',' stats snps-average > 1000g-popstats.csv

######## Delete and move excess files
rm -rf stats
mv snps-average additional-output/
mv neutrality-stats.csv additional-output/
mv vcfstats.csv additional-output/
rm -rf ALL.*.vcf.gz.tbi
rm -rf coordinates
rm -rf freq-afr.*

######## Zip VCF files for further analyses
zip -r FASTA FASTA &> /dev/null
rm -rf FASTA/

############################################################################################################################

# Calculating ncRNA secondary structure and protein coding potential features

############################################################################################################################

# Define function to convert MAF to Stockholm

maf_to_stockholm() {

######## Extract all the multiz100way alignment information 
grep -v "#\|score" mafOut > mafOut.txt

######## Split each file into individual alignments
awk -v RS= '{print > ("mafOut-" NR ".txt")}' mafOut.txt

IFS=$'\n'
max=$( ls mafOut-*.txt | wc -l )
count=1

while [ $count -le $max ]
do
    ######## Parse through each file, in order 
    for line in $( cat mafOut-$count.txt )
    do
        start=$( echo $line | tr -s ' ' | cut -d ' ' -f 1 )
        ######## If the line contains sequence information
        if (( $(echo "$start == s" | bc -l ) ))
        then
            ######## If this is the first file, then also extract sequence name
            if (( $(echo "$count == 1" | bc -l) ))
            then
                name=$( echo $line | tr -s ' ' | cut -d ' ' -f 2 | tr '.' '/' )
                seq=$( echo $line | tr -s ' ' | cut -d ' ' -f 7 )
                if [ -z "$seq" ]
                then
                    :
                else
                    ######## Ensures all the sequences start in the same position in the STK file
                    printf "%-40s\n" $name >> mafOut-0-seq
                    echo $seq >> mafOut-$count-seq
                fi      
            else
                ######## Only extract sequence information in subsequent files
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

######## Format STK file
echo "# STOCKHOLM 1.0" > RNA.stk
echo >> RNA.stk
echo "#=GF ID $1" >> RNA.stk
echo >> RNA.stk

list=$( ls mafOut-*-seq | sort -V )
paste -d'\0' $list >> alignment

######## Find the length of the sequence of the full sequence
human_len=$( grep "hg38" alignment | tr -s ' ' | cut -d ' ' -f 2 | wc -m )

for line in $( cat alignment )
do

    seq_len=$( echo $line | tr -s ' ' | cut -d ' ' -f 2 | wc -m )
    if [[ $seq_len -eq $human_len ]]     # Only include full sequences and removes partial alignments
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

############################################################################################################################

######## Create folder for generated r-scape output
rm -rf rscapedata && mkdir rscapedata

######## If multiz100way files already downloaded, use them instead of re-downloading 
if [ -f maf.zip ] ; then unzip maf.zip &> /dev/null ; else mkdir maf ; fi

######## Reformats chromosome coordinates for mafFetch
rm -rf coordinates
grep -v 'Chromosome' $initial_data | cut -d ',' -f 3,4,5 | tr ',' ' ' | perl -lane '{print "$F[0] $F[1] $F[2]"}' > coordinates

######## Obtain multiple alignment files and calculate covariance
var=$( head -2 $initial_data | tail -1 | cut -d ',' -f 1 | tr -d RNA )
rm -rf *.stk

echo RNAcode_score,RNAalifold_score > rnacode.csv
echo MeanPhyloP,MaxPhyloP,MeanPhastCons,MaxPhastCons > conservation.csv

for line in  $(cat coordinates)
do
        ######## Format data for phastCons and phyloP
        echo $line > overBed
        chr=$( echo $line | cut -d ' ' -f 1 )
        start=$( echo $line | cut -d ' ' -f 2 )
        end=$( echo $line | cut -d ' ' -f 3 )
        
        ######## Obtain phastCons and phyloP values for each set of chromosome coordinates
        mean_pp=$( $bigWigSummary_exe -type=mean $phylo_bw $chr $start $end 1 2>&1 )
        max_pp=$( $bigWigSummary_exe -type=max $phylo_bw $chr $start $end 1 2>&1 )
        mean_pc=$( $bigWigSummary_exe -type=mean $phast_bw $chr $start $end 1 2>&1 )
        max_pc=$( $bigWigSummary_exe -type=max $phast_bw $chr $start $end 1 2>&1 )

        ######## Convert any missing data to NAs
        test_mnp=$( echo $mean_pp | wc -w )
        if [[ "$test_mnp" -eq "1" ]] ; then : ; else mean_pp=NA ; fi
 
        test_mxp=$( echo $max_pp | wc -w )
        if [[ "$test_mxp" -eq "1" ]] ; then : ; else max_pp=NA ; fi

        test_mnc=$( echo $mean_pc | wc -w )
        if [[ "$test_mnc" -eq "1" ]] ; then : ; else mean_pc=NA ; fi

        test_mxc=$( echo $max_pc | wc -w )
        if [[ "$test_mxc" -eq "1" ]] ; then : ; else max_pc=NA ; fi

        echo $mean_pp,$max_pp,$mean_pc,$max_pc >> conservation.csv
        
        ######## Obtain MSA files from multiz100way, unless they've already been downloaded
        rna_id='RNA'$var
        if [ -f maf/$rna_id.maf ]
        then
            cp maf/$rna_id.maf mafOut
        else            
            until $mafFetch_exe hg38 multiz100way overBed mafOut 2>>errors.log ; do sleep 4 ; done
            cp mafOut $rna_id.maf && mv $rna_id.maf maf/
        fi

        file_len=$( cat mafOut | wc -l )
        
        ######## If MSA contains data, convert is to STK format
        if [ "$file_len" != "1" ]
        then
            maf_to_stockholm $rna_id

            ######## Run RNAcode using original mafOut file
            $RNAcode_exe mafOut -o rnacode_output 2>>errors.log &> /dev/null
            
            ######## Takes the largest score, but records zero if no significant hits found.
            rnacode=$( grep -v "No significant coding regions found.\|#\|=" rnacode_output | grep . | tr -s ' ' | cut -d ' ' -f 10 | grep . | sort -V | tail -1 )
            [ -z $rnacode ] && rnacode=NA

            ######## Generate secondary structure consensus sequence and associated MFE value
            timeout 60m $RNAalifold_exe -q -f S --noPS --aln-stk=$rna_id RNA.stk >/dev/null 2>>errors.log
            rna_score=$( $RNAalifold_exe -q -f S --noPS RNA.stk 2>>errors.log | tail -1 | cut -d ' ' -f 2- | tr -d "(" | cut -d "=" -f 1 )

            [ -z "$rna_score" ] && rna_score=NA
            echo $rnacode,$rna_score | tr -d ' ' >> rnacode.csv

            exit_status=$?

            ######## Run R-scape only if previous analyses didn't time out
            if [ "$exit_status" -ne "124" ]
            then
                $rscape_exe -E 100 -s $rna_id.stk >/dev/null 2>>errors.log
            else
                touch $rna_id.cov
            fi
        else
            ######## Generate blanks if no MAF available
            touch $rna_id.cov
            echo NA,NA >> rnacode.csv
        fi
        
        var=$(($var+1))
        rm -rf *.sorted.cov
        mv *.cov rscapedata &> /dev/null
        rm -rf overBed
        rm -rf blank.txt
        rm -rf *.pdf
        rm -rf *ss.ps
        rm -rf *.svg
        rm -rf *.surv
        rm -rf *.sto
        rm -rf *.power
        rm -rf rnacode_output

done

zip -r maf maf &> /dev/null
rm -rf maf/

############################################################################################################################

# Process R-scape results to obtain covariance features

############################################################################################################################

echo Max_covariance,Sum_covariance > rscape-dataset.csv

######## Counter to process each rscape file individually and in order
total_rna=$( tail -1 $initial_data | tail -1 | cut -d ',' -f 1 | tr -d RNA )
count_file=$( head -2 $initial_data | tail -1 | cut -d ',' -f 1 | tr -d RNA )

while [ $count_file -le $total_rna ]
do
        if [ -f rscapedata/RNA$count_file.cov ] && [ -s rscapedata/RNA$count_file.cov ]  # Check file exists and is not empty
        then
            f=./rscapedata/RNA$count_file.cov
            max=$(grep -r 'GTp' $f | cut -d '[' -f 2 | tr ']' ',' | cut -d ',' -f 2)
            grep -v "#" $f | cut -f 4 > score
            sum=0
            
            ######## Sum over all covariance
            for number in $( cat score ) ; do sum=$(calc $sum+$number) ; done
            echo $max,$sum > info

        else
            ######## If no covariance calculated/no data available
            echo NA,NA > info
        fi

        cat info >> rscape-dataset.csv
        count_file=$(( $count_file + 1 ))

done

rm -rf score
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

############################################################################################################################

# Use IntaRNA to determine RNA:RNA interactions

############################################################################################################################

######## Need to clear any previously set lib path, as otherwise the defined lib path will be appended onto the previous
echo InteractionMIN,InteractionAVE > interaction-intermediate.csv
[ -z "$lib_directory" ] || unset LD_LIBRARY_PATH

for seq in $( grep ">" $initial_fasta )
do
    grep -A 1 $seq $initial_fasta > target.fa
    if [ -z "$lib_directory" ]   # If Boost library didn't have to specified, run as normal.
    then
        $IntaRNA_exe -q target.fa -t $interaction_database > intaRNA-results 2>>errors.log
    else
        LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$lib_directory $IntaRNA_exe -q target.fa -t $interaction_database > intaRNA-results 2>>errors.log
    fi 
    
    ######## Grab all recorded interactions
    grep "interaction energy" intaRNA-results | cut -d '=' -f 2 | tr -d "kcal/mol" | tr -d ' ' > kcal_mol

    min=0    # Min interaction energy
    count=$( grep -c "interaction energy" intaRNA-results )   # Number of interaction energies recorded
    sum=0    # Sum of all energies calculated, prior to averaging

    for number in $( cat kcal_mol )
    do
        ######## If the interaction energy is smaller than the recorded min, update the min
        if (( $( echo "$number < $min" | bc -l ) ))
        then
            min=$number
            sum=$( calc $sum+$number )
        else
            sum=$( calc $sum+$number )
        fi
    done

    ######## If no interactions were recorded, set minimum and average as NA 
    if (( $( echo "$count == 0" | bc -l) ))
    then
        ave=NA
        min=NA
    else
        ave=$( calc $sum/$count )
    fi

    echo $min,$ave >> interaction-intermediate.csv

done

rm -rf target.fa
rm -rf intaRNA-results
rm -rf kcal_mol

############################################################################################################################

# Calculate genomic copy number

############################################################################################################################

######## Run blastn to identify number of matches within human genome
echo Genome_copy_number,Genome_complete_match > copy-number.csv

$blastn_exe -query $initial_fasta -db $human_genome -evalue 0.01 -out output.csv -outfmt "10 qaccver saccver pident" >/dev/null 2>>errors.log

count_file=$( head -2 $initial_data | tail -1 | cut -d ',' -f 1 | tr -d RNA )

while [ $count_file -le $total_rna ]
do
    ######## Process each RNA in the file (by ID) as the counter increases
    grep "RNA$count_file" output.csv > results
    total=$(cat results | wc -l)
    [ -z "$total" ] && echo NA >> copy-number.csv || echo $total >> copy-number.csv
    count_file=$(( $count_file + 1 ))
done

rm -rf results
rm -rf rna.fa
rm -rf output.csv
rm -rf sequences

############################################################################################################################

# Distance to nearest transposable element (Dfam)

############################################################################################################################

######## Format data for bedtools
grep -v Chromosome $initial_data |  cut -f 3,4,5 -d "," | tr ',' '\t' > rnacentraloutput_FINAL.bed
sort -k1,1 -k2,2n rnacentraloutput_FINAL.bed | tr ' ' '\t' > rnacentraloutput_FINAL_sorted.bed 

######## Upstrean hits
$bedtools_exe closest -a rnacentraloutput_FINAL_sorted.bed -b $dfam_hits -io -D ref -iu > dfam_downstream.bed 2>>errors.log

######## Downstream hits
$bedtools_exe closest -a rnacentraloutput_FINAL_sorted.bed -b $dfam_hits -io -D ref -id > dfam_upstream.bed 2>>errors.log

paste <( cut -f 1,2,3,7 dfam_downstream.bed ) <( cut -f 7 dfam_upstream.bed ) --delimiters '\t' > dfam_combined.bed

######## Calculate sum of upstream and downstream, and combine into one file
echo Chromosome,Start,End,Dfam_min,Dfam_sum > dfam-distance.csv

for line in $( cat dfam_combined.bed )
do
    chr=$( echo $line | cut -f 1 | cut -c 4- )
    start=$( echo $line | cut -f 2 )
    end=$( echo $line | cut -f 3 )
    downstream=$( echo $line | cut -f 4 )
    upstream=$( echo $line | cut -f 5 | tr -d '-' )
    [ -z "$downstream" ] && downstream=0    # Implies that no region up or downstream, rather than data missing.
    [ -z "$upstream" ] && upstream=0 
    sum=$( calc $downstream+$upstream )
    if [[ $upstream -lt $downstream ]]    # Taking into account forward and reverse strand
    then
        echo $chr,$start,$end,$upstream,$sum >> dfam-distance.csv
    elif [[ $downstream -lt $upstream ]]
    then
        echo $chr,$start,$end,$downstream,$sum >> dfam-distance.csv
    else
        echo $chr,$start,$end,$downstream,$sum >> dfam-distance.csv
    fi
done

mv dfam_downstream.bed additional-output/
mv dfam_upstream.bed additional-output/
rm -rf dfam_combined.bed

############################################################################################################################

######## Function for matching up dfam-distance to ncRNA in R (bedtools is unordered)

reformat_dfam() {

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

######## Reorder dfam-distance data
cat $initial_data > file1.csv
cat dfam-distance.csv > file2.csv
reformat_dfam >/dev/null 2>>errors.log
cat file3.csv | cut -d ',' -f 7,8 > ncrna-dfam-distance.csv

rm -rf dfam-distance.csv
rm -rf rnacentraloutput_FINAL_sorted1.bed
rm -rf rnacentraloutput_FINAL.bed
rm -rf rnacentraloutput_FINAL_sorted.bed
rm -rf file1.csv
rm -rf file2.csv
rm -rf file3.csv

mv snp-intersection.bed additional-output/
mv snp-intersection-average.csv additional-output/

############################################################################################################################

# Extracting expression for a region across 12 ENCODE datasets

############################################################################################################################

grep -v "Start" $initial_data | cut -d ',' -f 3,4,5 | tr ',' ' ' | perl -lane '{print "$F[0]:$F[1]-$F[2]"}' > locations

echo CountAverage,MaxReadDepth > encode-rnaseq.csv

######## Define location of downloaded ENCODE RNAseq data
encode_data=$( ls $additional_folder/*.bam )

for line in $(cat locations)
do
    ######## Record read counts across all RNAseq files
    count_sum=0 
    for file in $encode_data ; do read=$( $samtools_exe view -c $file $line 2>>errors.log ) ; count_sum=$(( $count_sum+$read )) ; done 

    ######## Determine max read depth
    max_depth=0
    for file in $encode_data
    do 
        max=$( $samtools_exe depth -r $line $file 2>>errors.log | cut -f 3 | sort -V | tail -1 )
        if [[ $max -gt $max_depth ]] ; then max_depth=$max ; else : ; fi
    done

    [ -z "$max_depth" ] && max_depth=0    # record MRD as zero rather than blank

    ######## Determine sequence length
    start=$( echo $line | cut -d ':' -f 2 | cut -d '-' -f 1 )
    end=$( echo $line | cut -d ':' -f 2 | cut -d '-' -f 2 )
    length=$(calc $end-$start)
    
    ######## Calculate average RNAseq counts
    count_average=$( calc $count_sum/12 )
    count_ave=$( calc $count_average/$length )

    if [ -z "$count_average" ] ; then count_average=NA ; max_depth=NA ; else : ; fi
    if [ -z "$count_ave" ] ; then count_ave=NA ; max_dpeth=NA ; else : ; fi

    echo $count_ave,$max_depth >> encode-rnaseq.csv

done

rm -rf locations

############################################################################################################################

# Calculate GC%, MFE, Accessibility and CPC2 Fickett Score

############################################################################################################################

######## GC% calculation

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
        percentage=NA
    else
        GC=$( calc $GC_count/$total)
        percentage=$( calc $GC*100)
    fi
    echo $percentage >> GC.csv
done

rm -rf sequence

######## MFE calculation

echo MFE > MFE.csv

$RNAfold_exe --outfile=rnafold-output $initial_fasta >/dev/null 2>>errors.log

grep "(" rnafold-output | rev | cut -d "(" -f 1 | rev | tr -d ")" | tr -d " " >> MFE.csv
awk '!NF{$0="0"}1' MFE.csv > MFE-final.csv
rm -rf *.ps

######## Accessibility calculation

echo Accessibility > access.csv

for sequence in $( grep -v ">" $initial_fasta )
do 
    access=$( timeout 60m python3 $additional_folder/access_py.py -s $sequence 2>>errors.log ) 
    exit_status=$?
    if [ -z "$access" ]  # If no value calculated, record NA
    then
        echo NA >> access.csv
    elif [[ "$access" == "THIS happened"* ]]   # If error occurred, record NA
    then
        echo NA >> access.csv
    elif [ "$exit_status" -eq "124" ]    # If calculation timed out, record NA
    then
        echo NA >> access.csv
    else
        echo $access >> access.csv
    fi
done

sed -i 's/nan/NA/g' access.csv  # make NA readable by R

rm -rf MFE.csv

######## Fickett score calculation (python script so cannot be added to $PATH)

python2 $cpc2_directory/CPC2.py -i $initial_fasta -o output >/dev/null 2>>errors.log    # Original version
#python3 $cpc2_directory/CPC2.py -i $initial_fasta -o output >/dev/null 2>>errors.log   # Uses updated biopython packages
cat output.txt | cut -f 4 > CPC2.csv
rm -rf output.txt

############################################################################################################################

# Combine all generated data into one file

############################################################################################################################

######## Update data parameter, in case script runs over multiple days
date=$(date +%y%m%d)
name=$( echo $initial_data | cut -d '-' -f 2,3 )

paste -d ',' $initial_data 1000g-freqsummary.csv 1000g-popstats.csv rscape-dataset.csv interaction-intermediate.csv copy-number.csv encode-rnaseq.csv GC.csv MFE-final.csv access.csv rnacode.csv ncrna-dfam-distance.csv conservation.csv CPC2.csv > $date-$name-final-dataset.csv

######## Move previous files to folder, in case there was an issue during feature calculation
mv 1000g-freqsummary.csv additional-output/
mv 1000g-popstats.csv additional-output/
mv rscape-dataset.csv additional-output/
mv interaction-intermediate.csv additional-output/
mv copy-number.csv additional-output/
mv encode-rnaseq.csv additional-output/
mv GC.csv additional-output/
mv MFE-final.csv additional-output/
mv access.csv additional-output/
mv rnacode.csv additional-output/
mv rnafold-output additional-output/
mv ncrna-dfam-distance.csv additional-output/
mv conservation.csv additional-output/
mv CPC2.csv additional-output/

#####################################################################

echo Finished calculating functionality traits.
echo

rm -rf id
rm -rf output
rm -rf test
