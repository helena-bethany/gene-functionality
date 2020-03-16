# Features of non-coding RNA functionality in humans

### Materials and Methods

* [Retrieval of functional ncRNA](#retrieval-of-functional-ncrna)

* [Positive and negative control sequences](#positive-and-negative-control-sequences)

* [Sequence conservation](#sequence-conservation)

* [Population genetics](#population-genetics)

* [Secondary structure conservation](#secondary-structure-conservation)

* [RNA-RNA Interactions](#rna-rna-interactions)

* [Genomic Copy Number](#genomic-copy-number)

* [Transcription](#transcription)

* [Analysis of functionality predictors](#analysis-of-functionality-predictors)

* [References](#references)

-----------------------------------------------------------------------

### Retrieval of functional ncRNA

The links for downloading the ncRNA and lncRNA FASTA sequences are available below:

[Link for downloading functional **ncRNA** from RNAcentral](https://rnacentral.org/search?q=HGNC%20AND%20NOT%20rna_type:%22lncRNA%22%20%20AND%20NOT%20rna_type:%22rRNA%22%20%20AND%20NOT%20rna_type:%22precursor%20RNA%22)

[Link for downloading functional **only lncRNA** from RNAcentral](https://rnacentral.org/search?q=HGNC%20AND%20rna_type:%22lncRNA%22)

To obtain the chromosome coordinates for each of the ncRNA from RNAcentral, the Homo sapiens GRCh38 bed file was obtained from the RNAcentral FTP directory, which contains the chromosome coordinates for all human ncRNA (The RNAcentral Consortium et al., 2017).

Download homo_sapiens.GRCh38.bed.gz from: ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/genome_coordinates/bed/homo_sapiens.GRCh38.bed.gz

### Positive and negative control sequences

**A bunch of writing and links related to a currently non-existent positive control.**

The human genome from NCBI was reformmated to a CSV file using the FASTA formatter from [FASTX-Toolkit version 0.0.13](http://hannonlab.cshl.edu/fastx_toolkit/) to allow for easier manipulation and parsing (Hannon, 2010). The X, Y, mictochondrial genome and scaffolds should then be removed from the CSV file.

[Link for downloading the NCBI GRCh38.p12 human genome](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.38/)

```
# Example of converting a genome file to a CSV

# Convert from FASTA to tabular
fasta_formatter -i GRCh38_latest_genomic.fna -o GRCh38_part1.csv -t

# Remove scaffolds
grep "NC_" GRCh38_part1.csv > GRCh38_genome.csv
```

**Talk about [01-generate-ncrna.sh](Scripts/01-generate-ncrna.sh)**

* The initial dataset was generated using the script 01-generate-ncrna.sh, which takes the RNAcentral FASTA file, chromosome coordinates file and CSV human genome file as input. In summary, this script matches up the chromosome coordinates for each of the functional ncRNA, generates the negative control sequences to represent non-functional ncRNA and calculates the bit score for each ncRNA against curated Rfam CMs. This then generates two output files, which are a CSV file of the functional and non-functional ncRNA training and test data, and a text file of the chromosome coordinates that have been correctly formatted for UCSC Table Browser. 

### Sequence conservation

The file snp151_trimmed_sorted.bed, which is a local copy of release 151 of the dbSNP database (Sherry et al., 2001) is provided, but the code for generating it is as follows:

Download snp151.txt.gz from: ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/

```
# Generating a parasable local copy of the dbSNP database

gunzip snp151.txt.gz
grep  -v chr[A-Z]_* I think > snp151_no_alts
grep -v Chromosome rnacentraloutput_FINAL.csv |  cut -f 3,4,5 -d "," | tr ',' '\t' > rnacentraloutput_FINAL.bed
cut -f 2,3,4 snp151_no_alts > snp151_trimmed.bed
sort -k1,1 -k2,2n snp151_trimmed.bed > snp151_trimmed_sorted.bed 
sort -k1,1 -k2,2n rnacentraloutput_FINAL.bed > rnacentraloutput_FINAL_sorted.bed 

awk '$1="chr"$1' rnacentraloutput_FINAL_sorted.bed | tr ' ' '\t' > rnacentraloutput_FINAL_sorted1.bed 
bedtools intersect -c -a rnacentraloutput_FINAL_sorted1.bed -b snp151_trimmed_sorted.bed -wa -sorted > rnacentraloutput_snp_intersection
```

**Talk about [02-ucsc-dbsnp.sh](Scripts/02-ucsc-dbsnp.sh)**

The phastCons and phyloP conservation scores can be obtained from the [UCSC table browser](http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=701273921_eML3CJEnXm4rsnQmqAuYnNHkZVkV) (Karolchik et al., 2004). Both the maximum and average columns should be obtained for each of phastCons and PhyloP.

```
# Reformat chromosome coordinates using dated output from 01-generate-ncrna.sh
grep -v 'Chromosome' file.csv | cut -d ',' -f 1,2,3 | tr ',' ' ' | perl -lane '{print "$F[0]:$F[1]-$F[2]"}'

# Go onto UCSC table browser and input the following (replace with an image???)
# group: Comparitive Genomes
# track: Conservation
# table: Cons 100 Verts (phyloP100 way) or Cons 100 Verts (phastCons100way)
# region: define regions > submit chromosome coordinates (note only 1000 at a time)
# button to click: Summary/Statistics
```

**This will not return all results, so use ##script name to join the returned SNP statitistics to their corresponding IDs - INSERT SCRAP OF CODE BELOW?**

* The parameters used for cmscan were chosen to allow the analysis to run the same as the Rfam web server, which are described in detail on pages 27-28 of the Infernal 1.1 user guide (Nawrocki and Eddy, 2013). Another parameter was included in addition to these, which was --cpu 20, to increase the number of parallel CPU workers and decrease the time taken for cmscan to run (Nawrocki and Eddy, 2013). An E-value of 10,000 was used during the generation of the negative control sequences, using the parameter -E 10000 in cmscan. This meant that the number of bit scores produced for the negative control sequences could be maximised, since not all the negative control sequences would match well to known ncRNA sequences.

### Population genetics

* The number of overlapping features was using bedtools intersect from BEDtools version 2.29.0 to retrieve the number of SNPs in each ncRNA (Quinlan, 2014). 

* Population data from phase 3 of the 1000 Genome Project was obtained in the script 03-1000genomes.sh, which used tabix from Sequence Alignment/Map Tools (SAMtools) to download VCF files for each set of chromosome coordinates in the dataset (Consortium and The 1000 Genomes Project Consortium, 2015; Li et al., 2009). Both -f and -h were specified with tabix which meant the generated index file was overwritten and not saved for each ncRNA (Li et al., 2009). VCFtools version 0.1.13 was used to analyse and determine the frequency of the SNPs present in the VCF files, with the VCF file specified using --vcf, allele frequency calculated using --freq and output specified using --out (Danecek et al., 2011). 

### Secondary structure conservation

* This was implemented using 05-secondary-structure.sh, with the multiple sequence alignment extracted using mafFetch with parameters hg38 and multiz100way, and then converted into a Stockholm alignment using modules Biopython version 1.73 (Cock et al., 2009; Kent, 2018). The consensus secondary structure predicted using RNAalifold from ViennaRNA Package 2.0, and a CM was generated using cmbuild from Infernal 1.1 once a Stockholm alignment and secondary structure consensus was produced for all ncRNA (Lorenz et al., 2011; Nawrocki and Eddy, 2013). 

* The E-value was set to 100 using the parameter -E 100, to allow covariance to be calculated for the negative control sequences, which were less likely to form a realistic secondary structure, impacting the calculations.

### RNA-RNA Interactions

### Genomic Copy Number

* For makeblastdb, -dbtype nucl was specified for the nucleotide input file and -parse_seqids was used to enable sequence ID parsing (Altschul et al., 1990). The output from blastn was customised to include the query accession, subject accession and percentage of identical matches using -outfmt "10 qaccver saccver pident" (Altschul et al., 1990). 

### Transcription

The links for the ENCODE total RNA-seq data (BAM files) for cell lines and differentiated cells are listed below (ENCODE Project Consortium, 2012).

Cell Lines: 
* H7-hESC ([ENCFF089EWC](https://www.encodeproject.org/files/ENCFF089EWC/))
* HepG2 ([ENCFF067CVP](https://www.encodeproject.org/files/ENCFF067CVP/))
* GM12878 ([ENCFF893HSY](https://www.encodeproject.org/files/ENCFF893HSY/))
* K562 ([ENCFF796BVP](https://www.encodeproject.org/files/ENCFF796BVP/))

Differentiated Cells: 
* Smooth muscle cell from H9 ([ENCFF369DYD](https://www.encodeproject.org/files/ENCFF369DYD/))
* Hepatocyte from H9 ([ENCFF907AIK](https://www.encodeproject.org/files/ENCFF907AIK/))
* Neural progenitor cell from H9 ([ENCFF766WSM](https://www.encodeproject.org/files/ENCFF766WSM/))
* Myocyte from LHCN-M2 ([ENCFF722EAR](https://www.encodeproject.org/files/ENCFF722EAR/))
* Bipolar neuron from GM23338 ([ENCFF713UNS](https://www.encodeproject.org/files/ENCFF713UNS/))
* Myotube from skeletal muscle myoblast ([ENCFF178TTA](https://www.encodeproject.org/files/ENCFF178TTA/))
* Hematopoietic multipotent progenitor cell ([ENCFF065MVD](https://www.encodeproject.org/files/ENCFF065MVD/)) 
* Cardiac muscle from RUES2 ([ENCFF475WLJ](https://www.encodeproject.org/files/ENCFF475WLJ/))

Prior to estimating the level of transcription for ncRNA using the downloaded ENCODE data, the BAM files need to be indexed using samtools (Li et al., 2009).

```
# Example of indexing an ENCODE RNA-seq BAM file

samtools index ENCFF089EWC.bam
```

**Talk about [08-encode-transcription.sh](Scripts/08-encode-transcription.sh)**

* The downloaded BAM files were indexed using samtools index from SAMtools, and then the counts for each set of chromosome coordinates were obtained using samtools view, which can both be implemented in 08-encode-transcription.sh (Li et al., 2009). 

### Analysis of functionality predictors

*  For randomForest, ntree=1000 was specified so each model was created using 1000 decision trees, and proximity=TRUE was specified so the proximity measure was calculated between rows (Liaw et al., 2002). 

---------------------------------------------------------------------------------------

### References

ENCODE Project Consortium (2012). An integrated encyclopedia of DNA elements in the human genome. Nature 489, 57–74.
  
Hannon, G.J. (2010). FASTX-Toolkit. http://<span>hannonlab.cshl.edu<span>/fastx_toolkit/. [November 2018].

Karolchik, D., Hinrichs, A.S., Furey, T.S., Roskin, K.M., Sugnet, C.W., Haussler, D., and Kent, W.J. (2004). The UCSC Table Browser data retrieval tool. Nucleic Acids Res. 32, D493–D496.

Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., and 1000 Genome Project Data Processing Subgroup (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics 25, 2078–2079.

Sherry, S.T., Ward, M.H., Kholodov, M., Baker, J., Phan, L., Smigielski, E.M., and Sirotkin, K. (2001). dbSNP: the NCBI database of genetic variation. Nucleic Acids Res. 29, 308–311.
  
The RNAcentral Consortium, Petrov, A.I., Kay, S.J.E., Kalvari, I., Howe, K.L., Gray, K.A., Bruford, E.A., Kersey, P.J., Cochrane, G., Finn, R.D., et al. (2017). RNAcentral: a comprehensive database of non-coding RNA sequences. Nucleic Acids Res. 45, D128–D134.

----------------------------------------------------------------------------------------

# UNEDITED

### Population genetics

1. Use **extractVCF.sh** to download VCF files (if available) from Phase 3 of the 1000 Genomes Project and puts them in a folder named FASTA (folder should be called this to make PopGenome run easier in R).

* [**extractVCF.sh:**](extractVCF.sh)
  * **Argument:** $1 is the name of the dataset file (including location to the file). 
  * **Additional Functions:** [tabix](http://www.htslib.org/doc/tabix.html).
  * **Output:** VCF files for each ncRNA.

2. Go into R and obtain the population statistics using PopGenome.

```
library(PopGenome)
genome.class <- readData("FASTA", format="VCF")
genome.class <- neutrality.stats(genome.class)
dataset <- get.neutrality(genome.class)[[1]]
write.csv(dataset, "name.csv")     
``` 

3. Attach the PopGenome output to the current dataset file, as well as calculate the snp average.

```
Pretty sure I calculated this in the excel file so I haven't got command line code for this.
snp average is number of snp divided by the ncRNA sequence length.
All NAs that were generated by PopGenome should also be changed to zeros.
```

4. Use **filter1000genomes.sh** to calculate the number of transitions, transversions and relevant MAF in the dataset using the downloaded VCF files.

* [**filter1000genomes.sh:**](filter1000genomes.sh) 
  * **Arguments:** $1 is the dataset file and location and $2 is the name of the updated file (without file extension).
  * **Additional Files:** VCF files in FASTA.
  * **Additional Functions:** [vcftools](http://vcftools.sourceforge.net/)
  * Transition:Transversion Ratio = (transition +1)/(transition + transversion + 2)
  * **Output:** $2.csv is the original dataset with the new data added in.

### Secondary structure conservation

1. Use **calculate_CM.sh** to obtain the multiple sequence alignment files from multiz100way (UCSC). These files are then used to calculate CM and HMM using cmbuild.

* [**calculate_CM.sh:**](calculate_CM.sh) 
  * **Arguments:** $1 is the name of the dataset file (including location to file) and $2 is the name of the final combined file (eg: 181218dataset3)
  * **Additional Scripts:** [conversion3.py](conversion3.py)
  * **Additional Functions:** [mafFetch](http://hgdownload.soe.ucsc.edu/admin/exe/), [RNAalifold](https://www.tbi.univie.ac.at/RNA/), [R-scape](http://eddylab.org/R-scape/) and [cmbuild](http://eddylab.org/infernal/).
  * **Output:** $2.csv is the original dataset with the new data added in.
  
### RNA-RNA Interactions

1. Using NCBI Gene, download the text file for the "most relevant" promoters for protein coding genes, "most relevant" mRNA and "most relevant" ncRNA (eg: type ncRNA into NCBI Gene and download the results, alternatively the links for these are in the script filtergenes.sh). Using **filtergenes.py**, extract the sequences for the first 250 sequences for each of these groups.
To geenerate the random sequences from your dataset, this python code can be run:

```
from Bio import SeqIO
from random import sample
import sys
with open(sys.argv[1]) as f:
    seqs = SeqIO.parse(f, "fasta")
    samps = ((seq.name, seq.seq) for seq in  sample(list(seqs),250))
    for samp in samps:
        print(">{}\n{}".format(*samp))
```

* [**filtergenes.sh:**](filtergenes.py) 
  * **Note:** this code uses a lot of memory.
  * **Arguments:** $1 is the downloaded text file from NCBI, $2 is the name of your output (fasta).
  * **Additional Files:** GRCh38_genome.csv.
  * **Output:** fasta file containing sequence for the first 250 genes.

2. Generate the test interaction database using the 750 genes from NCBI plus 250 random ncRNA from the test dataset. Will need to cat the four fasta files together into one before creating the interaction database.

```
cat promoter.fa protein-coding.fa ncrna-ncbi.fa ncrna-dataset.fa > interactionstest.fa 
./risearch2.x -c interactionstest.fa -o interaction.suf
```

3. Use **risearch.sh** to calculate the number of interactions between the ncRNA dataset and the test interaction database.

* [**risearch.sh:**](risearch.sh)
  * **Arguments:** $1 is the name of the dataset file (including location to file) and $2 is the name of the final combined file (eg: 181218dataset3)
  * **Additional Scripts:** [average.py](average.py)
  * **Additional Files:** interaction.suf
  * **Additional Functions:** [RIsearch2](https://rth.dk/resources/risearch/)
  * **Output:** $2.csv is the original dataset with the new data added in.

### Genomic Copy Number

1. Create human genome database that can be used by command line blastn. First you need to download [GCF_000001405.38_GRCh38.p12_genomic.fna.gz](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.38_GRCh38.p12) and unzip the file. Then format this file into a BLAST database.

```
makeblastdb -in GCF_000001405.38_GRCh38.p12_genomic.fna -dbtype nucl -parse_seqids -out human_genome
```

2. Use **filterblastn.sh** to take the ncRNA and blast them against the GRCh38.p12 reference genome.

* [**Filterblastn.sh:**](filterblastn.sh)
  * **Arguments:** $1 is the name of the dataset file (including location to file) and $2 is the name of the final combined file (eg: 181218dataset3)
  * **Additional Functions:** [blastn](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
  * **Additional Files:** human_genome.*
  * **Output:** $2.csv is the original dataset with the new data added in.



### Analysis of functionality predictors

1. Use **randomforest.R** to run RandomForest for x number of times and exports the importance, error and prediction values for each run.

* [**randomforestx.R:**](randomforest.R)
  * **Arguments:** location of original data, name of importance file, name of errors file, name of prediction file, number of time to repeat random forest.
  * **Output:** Three CSV files corresponding to the importance of predictors across x models, error rates for predicting functionality across x models and the predictions made summed for all x models.

2. Calculate the mean of each of these values (I did this in the original csv file but this could probably be incorporated into the R script) and input back into R to graph the results. 

3. Graph the results in R using ggplot2. Can also create a [correlation plot](prettyCorrelation.R) of all the raw predictor data.

```
#Example: Plot of mean importance values; geom_hline should correspond to the "Start" predictor which represents neutral. Make sure the data is ordered from smallest to largest so that the plot is ordered logically.

ggplot(RF1, aes(x=Importance, y=Predictor))+geom_point()+ylab('')+xlab('Mean Variable Importance for 100 runs')+theme(text=element_text(size=18))+geom_line()+geom_hline(yintercept=17,color="Red")
```

```
#Example: Plot of ncRNA for specific predictors as a boxplot and scatter, which is seperated into functional (blue) and non-functional (red) ncRNA.

allData <- read.csv(“190131nullRF1.csv, header=TRUE)
fun <-data.frame(allData$Functionality, stringsAsFactors=FALSE)
CMs <- scale(allData$CM.scan)
test1<-cbind(fun,CMs)
CMb <- scale(allData$CM.build)
test14 <- cbind(fun,CMb)
test1<-cbind(fun,CMs)
test1$newcolumn <- "CM (CMscan)"
test14$newcolumn <- "CM (CMbuild)"
colnames(test1) <- c("Functionality", "Zscore","Predictor")
colnames(test14) <- c("Functionality", "Zscore","Predictor")
CMcompare <- rbind(test1,test14)
CMcompare$Predictor <- factor(CMcompare$Predictor,  levels=unique(CMcompare$Predictor))

qplot(data=CMcompare, x=Predictor, y=Zscore, geom=c("boxplot"), outlier.shape = NA) + geom_jitter(shape = 1, aes(colour=Functionality))+scale_color_manual(values=c("indianred2","cyan3"))+theme(text=element_text(size=22))+ylim(low=-2.7,high=5)+theme(legend.position="none")
```

**OR**
  
1. Use **rftest.R** to run RandomForest for x number of times using only a subset of the dataset being used as the testData (ie: predict specific types of ncRNA).

* [**rftest.R:**](RFtest.R)
  * **Arguments:** location of original data that wasn’t part of the testData, location of data that is part of the testData (70% will also be included in the trainData), name of prediction file, number of times to repeat random forest.
  * **Output:** One CSV file corresponding to the summed predictions made for all x models.
  
