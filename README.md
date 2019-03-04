# ncRNA-functionality

### Formatting RNAcentral data

* [**RNAcentral.sh:**](RNAcentral.sh) reformats the RNAcentral FASTA file and obtains the chromosome coordinates for each sequence.
  * Also generates the random sequences for the null dataset.
  * **Arguments:** $1 is the downloaded FASTA file from RNAcentral and $2 is the name you want to give the dataset (eg: 181203dataset).
  * **Additional Scripts:** [coordinates2.py](coordinate2.py), [RfamCM.py](RfamCM.py) and [random_sequences.py](random_sequences.py).
  * **Additional Files:** Homo_sapiens.GRCh38.bed and GRCh38_genome.csv
  * **Additional Functions:** fasta_formatter and cmscan
  * **Output:** $2_final.csv is the converted dataset for RNAcentral data and $2_spare is the null dataset based off the RNAcentral data.
  * **Functional Values:** CM score and E value associated with CM score.

### UCSC Table Browser
 
* [**UCSCtables.sh:**](UCSCtables.sh) use mysql to obtain the number of snps for each set of coordinates from UCSC.
  * **Arguments**: $1 is the name of the dataset file including location and $2 is the name of the final combined file (eg: 181218dataset3)
  * **Additional Functions:** mysql (to access snp151 table)
  * **Functional Values:** number of snps (UCSC), average can be calculated by dividing my length of sequence, mean phyloP and phastCons score.
  
### NCBI

* [**Filterblastn.sh:**](filterblastn.sh) uses blastn on the command line to blast rna against the human genome.
  * **Arguments:** $1 is the name of the dataset file including location and $2 is the name of the final combined file (eg: 181218dataset3)
  * **Additional Functions:** blastn
  * **Additional Files:** human_genome.* files in /inst/ncbi-blast-2.8.1+/bin (check NCBI emails for how I got these)
  * **Functional Values:** Number of repeats, number of different chromosomes repeats appeared on and the number of repeats that had an approx 100% sequence match.
  
### Calculating CM and Covariance

* [**calculate_CM.sh:**](calculate_CM.sh) uses multiple alignment files from multiz100way (UCSC) to generate a covariance model and score it, as well as calculating the maximum and minimum covariance for the sequence.
  * **Arguments:** $1 is the name of the dataset file including location and $2 is the name of the final combined file (eg: 181218dataset3)
  * **Additional Scripts:** [conversion3.py](conversion3.py)
  * **Additional Files:** None
  * **Additional Functions:** mafFetch, RNAalifold, R-scape and cmbuild.
  * **Output:** $1.cmhmm.txt is the calculated CM and HMM scores and $1cov.csv is the maximum and minimum covariance scores.
  * **Functional Values:** maximum and minimum covariance, number of sequences in MAF (and normalised version), CM and HMM (cmbuild), CM-HMM score.
  
### Population Statistics

* [**extractVCF.sh:**](extractVCF.sh) obtained the VCF files from 1000 genomes
  * **Argument:** $1 is the name of the dataset file including location 
  * **Additional Functions:** tabix
  * **Additional Files:** downloaded VCF files should be in the folder FASTA.
  * **R code:** library(PopGenome)
                genome.class <- readData("FASTA", format="VCF")
                genome.class <- neutrality.stats(genome.class)
                dataset <- get.neutrality(genome.class)[[1]]
                write.csv(dataset, "name.csv")
  * **Functional Values:** Tajima’s D, Fu Li’s F and D, number of SNPs and average for 1000 genomes data (number of segregating sites).

### Transitions and Transversions

* [**filter1000genomes.sh:**](filter1000genomes.sh) loops through all the VCF files downloaded for 1000 genomes and prints the number of transitions and transversions in the data.
  * Also extracts the minimum, maximum and average minor allele frequency.
  * **Arguments:** $1 is the dataset file and location and $2 is the name of the updated file (without file extension).
  * **Additional Files:** VCF files that were downloaded for the population statistics.
  * **Additional Functions:** vcftools
  * Transition:Transversion Ratio = (transition +1)/(transition + transversion + 2)
  * **Functional Values:** Number of transitions, transversions and the ratio of the two, minimum, maximum and average MAF.

### RNA:RNA Interactions

* **interactiontest.fa:** Has about 25% of promoters for protein coding genes, mRNA, ncRNA from NCBI (bias towards miRNA) and ncRNA from the dataset (bias towards snoRNA). 
* Got the 75% of ncRNA not from the dataset from NCBI gene by searching promoter region, mRNA and filtering for non-coding genes.
* [**filtergenes.py:**](filtergenes.py) extract top 250 ncrna from the input file and extract its location and sequence. NOTE: this code uses a lot of memory.
  * **Arguments:** $1 is the downloaded text file from NCBI, $2 is the name of your output
  * **Additional Files:** GRCh38_genome.csv
* [**risearch.sh:**](risearch.sh) takes dataset and calculates the number of interactions with the test dataset.
  * **Arguments:** $1 is the location and name for the dataset and $2 is either RIsearch or RIblast.
  * **Additional Scripts:** [average.py](average.py)
  * **Additional Files:** interaction.suf, test_db
  * **Additional Functions:** RIsearch2, RIblast
  * **Output:** risearch.txt and/or riblast.txt
  * **Functional Values:** minimum and average Energy for RNA interactions for RIsearch2 and RIblast.

### Expression

* Cell Lines: H7-hESC (ENCFF089EWC), HepG2 (ENCFF067CVP), GM12878 (ENCFF893HSY), K562 (ENCFF796BVP).
* Differentiated Cells: smooth muscle cell from H9 (ENCFF369DYD), hepatocyte from H9 (ENCFF907AIK), neural progenitor cell from H9 (ENCFF766ESM), myocyte from LHCN-M2 (ENCFF722EAR), bipolar neuron from GM23338 (ENCFF713UNS), myotube from skeletal muscle myoblast (ENCFF178TTA), hematopoietic multipotent progenitor cell (ENCFF065MVD) and cardiac muscle from RUES2 (ENCFF475WLJ).
* [**expencode.sh:**](expencode.sh) will run through all the BAM files and extract the amount of reads for each set of chromosome coordinates.
  * **Arguments:** $1 is the location and name of the dataset, $2 is the name of the output (excludes file extension). 
  * **Additional Files:** the bam.bai files that have been generated for each BAM file.
  * **Command Line:** samtools index ENCFF089EWC.bam
  * **Additional Functions:** samtools
  * **Output:** $2 contains the number of reads that matched the chromosome coordinates, which has then been averaged out across the number of different cell types/lines used.
  * **Functional Values:** number of reads for each cell type or line as well as the average across all cell types or lines.

### RandomForest

* [**randomforestx.R:**](randomforestx.R) runs random forest for x number of times and exports the appropriate information from it.
  * **Arguments:** location of original data, name of importance file, name of errors file, name of prediction file, number of time to repeat random forest.
  * **Output:** Three CSV files corresponding to the importance of predictors across x models, error rates for predicting functionality across x models and the predictions made summed for all x models.
* [**rftest.R:**](rftest.R) runs random forest x number of times with only a subset of the dataset being used as the testData (ie: can split into predicting specific types of ncRNA).
  * **Arguments:** location of original data that isn’t part of the testData, location of data that is part of the testData (70% will also be included in the trainData), name of prediction file, number of times to repeat random forest.
  * **Output:** One CSV file corresponding to the summed predictions made for all x models.
