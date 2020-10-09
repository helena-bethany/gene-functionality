# Features of Functional Human Genes.

* [Required dependencies](#required-dependencies)

* [Retrieval of functional genes](#retrieval-of-functional-genes)

* [Negative control sequences](#negative-control-sequences)

* [Intrinsic sequence features](#intrinsic-sequence-features)

* [Sequence conservation features](#sequence-conservation-features)

* [Transcriptome expression features](#transcriptome-expression-features)

* [Genomic repeat associated features](#genomic-repeat-associated-features)

* [Protein and RNA specific features](#protein-and-rna-specific-features)

* [Population variation features](#population-variation-features)

* [Analysis of potential functionality features](#analysis-of-potential-functionality-features)

* [References](#references)

-----------------------------------------------------------------------

### Required dependencies

| Software | Program(s) Used | Reference |  
|:------------:|:------------:|:------------:|
|[access_py.py](https://github.com/bkb3/openen)|access_py.py| Bhandari et al., 2019 |
|[Bedtools 2.29.0](https://bedtools.readthedocs.io/en/latest/content/installation.html)|bedtools| Quinlan, 2014 |
|[Blast 2.10.0](https://ncbiinsights.ncbi.nlm.nih.gov/2019/12/18/blast-2-10-0/)|blastn ; makeblastdb| Altschul et al., 1990 |
|[CPC2 beta](http://cpc2.cbi.pku.edu.cn/)|CPC2.py| Kang et al., 2017|
|[IntaRNA 3.2.0](https://github.com/BackofenLab/IntaRNA/#install)|IntaRNA|Mann et al., 2017|
|[R-scape 1.4.0](http://eddylab.org/R-scape/)|r-scape|Rivas et al.,2017| 
|[Samtools version 1.10](http://www.htslib.org/download/)|bgzip ; samtools| Li et al., 2009 | 
|[Tabix 1.10](http://www.htslib.org/doc/tabix.html)|tabix|Li, 2011|
|[UCSC genome browser 'kent' bioinformatic utilities](http://hgdownload.soe.ucsc.edu/admin/exe/)|bigWigSummary ; bigWigToBedGraph ; liftOver ; mafFetch | Haeussler et al., 2019 | 
|[Variant Call Format 4.2](https://vcftools.github.io/downloads.html)|vcftools| Danecek et al., 2011 | 
|[ViennaRNA 2.4.14](https://www.tbi.univie.ac.at/RNA/)|RNAalifold ; RNAcode ; RNAfold ; RNAplfold | Lorenz et al., 2011 | 

-----------------------------------------------------------------------

### Retrieval of functional genes

The database IDs, chromosome coordinates and sequence for the final functionally assigned protein-coding, short ncRNA and lncRNA sequences are available in [data](data/).

-----------------------------------------------------------------------

##### Protein-coding (HGNC)

A text file for all protein-coding genes from HGNC (Braschi et al., 2019), called gene_with_protein_product.txt, was downloaded and processed as described below. The resulting text file is then used as one of the input files for [protein-coding-processing.sh](bin/protein-coding-processing.sh).

```
# Download text file for all protein-coding genes from HGNC
wget ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt

# Pipes one and two removes sequences encoded in the Y and Mt chromosomes or are no longer approved genes.
# Pipes three to five extracts the columns associated with RefSeq IDs, and removes sequences associated with more than one ID.

grep -v "mitochondrially encoded\|Entry Withdrawn\|Yq" gene_with_protein_product.txt | grep "Approved" | cut -f 24 | grep -v "|" | grep . > 200625-hgnc-PC-refseq-ids.txt
```
-----------------------------------------------------------------------

##### Short and long ncRNA (RNAcentral)

To obtain the chromosome coordinates for each RNAcentral ncRNA, the Homo sapiens GRCh38 bed file was obtained from the RNAcentral FTP directory, which contains the chromosome coordinates for all human ncRNAs (The RNAcentral Consortium, 2019). The [RNAcentral web interface](https://rnacentral.org/) was used to obtain the ncRNA data and used the following filtering steps:

* Short ncRNA: include HGNC database, exclude precursor miRNA/primary transcript, exclude rRNA, exclude lncRNA.

* Precursor miRNA: include HGNC database, include precursor miRNA/primary transcript only.

* lncRNA: include HGNC database, include lncRNA only.

The downloaded short ncRNA and precursor miRNA FASTA files are then processed using [short-ncrna-processing.sh](bin/short-ncrna-processing.sh) and the lncRNA FASTA file is processed using [long-ncrna-processing.sh](bin/long-ncrna-processing.sh). 

```
# Downloaded homo_sapiens.GRCh38.bed.gz from RNAcentral FTP directory
wget ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/15.0/genome_coordinates/bed/

# Line 91 of short-ncrna-processing.sh can be altered to choose a specific number of precursor miRNAs to include.
shuf -i 1-$max -n 89 > numbers
```

-----------------------------------------------------------------------

### Negative control sequences

The chromosome coordinates and sequences for the final negative control protein-coding, short ncRNA and lncRNA genes are available in [data](data/). To generate these negative control sequences, the [human genome](https://www.ncbi.nlm.nih.gov/genome/guide/human/) from NCBI (O'Leary et al., 2016) was converted to a tab delimited CSV file using the FASTA formatter from [FASTX-Toolkit version 0.0.13](http://hannonlab.cshl.edu/fastx_toolkit/), with scaffolds and alternative chromosomes removed, to allow for easier manipulation and parsing. 

```
# Download and unzip GFF file for human genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz   
gunzip GCF_000001405.39_GRCh38.p13_genomic.gff.gz

# Download and unzip FNA file for human genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz   
gunzip GCF_000001405.39_GRCh38.p13_genomic.fna.gz

# Reformat FNA to tab delimited CSV file
fasta_formatter -i GCF_000001405.39_GRCh38.p13_genomic.fna.gz -o GRCh38_interim.csv -t

# Remove any scaffolds or alternative chromosomes
grep "NC_" GRCh38_interim.csv > GRCh38.p13_genome.csv
```

Both Swiss-Prot and GENCODE annotations are used to filter negative control sequences which overlap with known protein-coding and ncRNAs, which was done using ```bedtools intersect``` (Haeussler et al., 2019; Harrow et al., 2012; Quinlan, 2014).

```
# Downloading Swiss-Prot annotations of known protein-coding genes
wget https://hgdownload.soe.ucsc.edu/gbdb/hg38/uniprot/unipAliSwissprot.bb

# Downloading GENCODE annotations of known ncRNAs
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gtf.gz

# Converting GENCODE GTF to bed file:
# Pipe two filters sequences so only known ncRNAs are included.
# Pipes three and four reformat the data into a bed file format.

cat gencode.v34.annotation.gtf | grep "Mt_rRNA\|Mt_tRNA\|miRNA\|misc_RNA\|rRNA\|scRNA\|snRNA\|snoRNA\|ribozyme\|sRNA\|scaRNA\|lncRNA" | awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$10,$16,$7}}' | tr -d '";' > gencode-ncrna-annotation.bed

```

-----------------------------------------------------------------------

### Intrinsic sequence features

All features were calculated using the custom script [feature-generation.sh](bin/feature-generation.sh). GC content was calculated by dividing the number of guanine and cytosine nucleotides by total sequence length.

-----------------------------------------------------------------------

### Sequence conservation features

The PhastCons and PhyloP scores for each sequence were extracted from the UCSC PhastCons and PhyloP 100-way bigWig files using ```bigWigSummary``` (Haeussler et al., 2019). 

```
# Download PhastCons and PhyloP 100way scores
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/hg38.phyloP100way.bw
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phastCons100way/hg38.phastCons100way.bw

##################################

# Parameters used for running bigWigSummary
# -type=mean reports the mean score
# -type=max reports the maximum score
# $chr $start $end is the chromosome coordiantes for the sequence
# 1 causes the score across the whole sequence to be reported

bigWigSummary -type=mean hg38.phyloP100way.bw $chr $start $end 1
bigWigSummary -type=max hg38.phyloP100way.bw $chr $start $end 1

```

GERP scores for each sequence were extracted from the Ensembl 111-mammal GERP bigWig file using ```bedtools map```, once the bigWig file had been converted into individual bed files for each chromosome (Quinlan, 2014; Yates et al., 2020).

```
# Download 111-mammal GERP scores
curl -O ftp://ftp.ensembl.org/pub/release-101/compara/conservation_scores/111_mammals.gerp_conservation_score/gerp_conservation_scores.homo_sapiens.GRCh38.bw

# Convert file from bigWig to bedGraph using bigWigToBedGraph (UCSC Kent Tools)
bigWigToBedGraph gerp_conservation_scores.homo_sapiens.GRCh38.bw gerp_111_mammals_GRCh38.bedgraph

# Order bedgraph file
sort -k1,1 -k2,2n gerp_111_mammals_GRCh38.bedgraph | tr ' ' '\t' > gerp_111_mammals_GRCh38_sorted.bedgraph

# Split bedgraph by chromosomes, and input into a separate folder
awk '{print $0 >> "gerp"$1".bed"}' gerp_111_amniotes_GRCh38_sorted.bedgraph
mkdir gerp-mammals-index && mv gerp*.bed gerp-mammals-index/

# bgzip each bed file and index with tabix to decrease runtime for bedtools
IFS=$'\n' && for file in $( ls gerp-mammals-index/*.bed ) ; do bgzip $file ; tabix -p bed $file.gz ; done

##################################

# Parameters for running bedtools map
# -a is the input file, containing the sequence coordinates
# -b is the gerp bed file for the chromosome being analysed
# -c 4,4 specifies that the mean and max should be calculated using the 4th column, which contains the per base GERP scores
# -o mean,max causes bedtools to output the max and mean GERP score observed for the sequence
# -sorted tells bedtools the bed file inputs are sorted, to reduce runtime

bedtools map -a input.bed -b gerp-mammals-index/gerp$chr.bed.gz -c 4,4 -o mean,max -sorted 

```
-----------------------------------------------------------------------

### Transcriptome Expression features

The links for the ENCODE total and small RNA-Seq data (BAM files) for [primary-cell](data/encode-primary-cell-rnaseq.txt) and [tissue](data/encode-tissue-rnaseq.txt) samples are available in the [data](data/) folder (ENCODE Project Consortium, 2012), with all the RNA-Seq data kept in sample specific folders. Prior to estimating the level of transcription, the BAM files need to be indexed using ```samtools index``` and the total reads per file calculated into a text file in the same folder (Li et al., 2009). Read counts were calculated using ```samtools view -c``` and maximum read depth was calculated using ```samtools depth -r``` (Li et al., 2009).

```
# Downloading ENCODE RNA-Seq datasets (must be done within a specified folder for each sample type)
mkdir ENCODE-tissue && mkdir ENCODE-primary-cell
xargs -L 1 curl -O -J -L < encode-tissue-rnaseq.txt         # Done within ENCODE-tissue/
xargs -L 1 curl -O -J -L < encode-primary-cell-rnaseq.txt   # Done within ENCODE-primary-cell/

# Index files and calculate total reads per RNA-Seq dataset (done within each specified sample type folder)
for file in $( ls *.bam ) ; do samtools index $file ; total=$( samtools view $file | wc -l ) ; echo $file,$total >> total-read-count.txt ; echo $file ; done
```
-----------------------------------------------------------------------

### Genomic repeat associated features

Prior to using ```blastn``` to estimate the number of sequence copies in the human genome, as local blast database was created from the previously downloaded [GRCh38.p13 human genome](https://www.ncbi.nlm.nih.gov/genome/guide/human/) (Altschul et al., 1990; O’Leary et al., 2016). 

```
# Create blast database for blastn
makeblastdb -in GCF_000001405.39_GRCh38.p13_genomic.fna -dbtype nucl -parse_seqids -out human_genome

##################################

# Parameters used for running blastn:
# -query $initial_fasta specifies the sequence input
# -db human_genome specifies the previously generated blastn database
# -evalue 0.01 sets the e-value cut-off to 0.01
# -out output.csv names the generated output file
# -outfmt "10 qaccver saccver pident" specifies that the output should include the query accession, subject accession and percentage of identical matches

blastn -query $initial_fasta -db human_genome -evalue 0.01 -out output.csv -outfmt "10 qaccver saccver pident"
```

Distances to the nearest repetitive elements were calculated from a file of non-redundant hits of repetitive DNA elements in the human genome obtained from Dfam v3.1 (Hubley et al., 2016). Once this was converted to bed file format, the distance to the nearest hit was calculated using ```bedtools closest```  (Quinlan, 2014). 

```
# Download Dfam non-redundent hg38 repetitive DNA elements
wget https://www.dfam.org/releases/Dfam_3.1/annotations/hg38/hg38_dfam.nrph.hits.gz

# Convert Dfam hits into bed format

# Pipe one removes the header
# Pipe two obtains Chr, Start and End columns
# Pipe three removes alternate chromosome, incomplete scaffolds and chromosomes Y/Mt
grep -v "#seq_name" hg38_dfam.nrph.hits | cut -f 1,10,11 | grep -v "_random\|chrY\|chrM" > dfam-hg38-nrph.bed

# For loop reformats into bed format, and sorts the file so it can be used by bedtools
for line in $( cat dfam-hg38-nrph.bed ) ; do chr=$( echo $line | cut -f 1 ) ; start=$( echo $line | cut -f 2 ) ; end=$( echo $line | cut -f 3 ) ; if [ "$start" -gt "$end" ] ; then echo -e $chr'\t'$end'\t'$start >> dfam-hg38.bed ; else echo $line >> dfam-hg38.bed ; fi ; done && sort -k1,1 -k2,2n dfam-hg38.bed > dfam-hg38-sorted.bed

##################################

# Parameters used for running bedtools:
# bedtools closest reports the either overlapping or the nearest sequence in input 
# -a $input_bed is the sequences of interest, converted to a bed file
# -b dfam-hg38-nrph.bed is the Dfam non-redundent repetitive DNA elements
# -io ignores overlaps, as some functional ncRNAs are considered repetitive DNA elements
# -D ref reports distance with respect to the reference genome
# -iu ignores upstream to report closest downstream element
# -id ignores downstream to report closest upstream element

bedtools closest -a $input_bed -b dfam-hg38-sorted.bed -io -D ref -iu 
bedtools closest -a $input_bed -b dfam-hg38.sorted.bed -io -D ref -id 
```
-----------------------------------------------------------------------

### Protein and RNA specific features

##### Coding potential

Coding potential scores were either calculated from sequence alignments using ```RNAcode``` or individual sequences using ```CPC2.py```, with default parameters being used for both (Kang et al., 2017; Lorenz et al., 2011). 

-----------------------------------------------------------------------

##### RNA structure

Covariance scores and associated E-values were calculated used ```rscape```, RNAalifold scores were calculated using ```RNAalifold```, MFE was calculated with ```RNAfold``` and accessibility was calculated using ```access_py.py``` (Bhandari et al., 2019; Lorenz et al., 2011). Multiple sequence alignments for each sequence were obtained from the UCSC multiz100way alignment using ```mafFetch``` (Haeussler et al., 2019). Unless described below, default parameters were used for each executable.

```
# Parameters used for running mafFetch
# hg38 specifies the original species genome
# multiz100way specifies the original multiple sequence alignment to be used
# overBed is the input
# mafOut is the name of the output file

mafFetch hg38 multiz100way overBed mafOut

##################################

# Parameters used for running R-scape
# -E 100 sets the E-value to 100, which was used to ensure as many hits were obtained as possible for the negative control sequences or sequences with little covariance.
# -s $rna_id.stk specifies that the input is a Stockholm file with $rna_id being the ID for a specific sequence

rscape -E 100 -s $rna_id.stk 

##################################

# Parameters used for running RNAalifold
# -q prevents output from being printed
# -f S specifies the input is in Stockholm format
# --noPS prevents postscript files from being created
# RNA.stk is the input file

RNAalifold_exe -q -f S --noPS RNA.stk 
```
-----------------------------------------------------------------------

##### RNA:RNA interactions

Interaction energies were calculated using ```IntaRNA```, which were run against a interaction database containing [34 RNAcentral v15 ncRNAs](https://rnacentral.org/search?q=URS000013F331_9606%20OR%20URS00003D2CC9_9606%20OR%20URS000038803E_9606%20OR%20URS0000735371_9606%20OR%20URS000072E1AF_9606%20OR%20URS00006A14B6_9606%20OR%20URS000065DEBF_9606%20OR%20URS0000C8E9D4_9606%20OR%20URS0000233681_9606%20OR%20URS000047B05D_9606%20OR%20URS00001A72CE_9606%20OR%20URS0000639DBE_9606%20OR%20URS00006D74B2_9606%20OR%20URS0000659172_9606%20OR%20URS000074C9DF_9606%20OR%20URS0000733374_9606%20OR%20URS000067A424_9606%20OR%20URS000007E37F_9606%20OR%20URS000034AAC2_9606%20OR%20URS0000744456_9606%20OR%20URS0000C8E9CE_9606%20OR%20URS00000F9D45_9606%20OR%20URS00006BDF17_9606%20OR%20URS00003EE995_9606%20OR%20URS00003F07BD_9606%20OR%20URS000075BAAE_9606%20OR%20URS000035C796_9606%20OR%20URS000000A142_9606%20OR%20URS000075EF5D_9606%20OR%20URS000075ADBA_9606%20OR%20URS0000443498_9606%20OR%20URS0000103047_9606%20OR%20URS00005CF03F_9606%20OR%20URS0000007D24_9606) known to interact with a variety of RNAs (Mann et al., 2017; The RNAcentral Consortium, 2019). Default parameters ```-q``` query sequence and ```-t``` for target sequences were used for ```IntaRNA```.

```
# Remove newlines from downloaded RNAcentral fasta file

n=$( grep -c ">" URS000013F331_9606_OR_URS00003D2CC9_9606_OR_URS000038803E_9606_OR_URS0000735371_9606_OR_URS000072E1AF_9606_OR_URS00006A14B6_9606_OR_URS000065DEBF_9606_OR_URS0000C8E9D4_9606_OR_URS0000233681_96_etc.fasta )
num=$(( n*2 ))
cat URS000013F331_9606_OR_URS00003D2CC9_9606_OR_URS000038803E_9606_OR_URS0000735371_9606_OR_URS000072E1AF_9606_OR_URS00006A14B6_9606_OR_URS000065DEBF_9606_OR_URS0000C8E9D4_9606_OR_URS0000233681_96_etc.fasta | awk '/^>/ {printf("\n%s\n",$1);next; } { printf("%s",$1);}  END {printf("\n");}' | tail -"$num" > curated-interaction-database.fa
```
-----------------------------------------------------------------------

### Population variation features

Prior to downloading the 1,000 Genomes Project (1kGP) data, the chromosome coordinates were converted from hg38 to hg19 using ```liftOver``` and the hg38ToHg19.over.chain conversion file from UCSC (Haeussler et al., 2019). Population data from phase 3 of the 1kGP was obtained by downloading VCF files using ```tabix -f -h``` (The 1000 Genomes Project Consortium, 2015; Li et al., 2009). SNP frequencies, which were used to caluculate the average minor allele frequency, were calculated using ```vcftools``` (Danecek et al., 2011). All chromosome SNP data from the Genome Aggregation Database (gnomAD) v3 was downloaded as a local copy, with the SNPs for each region extracted using ```tabix -f -h``` (Karczewski et al. 2020; Li et al., 2009).

```
# Download UCSC hg38ToHg19.over.chain conversion file
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz

# Download all chromosome gnomAD v3 VCF file
curl -O https://storage.googleapis.com/gnomad-public/release/3.0/vcf/genomes/gnomad.genomes.r3.0.sites.vcf.bgz

##################################

# Parameters used for running vcftools
# --vcf specifies that the file input is a VCF
# FASTA/$f is the input VCF file
# --freq specifies that allele frequencies should be calculated
# --out freq-afr is the name of the output file

vcftools --vcf FASTA/$f --freq --out freq-afr
```
-----------------------------------------------------------------------

### Analysis of potential functionality features

Spearman's Rho was calculated in R 3.6.0 using the ```cor``` function and the confidence intervals were calculated using ```spearman.ci``` from the package [```spearmanCI```](https://cran.r-project.org/package=spearmanCI) 1.0. The 95% confidence intervals were specified using ```conf.level = 0.95``` and calculated by boot-strapping 1,000 times, which was specified using ```nrep=1000```.  

RandomForest was run using [random-forest-analysis.R](bin/random-forest-analysis.R), which requires the dependecies [```randomForest```](https://cran.r-project.org/package=randomForest) 4.6-14, [```epiR```](https://cran.r-project.org/package=epiR) 1.0-15, [```pROC```](https://cran.r-project.org/package=pROC) 1.16.2 and [```mltools```](https://cran.r-project.org/package=mltools) 0.3.5 (Liaw et al., 2002; Robin et al., 2011).

```
# Parameters used for running randomForest
# Functional~ specifies that randomForest is classifying on whether a sequence is functional or not
# data=trainData specifies the training dataset
# ntree=1000 specifies that 1000 trees should be generated
# proximity=TRUE calculates a proximity matrix
# na.action=na.roughfix allows for sequences with NA values to be used. These missing values are approximated using the median feature values. 

randomForest(Functional~.,data=trainData,ntree=1000,proximity=TRUE, na.action=na.roughfix) 
```

---------------------------------------------------------------------------------------

### References

Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. 1990. Basic local alignment search tool. J Mol Biol 215: 403–410.

Bhandari BK, Lim CS, Gardner PP. 2019. Highly Accessible Translation Initiation Sites Are Predictive of Successful Heterologous Protein Expression.” bioRxiv.

Braschi B, Denny P, Gray K, Jones T, Seal R, Tweedie S, Yates B, Bruford E. 2019. Genenames.org: the HGNC and VGNC resources in 2019. Nucleic Acids Res 47: D786–D792.

Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA, Handsaker RE, Lunter G, Marth GT, Sherry ST, et al. 2011. The variant call format and VCFtools. Bioinformatics 27: 2156–2158.

ENCODE Project Consortium. 2012. An integrated encyclopedia of DNA elements in the human genome. Nature 489: 57–74.

Haeussler M, Zweig AS, Tyner C, Speir ML, Rosenbloom KR, Raney BJ, Lee CM, Lee BT, Hinrichs AS, Gonzalez JN, et al. 2019. The UCSC Genome Browser database: 2019 update. Nucleic Acids Res 47: D853–D858.

Harrow J, Frankish A, Gonzalez JM, Tapanari E, Diekhans M, Kokocinski F, Aken BL, Barrell D, Zadissa A, Searle S, et al. 2012. GENCODE: the reference human genome annotation for The ENCODE Project. Genome Res 22: 1760–1774.

Hubley R, Finn RD, Clements J, Eddy SR, Jones TA, Bao W, Smit AFA, Wheeler TJ. 2016. The Dfam database of repetitive DNA families. Nucleic Acids Res 44: D81–9.

Kang Y-J, Yang D-C, Kong L, Hou M, Meng Y-Q, Wei L, Gao G. 2017. CPC2: a fast and accurate coding potential calculator based on sequence intrinsic features. Nucleic Acids Res 45: W12–W16.

Karczewski KJ, Francioli LC, Tiao G, Cummings BB, Alföldi J, Wang Q, Collins RL, Laricchia KM, Ganna A, Birnbaum DP, et al. 2020. The mutational constraint spectrum quantified from variation in 141,456 humans. Nature 581: 434–443.

Li H. 2011. Tabix: fast retrieval of sequence features from generic TAB-delimited files. Bioinformatics 27: 718–719.

Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, 1000 Genome Project Data Processing Subgroup. 2009. The Sequence Alignment/Map format and SAMtools. Bioinformatics 25: 2078–2079.

Lorenz R, Bernhart SH, Höner Zu Siederdissen C, Tafer H, Flamm C, Stadler PF, Hofacker IL. 2011. ViennaRNA Package 2.0. Algorithms Mol Biol 6: 26.

Mann M, Wright PR, Backofen R. 2017. IntaRNA 2.0: enhanced and customizable prediction of RNA-RNA interactions. Nucleic Acids Res 45: W435–W439.

O’Leary NA, Wright MW, Brister JR, Ciufo S, Haddad D, McVeigh R, Rajput B, Robbertse B, Smith-White B, Ako-Adjei D, et al. 2016. Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. Nucleic Acids Res 44: D733–45.

Quinlan AR. 2014. BEDTools: The Swiss-Army Tool for Genome Feature Analysis. Curr Protoc Bioinformatics 47: 11.12.1–34.

Rivas E, Clements J, Eddy SR. 2017. A statistical test for conserved RNA structure shows lack of evidence for structure in lncRNAs. Nat Methods 14: 45–48.

Robin X, Turck N, Hainard A, Tiberti N, Lisacek F, Sanchez JC, Müller M. 2011. pROC: an open-source package for R and S+ to analyze and compare ROC curves. BMC Bioinformatics 12: 77.
  
The 1000 Genomes Project Consortium. 2015. A global reference for human genetic variation. Nature 526: 68–74. http://dx.doi.org/10.1038/nature15393.
  
The RNAcentral Consortium. 2019. RNAcentral: a hub of information for non-coding RNA sequences. Nucleic Acids Res 47: D221–D229.

Yates AD, Achuthan P, Akanni W, Allen J, Allen J, Alvarez-Jarreta J, Amode MR, Armean IM, Azov AG, Bennett R, et al. 2020. Ensembl 2020. Nucleic Acids Res 48: D682–D688.
