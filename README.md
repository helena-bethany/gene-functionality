# Features of functional human genes.

### Supplementary Materials and Methods

* [Required dependencies](#required-dependencies)

* [Retrieval of functional genes](#retrieval-of-functional-genes)

* [Negative control sequences](#negative-control-sequences)

* [Intrinsic sequence features](#intrinsic-sequence-features)

* [Sequence conservation features](#sequence-conservation-features)

* [Transcription features](#transcription-features)

* [Genomic repeat associated features](#genomic-repeat-associated-features)

* [Protein and RNA specific features](#protein-and-rna-specific-features)

* [Population variation features](#population-variation-features)

* [Analysis of potential functionality predictors](#analysis-of-potential-functionality-predictors)

* [References](#references)

-----------------------------------------------------------------------

### Required dependencies

| Software | Program(s) Used | Reference |  
|:------------:|:------------:|:------------:|
|[Bedtools 2.29.0](https://bedtools.readthedocs.io/en/latest/content/installation.html)|bedtools| Quinlan, 2014 |
|[Blast 2.10.0](https://ncbiinsights.ncbi.nlm.nih.gov/2019/12/18/blast-2-10-0/)|blastn ; makeblastdb| Altschul et al., 1990 |
|[CPC2 beta](http://cpc2.cbi.pku.edu.cn/)|CPC2.py| Kang et al., 2017|
|[IntaRNA 3.2.0](https://github.com/BackofenLab/IntaRNA/#install)|IntaRNA|Mann et al., 2017|
|[R-scape 1.4.0](http://eddylab.org/R-scape/)|r-scape|Rivas et al.,2017| 
|[Samtools version 1.10](http://www.htslib.org/download/)|samtools| Li et al., 2009 | 
|[Tabix 1.10](http://www.htslib.org/doc/tabix.html)|tabix|Li, 2011|
|[TIsigner](https://github.com/Gardner-BinfLab/TIsigner_paper_2019)|access_py.py| Bhandari et al., 2019 |
|[UCSC genome browser 'kent' bioinformatic utilities](http://hgdownload.soe.ucsc.edu/admin/exe/)|bigWigSummary ; liftOver ; mafFetch | Haeussler et al., 2019 | 
|[Variant Call Format 4.2](https://vcftools.github.io/downloads.html)|vcftools| Danecek et al., 2011 | 
|[ViennaRNA 2.4.14](https://www.tbi.univie.ac.at/RNA/)|RNAalifold ; RNAcode ; RNAfold ; RNAplfold | Lorenz et al., 2011 | 

-----------------------------------------------------------------------

### Retrieval of functional genes

The links for downloading the functional protein-coding genes, short ncRNA, precursor miRNA and lncRNA FASTA sequences are available below. Additional files and processing instructions reqired are also listed for each dataset.

Link for downloading functional **protein-coding genes** from HGNC:
ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt

* Download text file for all protein-coding genes from HGNC (Braschi et al., 2019), called gene_with_protein_product.txt, then process as described below. The generated text file is then used as one of the input files for [protein-coding-processing.sh](bin/protein-coding-processing.sh).

```
# Pipes one and two removes sequences encoded in the Y and Mt chromosomes or are no longer approved genes.
# Pipes three to five extracts the columns associated with RefSeq IDs, and removes sequences associated with more than one ID.

grep -v "mitochondrially encoded\|Entry Withdrawn\|Yq" gene_with_protein_product.txt | grep "Approved" | cut -f 24 | grep -v "|" | grep . > 200625-hgnc-PC-refseq-ids.txt
```

* To obtain the chromosome coordinates for each of the ncRNA from RNAcentral, the Homo sapiens GRCh38 bed file was obtained from the RNAcentral FTP directory, which contains the chromosome coordinates for all human ncRNA (The RNAcentral Consortium et al., 2017).

Link to RNAcentral FTP directory containing homo_sapiens.GRCh38.bed.gz: ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/genome_coordinates/bed/homo_sapiens.GRCh38.bed.gz

[Link for downloading functional **short ncRNA** from RNAcentral](https://rnacentral.org/search?q=HGNC%20AND%20NOT%20rna_type:%22lncRNA%22%20%20AND%20NOT%20rna_type:%22rRNA%22%20%20AND%20NOT%20rna_type:%22precursor%20RNA%22)

* The downloaded FASTA file is then processed using [short-ncrna-processing.sh](bin/short-ncrna-processing.sh)

[Link for downloading functional **precursor miRNA** from RNAcentral](https://rnacentral.org/search?q=rna_type:%22precursor%20RNA%22%20AND%20expert_db:%22HGNC%22%20AND%20TAXONOMY:%229606%22)

* The downloaded FASTA file is then processed using [short-ncrna-processing.sh](bin/short-ncrna-processing.sh). To create a short ncRNA dataset with 1000 sequences, line 91 of this script can be altered to with the required number.

```
# Line 91 of short-ncrna-processing.sh, to create a dataset of 89 randomly selected precursor miRNAs
shuf -i 1-$max -n 89 > numbers
```

[Link for downloading functional **only lncRNA** from RNAcentral](https://rnacentral.org/search?q=HGNC%20AND%20rna_type:%22lncRNA%22)

* The downloaded FASTA file is then processed using [long-ncrna-processing.sh](bin/long-ncrna-processing.sh)

### Negative control sequences

The links for downloading the required additional files for generating and filtering the negative control sequences, as well as processing steps, are available below. 

[Link for downloading the NCBI GRCh38.p13 human genome FASTA and gff3 files from webserver](https://www.ncbi.nlm.nih.gov/genome/guide/human/)

* The human genome from NCBI (O'Leary et al., 2016) was reformmated to a tab delimited CSV file using the FASTA formatter from [FASTX-Toolkit version 0.0.13](http://hannonlab.cshl.edu/fastx_toolkit/) to allow for easier manipulation and parsing. 

```
# Download and unzip FNA file for human genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz   
gunzip GCF_000001405.39_GRCh38.p13_genomic.fna.gz

# Reformat FNA to tab delimited CSV file
fasta_formatter -i GCF_000001405.39_GRCh38.p13_genomic.fna.gz -o GRCh38_interim.csv -t

# Remove any scaffolds or alternative chromosomes
grep "NC_" GRCh38_interim.csv > GRCh38.p13_genome.csv
```

* Both SwissProt and GENCODE annotations are used to filter negative control sequences which overlap with known protein-coding and ncRNAs, which was done using ```bedtools intersect``` (Haeussler et al., 2019; Harrow et al., 2012; Quinlan, 2014).

Link for downloading SwissProt annotations of known protein-coding genes:
https://hgdownload.soe.ucsc.edu/gbdb/hg38/uniprot/unipAliSwissprot.bb

Link for downloading GENCODE annotations of known ncRNAs:
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gtf.gz

```
# Converting GENCODE GTF to bed file

# Pipe two filters sequences so only known ncRNAs are included.
# Pipes three and four reformat the data into a bed file format.

cat gencode.v34.annotation.gtf | grep "Mt_rRNA\|Mt_tRNA\|miRNA\|misc_RNA\|rRNA\|scRNA\|snRNA\|snoRNA\|ribozyme\|sRNA\|scaRNA\|lncRNA" | awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$10,$16,$7}}' | tr -d '";' > gencode-ncrna-annotation.bed
```
-----------------------------------------------------------------------

### Intrinsic sequence features

GC content was calculated on command line by dividing the number of guanine and cytosine nucleotides by total sequence length.

### Sequence conservation features

The phastCons and phyloP scores for each sequence were calculated from the UCSC phastCons and phyloP 100-way bigWig files using either ```bigWigSummary``` (Haeussler et al., 2019).

Link for downloading phyloP100way scores:
http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/hg38.phyloP100way.bw

Link for downloading phastCons100way scores:
http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phastCons100way/hg38.phastCons100way.bw

```
# Parameters used for running bigWigSummary
# -type=mean reports the mean score
# -type=max reports the maximum score
# $chr $start $end is the chromosome coordiantes for the sequence
# 1 causes the score across the whole sequence to be reported

bigWigSummary -type=mean hg38.phyloP100way.bw $chr $start $end 1
bigWigSummary -type=max hg38.phyloP100way.bw $chr $start $end 1

```

### Transcription features

The links for the ENCODE total RNA-seq data (BAM files) for cell lines and differentiated cells are available below (ENCODE Project Consortium, 2012). Prior to estimating the level of transcription for ncRNA using the downloaded ENCODE data, the BAM files need to be indexed using ```samtools index``` (Li et al., 2009). Read counts were calculated using ```samtools view -c``` and maximum read depth was calculated using ```samtools depth -r``` (Li et al., 2009).

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

### Genomic repeat associated features

Prior to using ```blastn``` to estimate the number of sequence copies in the human genome, as blast database needs to be created from the previously downloaded [GRCh38.p13 human genome](https://www.ncbi.nlm.nih.gov/genome/guide/human/) (Altschul et al., 1990; O’Leary et al., 2016). 

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

Distances to the nearest repetitive element were calculated from a file of non-redundant hits of repetitive DNA elements in the human genome obtained from Dfam v3.1 (Hubley et al., 2016). Once this had been converted to bed file format, the distance to the nearest hit was calculated using ```bedtools closest```  (Quinlan, 2014). 

Link for downloading Dfam non-redundent hg38 repetitive DNA elements:
https://www.dfam.org/releases/Dfam_3.1/annotations/hg38/hg38_dfam.nrph.hits.gz

```
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

### Protein and RNA specific features

##### Coding potential

Coding potential scores were either calculated using from sequence alignments using ```RNAcode``` or individual sequences using ```CPC2.py```, with default parameters being used for both (Kang et al., 2017; Lorenz et al., 2011). 

##### RNA structure

Covariance scores were calculated used ```rscape```, RNAalifold score was calculated using ```RNAalifold```, MFE was calculated with ```RNAfold``` and accessibility was calculated using ```access_py.py``` (Bhandari et al., 2019; Lorenz et al., 2011). Multiple sequence alignments for each sequence were obtained from the UCSC multiz100way alignment using ```mafFetch``` (Haeussler et al., 2019). Unless described below, default parameters were used for each executable.

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

##### RNA:RNA interactions

Interaction energies were calculated using ```IntaRNA```, which were run against a interaction database containing 34 ncRNAs known to interact with a variety of RNAs (Mann et al., 2017; The RNAcentral Consortium et al., 2017). Default parameters ```-q``` query sequence and ```-t``` for target sequences were used for ```IntaRNA```.

[Link for downloading curated interaction database sequences from RNAcentral](https://rnacentral.org/search?q=URS000013F331_9606%20OR%20URS00003D2CC9_9606%20OR%20URS000038803E_9606%20OR%20URS0000735371_9606%20OR%20URS000072E1AF_9606%20OR%20URS00006A14B6_9606%20OR%20URS000065DEBF_9606%20OR%20URS0000C8E9D4_9606%20OR%20URS0000233681_9606%20OR%20URS000047B05D_9606%20OR%20URS00001A72CE_9606%20OR%20URS0000639DBE_9606%20OR%20URS00006D74B2_9606%20OR%20URS0000659172_9606%20OR%20URS000074C9DF_9606%20OR%20URS0000733374_9606%20OR%20URS000067A424_9606%20OR%20URS000007E37F_9606%20OR%20URS000034AAC2_9606%20OR%20URS0000744456_9606%20OR%20URS0000C8E9CE_9606%20OR%20URS00000F9D45_9606%20OR%20URS00006BDF17_9606%20OR%20URS00003EE995_9606%20OR%20URS00003F07BD_9606%20OR%20URS000075BAAE_9606%20OR%20URS000035C796_9606%20OR%20URS000000A142_9606%20OR%20URS000075EF5D_9606%20OR%20URS000075ADBA_9606%20OR%20URS0000443498_9606%20OR%20URS0000103047_9606%20OR%20URS00005CF03F_9606%20OR%20URS0000007D24_9606)

```
# Remove newlines from downloaded RNAcentral fasta file

n=$( grep -c ">" URS000013F331_9606_OR_URS00003D2CC9_9606_OR_URS000038803E_9606_OR_URS0000735371_9606_OR_URS000072E1AF_9606_OR_URS00006A14B6_9606_OR_URS000065DEBF_9606_OR_URS0000C8E9D4_9606_OR_URS0000233681_96_etc.fasta )
num=$(( n*2 ))
cat URS000013F331_9606_OR_URS00003D2CC9_9606_OR_URS000038803E_9606_OR_URS0000735371_9606_OR_URS000072E1AF_9606_OR_URS00006A14B6_9606_OR_URS000065DEBF_9606_OR_URS0000C8E9D4_9606_OR_URS0000233681_96_etc.fasta | awk '/^>/ {printf("\n%s\n",$1);next; } { printf("%s",$1);}  END {printf("\n");}' | tail -"$num" > curated-interaction-database.fa
```

### Population variation features

Prior to downloading the 1000 Genomes Project data, the sequence chromosomal coordinates were converted from hg38 to hg19 using ```liftOver``` and a hg38ToHg19.over.chain conversion file from UCSC (Haeussler et al., 2019). Population data from phase 3 of the 1000 Genomes Project was obtained by downloading VCF files using ```tabix -f -h``` (The 1000 Genomes Project Consortium, 2015; Li et al., 2009). SNP type and frequencies, which were used to caluculate the average minor allele frequency and transitions:transversions ratio, was calculated using ```vcftools``` (Danecek et al., 2011). Tajima's D and Fu and Li's D were calculated using [```popGenome```](https://cran.r-project.org/package=PopGenome) 2.7.5 in R (Pfeifer et al., 2014). Unless described below, default parameters were used for each executable.

Link for downloading the UCSC hg38ToHg19.over.chain conversion file:
http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz

```
# Parameters used for running vcftools
# --vcf specifies that the file input is a VCF
# FASTA/$f is the input VCF file
# --freq specifies that allele frequencies should be calculated
# --out freq-afr is the name of the output file

vcftools --vcf FASTA/$f --freq --out freq-afr
```

### Analysis of potential functionality predictors

Spearman's Rho was calculated in R using the ```cor``` function and the confidence intervals were calculated using ```spearman.ci``` from the package [```spearmanCI```](https://cran.r-project.org/package=spearmanCI) 1.0. The 95% confidence intervals were specified using ```conf.level = 0.95``` and were calculated by boot-strapping 1000 times, which was specified using ```nrep=1000```.  

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

Altschul, S.F., Gish, W., Miller, W., Myers, E.W., and Lipman, D.J. (1990). Basic local alignment search tool. J. Mol. Biol. 215, 403–410.

Bhandari, B.K., Lim, C.S., and Gardner, P.P. (2019). Highly Accessible Translation Initiation Sites Are Predictive of Successful Heterologous Protein Expression.” bioRxiv.

Braschi, B., Denny, P., Gray, K., Jones, T., Seal, R., Tweedie, S., Yates, B., and Bruford, E. (2019). Genenames.org: the HGNC and VGNC resources in 2019. Nucleic Acids Res. 47, D786–D792.

Danecek, P., Auton, A., Abecasis, G., Albers, C.A., Banks, E., DePristo, M.A., Handsaker, R.E., Lunter, G., Marth, G.T., Sherry, S.T., et al. (2011). The variant call format and VCFtools. Bioinformatics 27, 2156–2158.

ENCODE Project Consortium (2012). An integrated encyclopedia of DNA elements in the human genome. Nature 489, 57–74.

Haeussler, M., Zweig, A.S., Tyner, C., Speir, M.L., Rosenbloom, K.R., Raney, B.J., Lee, C.M., Lee, B.T., Hinrichs, A.S., Gonzalez, J.N., et al. (2019). The UCSC Genome Browser database: 2019 update. Nucleic Acids Res. 47, D853–D858.

Harrow, J., Frankish, A., Gonzalez, J.M., Tapanari, E., Diekhans, M., Kokocinski, F., Aken, B.L., Barrell, D., Zadissa, A., Searle, S., et al. (2012). GENCODE: the reference human genome annotation for The ENCODE Project. Genome Res. 22, 1760–1774.

Hubley, R., Finn, R.D., Clements, J., Eddy, S.R., Jones, T.A., Bao, W., Smit, A.F.A., and Wheeler, T.J. (2016). The Dfam database of repetitive DNA families. Nucleic Acids Res. 44, D81–D89.

Kang, Y.-J., Yang, D.-C., Kong, L., Hou, M., Meng, Y.-Q., Wei, L., and Gao, G. (2017). CPC2: a fast and accurate coding potential calculator based on sequence intrinsic features. Nucleic Acids Res. 45, W12–W16.

Karolchik, D., Hinrichs, A.S., Furey, T.S., Roskin, K.M., Sugnet, C.W., Haussler, D., and Kent, W.J. (2004). The UCSC Table Browser data retrieval tool. Nucleic Acids Res. 32, D493–D496.

Li, H. (2011). Tabix: fast retrieval of sequence features from generic TAB-delimited files. Bioinformatics 27, 718–719.

Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., and 1000 Genome Project Data Processing Subgroup (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics 25, 2078–2079.

Lorenz, R., Bernhart, S.H., Höner Zu Siederdissen, C., Tafer, H., Flamm, C., Stadler, P.F., and Hofacker, I.L. (2011). ViennaRNA Package 2.0. Algorithms Mol. Biol. 6, 26.

Mann, M., Wright, P.R., and Backofen, R. (2017). IntaRNA 2.0: enhanced and customizable prediction of RNA-RNA interactions. Nucleic Acids Res. 45, W435–W439.

O’Leary, N.A., Wright, M.W., Brister, J.R., Ciufo, S., Haddad, D., McVeigh, R., Rajput, B., Robbertse, B., Smith-White, B., Ako-Adjei, D., et al. (2016). Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. Nucleic Acids Res. 44, D733–D745.

Pfeifer, B., Wittelsbürger, U., Ramos-Onsins, S.E., and Lercher, M.J. (2014). PopGenome: an efficient Swiss army knife for population genomic analyses in R. Mol. Biol. Evol. 31, 1929–1936.

Quinlan, A.R. (2014). BEDTools: The Swiss-Army Tool for Genome Feature Analysis. Curr. Protoc. Bioinformatics 47, 11.12.1–34.

Rivas, E., Clements, J., and Eddy, S.R. (2017). A statistical test for conserved RNA structure shows lack of evidence for structure in lncRNAs. Nat. Methods 14, 45–48.

Robin, X., Turck, N., Hainard, A., Tiberti, N., Lisacek, F., Sanchez, J.-C., and Müller, M. (2011). pROC: an open-source package for R and S+ to analyze and compare ROC curves. BMC Bioinformatics 12, 77.

Sherry, S.T., Ward, M.H., Kholodov, M., Baker, J., Phan, L., Smigielski, E.M., and Sirotkin, K. (2001). dbSNP: the NCBI database of genetic variation. Nucleic Acids Res. 29, 308–311.
  
The 1000 Genomes Project Consortium (2015). A global reference for human genetic variation. Nature 526, 68–74.  
  
The RNAcentral Consortium, Petrov, A.I., Kay, S.J.E., Kalvari, I., Howe, K.L., Gray, K.A., Bruford, E.A., Kersey, P.J., Cochrane, G., Finn, R.D., et al. (2017). RNAcentral: a comprehensive database of non-coding RNA sequences. Nucleic Acids Res. 45, D128–D134.
