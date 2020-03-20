# Set-up of local databases and script usage

* [Required software and programs](#required-software-and-programs)

* [Generation of local databases](#generation-of-local-databases)

* [Script descriptions and parameters](#script-descriptions-and-parameters)

* [References](#references)

-----------------------------------------------------------------------

## Required software and programs

| Software | Program | Additional Notes |    
|:------------:|:------------:|:-----------:|
|[Bedtools 2.29.0](https://bedtools.readthedocs.io/en/latest/content/installation.html)|bedtools|&nbsp;|
|[Blast 2.10.0](https://ncbiinsights.ncbi.nlm.nih.gov/2019/12/18/blast-2-10-0/)|blastn|&nbsp;|
|[Dfam 3.1](https://www.dfam.org/help/tools)|dfamscan.pl| HMMER is required for perl script to run. |
|[EDirect 13.5](https://www.ncbi.nlm.nih.gov/books/NBK179288/)|esearch;efetch|&nbsp;| 
|[EMBOSS 6.6.0.0](https://www.ebi.ac.uk/seqdb/confluence/display/THD/EMBOSS+Transeq)|transeq|&nbsp;|
|[FASTX-Toolkit 0.0.13](http://hannonlab.cshl.edu/fastx_toolkit/)|fasta_formatter|&nbsp;|
|[HMMER 3.3](http://hmmer.org/download.html)|hmmscan; hmmpress|&nbsp;|
|[Infernal 1.1](https://github.com/EddyRivasLab/infernal)|cmscan; cmpress|&nbsp;|
|[PopGenome 2.7.5](https://cran.r-project.org/web/packages/PopGenome/index.html)|neutrality.stats; get.neutrality|This is a R library and should be installed accordingly.|
|[RIsearch 2.1](https://rth.dk/resources/risearch/#download)|risearch2.x|&nbsp;|
|[R-scape 1.4.0](http://eddylab.org/R-scape/)|r-scape|&nbsp;|
|[Samtools version 1.10](http://www.htslib.org/download/)|samtools; tabix|A previous version of tabix can't be used as the downloaded VCF files won't be formatted correctly.|
|[UCSC genome browser 'kent' bioinformatic utilities](http://hgdownload.soe.ucsc.edu/admin/exe/)|mafFetch|Requires a [.hg.conf file](http://genome.ucsc.edu/goldenPath/help/mysql.html) to be set up in the home directory.|
|[Variant Call Format 4.2](https://vcftools.github.io/downloads.html)|vcftools|&nbsp;|
|[ViennaRNA 2.4.14](https://www.tbi.univie.ac.at/RNA/)|RNAalifold|&nbsp;|

-----------------------------------------------------------------------

## Generation of local databases 

#### Converting human genome from .fna to .csv: 

The conversion of the human genome to a .csv file allows for incomplete scaffolds to be removed and for easier parsing when generating the null sequences.

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz   
gunzip GCF_000001405.39_GRCh38.p13_genomic.fna.gz
fasta_formatter -i GCF_000001405.39_GRCh38.p13_genomic.fna.gz -o GRCh38_interim.csv -t
grep "NC_" GRCh38_interim.csv > GRCh38.p13_genome.csv
```

#### Formatting of Rfam, Pfam and Dfam databases into binary

In order for downloaded CM and HMM files to run locally with ```cmscan``` or ```hmmscan```, the databases need to be converted into binary using ```cmpress``` or ```hmmpress``` respectively.

```bash
# Rfam
wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/14.1/Rfam.cm.gz
gunzip Rfam.cm.gz
cmpress Rfam.cm

# Pfam
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz
gunzip Pfam-A-hmm.gz
hmmpress Pfam-A.hmm

# Dfam
wget https://www.dfam.org/releases/current/families/Dfam.hmm.gz
gunzip Dfam.hmm.gz
hmmpress Dfam.hmm
```

#### Local copy of the dbSNP build 151

A local copy of the dbSNP database from UCSC was downloaded as it was faster to parse than obtaining the data for each region of interest from the UCSC SQL database. The SNP data was trimmed to remove any scaffolds regions and then ordered. 

```bash
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp151.txt.gz
gunzip snp151.txt.gz
grep -w "chr1\|chr2\|chr3\|chr4\|chr5\|chr6\|chr7\|chr8\|chr9\|chr10\|chr11\|chr12\|chr13\|chr14\|chr15\|chr16\|chr17\|chr18\|chr19\|chr20\|chr21\|chr22\|chrX" snp151.txt > snp151_no_alts
cut -f 2,3,4 snp151_no_alts > snp151_trimmed.bed
sort -k1,1 -k2,2n snp151_trimmed.bed > snp151_trimmed_sorted.bed 
```

#### Generating a test interaction database for RIsearch2

The code below is only an example, meaning the following file names will need to altered to your accordingly named files.

```bash
IFS=$'\n'

# Function for getting sequence data from downloaded NCBI files
get_ncbi_sequence() {

        grep "NC_" $1 > file
        max=$( cat file | wc -l )
        shuf -i 1-$max -n 250 > numbers
        id=$( echo $1 | cut -d '.' -f 1 )

        for line in $( cat numbers )
        do
                info=$( sed -n "$line"p file )
                chr=$( echo $info | cut -f 11 )
                start=$( echo $info | cut -f 13 )
                end=$( echo $info | cut -f 14 )
                sequence=$( grep -w "chromosome $chr" $2 | cut -f 2 | cut -c$start-$end )
                echo ">$chr-$start-$end" >> $id.fa
                echo $sequence >> $id.fa
        done

        rm -rf file
        rm -rf number

}

# The two inputs are the downloaded NCBI text file and .csv of the human genome
get_ncbi_sequence mrna.txt GRCh38.p13_genome.csv
get_ncbi_sequence promoter_region.txt GRCh38.p13_genome.csv
get_ncbi_sequence ncrna.txt GRCh38.p13_genome.csv

# Obtain 250 RNA sequences from the short ncRNA and lncRNA datasets (originally from HGNC/RNAcentral)
grep -v ">" short-ncrna/200313-functional-ncrna-seq.fa | shuf -n 250 > ncrna-250
for line in $( cat ncrna-250) ; do grep -B 1 $line short-ncrna/200313-functional-ncrna-seq.fa  >> short-ncrna.fa ; done

grep -v ">" long-ncrna/200316-functional-ncrna-seq.fa | shuf -n 250 > ncrna-250
for line in $( cat ncrna-250) ; do grep -B 1 $line long-ncrna/200316-functional-ncrna-seq.fa  >> long-ncrna.fa ; done

rm -rf ncrna-250

# Combine all sequences into one .fa file and create the interaction database
cat mrna.fa promoter_region.fa ncrna.fa short-ncrna.fa long-ncrna.fa > interaction_test.fa

./risearch2.x -c interaction_test.fa -o interaction_database.suf

```

#### Generating a local BLAST database for the human genome

Note that .fna file would have already been downloaded when creating the human genome .csv file

```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz   
gunzip GCF_000001405.39_GRCh38.p13_genomic.fna.gz
makeblastdb -in GCF_000001405.39_GRCh38.p13_genomic.fna -dbtype nucl -parse_seqids -out human_genome
```

#### Download and indexing RNA-seq BAM files from ENCODE

All .bam files need to be indexed before they can be parsed by samtools to determine the level of transcription for a region of interest.

```bash
# Example of indexing an ENCODE RNA-seq .bam file

wget https://www.encodeproject.org/files/ENCFF089EWC/@@download/ENCFF089EWC.bam
samtools index ENCFF089EWC.bam
```

-----------------------------------------------------------------------

## Script descriptions and parameters

#### [generate-ncrna-filter.sh](Scripts/generate-ncrna-filter.sh)

This script takes .fasta file of functional ncRNA directly downloaded from RNAcentral, the RNAcentral chromosome coordinates file, the .csv file of the human genome and the UCSC uniprot bed file as input. The output will be a .csv and .fa file of at most 1000 functional ncRNA and null sequences, which make up the negative control dataset.

#### [download-refseq.sh](Scripts/download-refseq.sh)

This script takes the text file of RefSeq IDs from HGNC and the .gff file from the lastest version of the human genome as input. The output will be a .csv and .fa file of 1000 protein-coding genes, creating the positive control dataset.

#### [generate-functionality-features.sh](Scripts/generate-functionality-features.sh)

This script prompts the input of the name of the .csv and .fa files generated by the previous two scripts, in addition to named local databases and the paths to the folders containing the local databases and installed dependencies. 

```bash
Initial parameter set-up:

Enter CSV file for initial dataset: initial_data

Enter FASTA file for initial dataset: initial_fasta

________________________________________________________________________________________

Enter folder location containing additional input: additional_folder

Type in the names of the following input within $additional_folder:

Rfam database (cm): rfam_name

Pfam database (hmm): pfam_name

Dfam database (hmm): dfam_name

RNA:RNA interaction database (suf): interaction_database

Human genome blastn database: human_genome

Processed snp151 database (bed): snp_database

________________________________________________________________________________________

Enter folder location containing additional dependencies: additional_dependencies

```

-----------------------------------------------------------------------

## References

Alkan, F., Wenzel, A., Palasca, O., Kerpedjiev, P., Rudebeck, A.F., Stadler, P.F., Hofacker, I.L., and Gorodkin, J. (2017). RIsearch2: suffix array-based large-scale prediction of RNA–RNA interactions and siRNA off-targets. Nucleic Acids Res. 45, e60–e60.

Altschul, S.F., Gish, W., Miller, W., Myers, E.W., and Lipman, D.J. (1990). Basic local alignment search tool. J. Mol. Biol. 215, 403–410.

Danecek, P., Auton, A., Abecasis, G., Albers, C.A., Banks, E., DePristo, M.A., Handsaker, R.E., Lunter, G., Marth, G.T., Sherry, S.T., et al. (2011). The variant call format and VCFtools. Bioinformatics 27, 2156–2158.

Hannon, G.J. (2010). FASTX-Toolkit. http://<span>hannonlab.cshl.edu<span>/fastx_toolkit/. 

Kans, J. (2020). Entrez Direct: E-utilities on the UNIX Command Line (National Center for Biotechnology Information (US)).

Karolchik, D., Hinrichs, A.S., Furey, T.S., Roskin, K.M., Sugnet, C.W., Haussler, D., and Kent, W.J. (2004). The UCSC Table Browser data retrieval tool. Nucleic Acids Res. 32, D493–D496.

Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., and 1000 Genome Project Data Processing Subgroup (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics 25, 2078–2079.

Lorenz, R., Bernhart, S.H., Höner Zu Siederdissen, C., Tafer, H., Flamm, C., Stadler, P.F., and Hofacker, I.L. (2011). ViennaRNA Package 2.0. Algorithms Mol. Biol. 6, 26.

Nawrocki, E.P., and Eddy, S.R. (2013). Infernal 1.1: 100-fold faster RNA homology searches. Bioinformatics 29, 2933–2935.

Pfeifer, B., Wittelsbürger, U., Ramos-Onsins, S.E., and Lercher, M.J. (2014). PopGenome: an efficient Swiss army knife for population genomic analyses in R. Mol. Biol. Evol. 31, 1929–1936.

Quinlan, A.R. (2014). BEDTools: The Swiss-Army Tool for Genome Feature Analysis. Curr. Protoc. Bioinformatics 47, 11.12.1–34.

Rivas, E., Clements, J., and Eddy, S.R. (2017). A statistical test for conserved RNA structure shows lack of evidence for structure in lncRNAs. Nat. Methods 14, 45–48.
