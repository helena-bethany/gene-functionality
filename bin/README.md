# Script usage

* [Sequence processing scripts](#sequence-processing-scripts)

* [Feature generation](#feature-generation)

* [Random forest](#random-forest)

-----------------------------------------------------------------------

## Sequence processing scripts

There are three scripts that are used to process the functional sequences and generate the corresponding negative control sequences. Short ncRNAs are processed using [short-ncrna-processing.sh](./bin/short-ncrna-processing.sh), lncRNAs are processed using [long-ncrna-processing.sh](./bin/long-ncrna-processing.sh) and protein-coding genes are processed using [protein-coding-processing.sh](./bin/protein-coding-processing.sh).

These scripts have six required input:
  1. The location of the folder which contains all generated local databases and downloaded functional sequence files, as described [previously](./README.md).
  2. Name of the downloaded functional sequence file. Eg: ```200625-hgnc-PC-refseq-ids.txt``` OR ```RNAcentral-short-ncrna.fa```
  3. This input is different depending on whether the sequences analysed are short/long ncRNAs or protein-coding genes.
  * Short/Long ncRNAs: the name of RNAcentral chromosome coordinates file. Eg: ```Homo_sapiens.GRCh38.bed```
  * Protein-coding genes: the name of the human genome GFF. Eg: ```GCF_000001405.39_GRCh38.p13_genomic.gff```
  4. Name of the converted human genome CSV file. Eg: ```GRCh38.p13_genome.csv```
  5. Name of the uniprot bed file. Eg: ```uni-prot-ucsc.bed```
  6. Name of the GENCODE bed file. Eg: ```gencode-ncrna-annotation.bed```

```
# Example of short/long-ncrna-processing.sh usage
bash short-ncrna-processing.sh /location/to/local-databases RNAcentral-short-ncrna.fa Homo_sapiens.GRCh38.bed GRCh38.p13_genome.csv uni-prot-ucsc.bed gencode-ncrna-annotation.bed
bash long-ncrna-processing.sh /location/to/local-databases RNAcentral-lncrna.fa Homo_sapiens.GRCh38.bed GRCh38.p13_genome.csv uni-prot-ucsc.bed gencode-ncrna-annotation.bed

# Example of protein-coding-processing.sh usage
bash protein-coding-processing.sh /location/to/local-databases 200625-hgnc-PC-refseq-ids.txt GCF_000001405.39_GRCh38.p13_genomic.gff GRCh38.p13_genome.csv uni-prot-ucsc.bed gencode-ncrna-annotation.bed
```

The outputs from these scripts are the CSV and FASTA files for the functional sequences and generated negative control sequences, with examples of the naming conventions listed below.
* Functional ncRNA sequences: ```200702-functional-ncrna-dataset.csv``` and ```200702-functional-ncrna-seq.fa```
* Functional lncRNA exons: ```200702-functional-ncrna-exon2-dataset.csv``` and ```200702-functional-ncrna-exon2-seq.fa```
* Functional protein-coding exons: ```200702-protein-exon2-dataset.csv``` and ```200702-protein-exon2-seq.fa```
* Negative control sequences: ```200702-negative-control-dataset.csv``` and ```200702-negative-control-seq.fa```

-----------------------------------------------------------------------------------------------------------------------

## Feature generation

The script [feature-generation.sh](./bin/feature-generation.sh) calculates all chosen features of functions for protein-coding exons, lncRNA exons, short ncRNAs and negative-control regions. 

This has two required input:
 1. The file identifer for the sequence dataset and fasta file. Eg: The files ```200702-functional-ncrna-dataset.csv``` and ```200702-functional-ncrna-seq.fa``` have the file identifier ```200702-functional-ncrna```.
 2. The location of the folder which contains all generated local databases, as described [previously](./README.md).

If downloaded or reformatted files have been renamed differently to what was described [previously](./README.md), then the script will prompt for the custom name of each file. The script will exit prematurely if a required dependency is not available.

```
# Example of script usage
bash feature-generation 200702-functional-ncrna /location/to/local-databases

# Example of additional prompt
Please enter custom name of RNA:RNA interaction database (fa):
```

The main output for this script is a CSV file that contains all the calculated feature data, with the file name being based on the initial dataset and the date that the script finished running. Eg: ```200703-functional-ncrna-final-dataset.csv``` 

-----------------------------------------------------------------------

## Random forest

The script [random-forest-analysis.R](./bin/random-forest-analysis.R) generates the random forest models and outputs the performance metrics calculated across all models. Prior to running, a combined file should be generated which contains the final dataset files from [feature-generation.sh](#feature-generation) for the functional sequences and corresponding negative control sequences. This combined file should also have the Sequence ID, Chromosome, End position and Sequence columns removed, so only feature data is included for the random forest analysis. 

```
# Generate combined features dataset:
# Pipe two removes the duplicate header row
# Pipe three removes the Sequence ID (1), Chromosome (3), End Position (5) and Sequence (6) columns
cat 200703-functional-ncrna-final-dataset.csv 200703-negative-control-final-dataset.csv | awk '!a[$0]++' | cut -d ',' -f 2,4,7- > 200703-short-ncrna-final-dataset-RF.csv

# Load R and run random forest analysis
source("random-forest-analysis.R")
randomforestx("200703-short-ncrna-final-dataset-RF.csv")
```

There are four files generated as output, which use the date the analysis was function was run and input CSV file to generate the output file names.
* The "stats" file is the minimum, average and maximum performance metrics across all 100 random forest models. 
  * Eg: ```2020-07-03-stats-short-ncrna-short-ncrna-final-dataset-RF.csv```
* The "importance" file is the variable importances recorded for each model, in addition to the average. 
  * Eg: ```2020-07-03-importance-short-ncrna-short-ncrna-final-dataset-RF.csv```
* The "error" file contains the OOB values for each model, in addition to the average. 
  * Eg: ```2020-07-03-error-short-ncrna-short-ncrna-final-dataset-RF.csv```
* The "predictions" file contains the sum of all confusion matrices generated by the models. 
  * Eg: ```2020-07-03-predictions-short-ncrna-short-ncrna-final-dataset-RF.csv```

