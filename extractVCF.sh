grep -v 'Chromosome' $1 | cut -d ',' -f 3,4,5 | tr ',' ' ' | perl -lane '{print "$F[0]:$F[1]-$F[2]"}' > coordinates

IFS=$'\n'
var=1
rm -rf \FASTA
mkdir FASTA

for line in $(cat coordinates)
do
        name='RNA'$var
	echo $line | cut -d ':' -f 1 > text.file
	chr=$(cat text.file)
	chr="${chr%\"}"
	chr="${chr#\"}"
	tabix -f -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz $line > $name.vcf
	mv $name.vcf \FASTA
	var=$((var+1))
	rm -rf text.file
done

