# Extracting expression for a region across 12 ENCODE datasets

#####################################################################

# $1 is the input ncRNA dataset file

IFS=$'\n'
d=$(date +%y%m%d) # date identifier
id=$( echo $1 | cut -d '.' -f 1 | cut -d '-' -f 4 )
count=$(( $id + 1 ))
calc() { awk "BEGIN{print $*}"; }

#####################################################################

# Create indexed bam files

read -p "Make indexed BAM files? (y/n)" response

if [[ $response == "y" ]]
then
    
    samtools index ENCFF067CVP.bam
    samtools index ENCFF089EWC.bam
    samtools index ENCFF796BVP.bam
    samtools index ENCFF893HSY.bam
    
    samtools index ENCFF369DYD.bam
    samtools index ENCFF907AIK.bam
    samtools index ENCFF766WSM.bam 
    samtools index ENCFF722EAR.bam
    samtools index ENCFF713UNS.bam
    samtools index ENCFF178TTA.bam
    samtools index ENCFF065MVD.bam 
    samtools index ENCFF475WLJ.bam 
	
else
    :
fi

#####################################################################

# Cell line datasets

grep -v "Start" $1 | cut -d ',' -f 3,4,5 | tr ',' ' ' | perl -lane '{print "chr$F[0]:$F[1]-$F[2]"}' > locations

echo HepG2,hESC,K562,GM12878,AverageLINE,AveragebpLINE > $d-line.csv

for line in $(cat locations)
do
	HepG2=$(samtools view ENCFF067CVP.bam $line | wc -l)
	hESC=$(samtools view ENCFF089EWC.bam $line | wc -l)
	K562=$(samtools view ENCFF796BVP.bam $line | wc -l)
	GM12878=$(samtools view ENCFF893HSY.bam $line | wc -l)
	total=$(calc $HepG2+$hESC+$K562+$GM12878)
	average=$(calc $total/4)
	echo $line > line.txt
	start=$(cat line.txt | tr ":" "," | cut -d ',' -f 2 | tr '-' ' ' | perl -lane '{print "$F[0]"}')
	end=$(cat line.txt | tr ":" "," | cut -d ',' -f 2 | tr '-' ' ' | perl -lane '{print "$F[1]"}')
	length=$(calc $end-$start)
	bp=$(calc $average/$length)
	echo $HepG2, $hESC, $K562, $GM12878, $average, $bp > info
	cat info >> $d-line.csv
done

#####################################################################

# Differentiated cell type datasets

echo SmoothMuscleCell,Hepatocyte,NeuralProgenitorCell,Myocyte,BipolarNeuron,Myotube,HematopoieticMultipotent,CardiacMuscle,AverageTYPE,AveragebpTYPE > $d-type.csv

for line in $(cat locations)
do
	smoothmuscle=$(samtools view ENCFF369DYD.bam $line | wc -l)
	hepatocyte=$(samtools view ENCFF907AIK.bam $line | wc -l)
	progenitor=$(samtools view ENCFF766WSM.bam $line | wc -l)
	myocyte=$(samtools view ENCFF722EAR.bam $line | wc -l)
	neuron=$(samtools view ENCFF713UNS.bam $line | wc -l)
	myotube=$(samtools view ENCFF178TTA.bam $line | wc -l)
	multipotent=$(samtools view ENCFF065MVD.bam $line | wc -l)
	cardiac=$(samtools view ENCFF475WLJ.bam $line | wc -l)
	total=$(calc $smoothmuscle+$hepatocyte+$progenitor+$myocyte+$neuron+$myotube+$multipotent+$cardiac)
	average=$(calc $total/8)
	echo $line > line.txt
        start=$(cat line.txt | tr ":" "," | cut -d ',' -f 2 | tr '-' ' ' | perl -lane '{print "$F[0]"}')
        end=$(cat line.txt | tr ":" "," | cut -d ',' -f 2 | tr '-' ' ' | perl -lane '{print "$F[1]"}')
        length=$(calc $end-$start)
        bp=$(calc $average/$length)
	echo $smoothmuscle, $hepatocyte, $progenitor, $myocyte, $neuron, $myotube, $multipotent, $cardiac, $average, $bp > info
	cat info >> $d-type.csv
done

paste -d ',' $1 $d-line.csv $d-type.csv > $d-ncrna-dataset-$count.csv

#####################################################################

rm -rf line.txt
rm -rf locations
rm -rf info
