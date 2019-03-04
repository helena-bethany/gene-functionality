#Extracting expression for a region across 4 ENCODE datasets

grep -v "Start" $1 | cut -d ',' -f 3,4,5 | tr ',' ' ' | perl -lane '{print "chr$F[0]:$F[1]-$F[2]"}' > locations

IFS=$'\n'
calc() { awk "BEGIN{print $*}"; }

echo HepG2,hESC,K562,GM12878,AverageLINE,AveragebpLINE > $2line.csv

for line in $(cat locations)
do
	HepG2=$(samtools view ENCODE/ENCFF067CVP.bam $line | wc -l)
	hESC=$(samtools view ENCODE/ENCFF089EWC.bam $line | wc -l)
	K562=$(samtools view ENCODE/ENCFF796BVP.bam $line | wc -l)
	GM12878=$(samtools view ENCODE/ENCFF893HSY.bam $line | wc -l)
	total=$(calc $HepG2+$hESC+$K562+$GM12878)
	average=$(calc $total/4)
	echo $line > line.txt
	start=$(cat line.txt | tr ":" "," | cut -d ',' -f 2 | tr '-' ' ' | perl -lane '{print "$F[0]"}')
	end=$(cat line.txt | tr ":" "," | cut -d ',' -f 2 | tr '-' ' ' | perl -lane '{print "$F[1]"}')
	length=$(calc $end-$start)
	bp=$(calc $average/$length)
	echo $HepG2, $hESC, $K562, $GM12878, $average, $bp > info
	cat info >> $2line.csv
done

echo SmoothMuscleCell,Hepatocyte,NeuralProgenitorCell,Myocyte,BipolarNeuron,Myotube,HematopoieticMultipotent,CardiacMuscle,AverageTYPE,AveragebpTYPE > $2type.csv

for line in $(cat locations)
do
	smoothmuscle=$(samtools view ENCODE/ENCFF369DYD.bam $line | wc -l)
	hepatocyte=$(samtools view ENCODE/ENCFF907AIK.bam $line | wc -l)
	progenitor=$(samtools view ENCODE/ENCFF766WSM.bam $line | wc -l)
	myocyte=$(samtools view ENCODE/ENCFF722EAR.bam $line | wc -l)
	neuron=$(samtools view ENCODE/ENCFF713UNS.bam $line | wc -l)
	myotube=$(samtools view ENCODE/ENCFF178TTA.bam $line | wc -l)
	multipotent=$(samtools view ENCODE/ENCFF065MVD.bam $line | wc -l)
	cardiac=$(samtools view ENCODE/ENCFF475WLJ.bam $line | wc -l)
	total=$(calc $smoothmuscle+$hepatocyte+$progenitor+$myocyte+$neuron+$myotube+$multipotent+$cardiac)
	average=$(calc $total/8)
	echo $line > line.txt
        start=$(cat line.txt | tr ":" "," | cut -d ',' -f 2 | tr '-' ' ' | perl -lane '{print "$F[0]"}')
        end=$(cat line.txt | tr ":" "," | cut -d ',' -f 2 | tr '-' ' ' | perl -lane '{print "$F[1]"}')
        length=$(calc $end-$start)
        bp=$(calc $average/$length)
	echo $smoothmuscle, $hepatocyte, $progenitor, $myocyte, $neuron, $myotube, $multipotent, $cardiac, $average, $bp > info
	cat info >> $2type.csv
done

paste -d ',' $1 $2line.csv $2type.csv > $2.csv
rm -rf line.txt
