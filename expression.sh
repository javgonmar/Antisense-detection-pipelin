#! /bin/bash

SAMPLE_PATH=$1
SAMPLE=$2
CONDITION=$3
INSDIR=$4
EXP=$5
NUMBER_SAMPLES=$6
NUMBER_CONDITIONS=$7

cd $INSDIR

for ((i=1; i <=${NUMBER_CONDITIONS}; i++));
do
	rm intersect_c${i}.sh
	rm sample_proc_condition${i}.sh
done


cd $SAMPLE_PATH
stringtie -G ../../../results/annot_with_antisense.gtf -o sample_${SAMPLE}.gtf sample_${SAMPLE}_antisense_reads.bam
stringtie -e -G ../../../results/annot_with_antisense.gtf -o sample_c${CONDITION}_${SAMPLE}.gtf -A sample_c${CONDITION}_${SAMPLE}.tsv sample_${SAMPLE}_antisense_reads.bam

featureCounts -s 0 -g gene_id -a ../../../results/annot_with_antisense.gtf -o sample_c${CONDITION}_${SAMPLE}.txt sample_${SAMPLE}_antisense_reads.bam
realpath sample_c${CONDITION}_${SAMPLE}.txt >> $INSDIR/$EXP/scripts/expression_count.txt

NUM_LINES_EXPRESSION=$(wc -l < $INSDIR/$EXP/scripts/expression_count.txt)
if [ $NUM_LINES_EXPRESSION -eq $NUMBER_SAMPLES ]
then
	echo "Expression for all files is done"
else
	echo "Still processing for other files"
fi
