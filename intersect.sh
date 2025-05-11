#! /bin/bash

PARAMS=$1
CONDITION=$2
INSDIR=$3
EXP=$4
NUMBER_CONDITIONS=$5
NUMBER_SAMPLES=$6
PATHS=$7

NUM_FASTA=$(wc -l < $PARAMS)
echo "Number of fasta for condition $CONDITION is $NUM_FASTA"

cd
cd $INSDIR/$EXP

bedtools intersect -a $(head -n 1 "$PARAMS" | cut -d':' -f2 | sed 's/^ *//') $(tail -n +2 "$PARAMS" | cut -d':' -f2 | sed 's/^ *//; s/^/-b /') > results/transcripts_c${CONDITION}.gtf
echo "Transcripts intersected for condition $CONDITION done"
realpath "results/transcripts_c${CONDITION}.gtf" >> scripts/intersection_count.txt

COUNT_INTERSECT=$(wc -l < scripts/intersection_count.txt)
if [ $COUNT_INTERSECT -eq $NUMBER_CONDITIONS ]
then
	echo "Starting generation combined gtf file for all conditions"
	stringtie --merge -o results/annot_merged.gtf scripts/intersection_count.txt
	cd results
	Rscript $INSDIR/filter_antisense_script.R annot_merged.gtf ../annotation/annot.gtf
	while [ ! -f "annot_merged_filtered.gtf" ]
	do
		echo "Waiting for merging..."
		sleep 1
	done
	echo "Filtered gtf file was created"
	cat annot_merged_filtered.gtf ../annotation/annot.gtf >> annot_with_antisense.gtf
	cuffcompare -r ../annotation/annot.gtf -o cuffcompare annot_with_antisense.gtf
	gffcompare -r ../annotation/annot.gtf -o gffcompare annot_with_antisense.gtf
	paste cuffcompare.tracking gffcompare.tracking > fusion_compare.tracking
	cd ..

	for ((a=1; a <=${NUMBER_CONDITIONS}; a++));
	do
        	NUMBER_FASTA_CONDITION=$(grep number_fasta_condition_${a} $INSDIR/$PATHS | awk '{print($2)}')
        	for ((i=1; i <=${NUMBER_FASTA_CONDITION}; i++));
        	do
			echo "Launching expression scripts"
			sbatch --job-name=exp_c${a}_${i} --output=$INSDIR/$EXP/scripts/out_exp_c${a}_${i} --error=$INSDIR/$EXP/scripts/er_exp_c${a}_${i} $INSDIR/expression.sh $INSDIR/$EXP/samples/condition_${a}/sample_${i} $i $a $INSDIR $EXP $NUMBER_SAMPLES $NUMBER_CONDITIONS

        	done
	done

else
	echo "Still processing for other conditions"
fi
