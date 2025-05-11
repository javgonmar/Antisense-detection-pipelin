#! /bin/bash

WD=$1
INSDIR=$2
SAMPLE_PATH=$3
INDEX=$4
SAMPLE_DIRECTORY=$5
I=$6
NUMBER_SAMPLES=$7
EXP=$8
NUMBER_FASTA_CONDITION=$9
CONDITION=${10}
NUMBER_CONDITIONS=${11}
PATHS=${12}

cd $SAMPLE_DIRECTORY

fastqc sample_${I}.fq.gz
mv sample_${I}_fastqc.html sample_c${CONDITION}_${I}.html
mv sample_${I}_fastqc.zip sample_c${CONDITION}_${I}.zip


STAR --genomeDir $INDEX --readFilesIn $SAMPLE_PATH --readFilesCommand "gunzip -c" --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --alignIntronMax 1000 --outFileNamePrefix sample_$I 
samtools index sample_${I}Aligned.sortedByCoord.out.bam
bamCoverage -b sample_${I}Aligned.sortedByCoord.out.bam --filterRNAstrand reverse -o bamCoverage_${I}.bw
bedtools intersect -abam sample_${I}Aligned.sortedByCoord.out.bam -b ../../../annotation/annot.bed -s > sample_${I}_antisense_reads.bam
samtools index sample_${I}_antisense_reads.bam
bamCoverage -b sample_${I}_antisense_reads.bam --filterRNAstrand reverse -o antisense_reads_${I}.bw

stringtie -G ../../../annotation/annot.gtf -o sample_${I}.gtf sample_${I}_antisense_reads.bam
cufflinks -o antisense_assembly_${I} sample_${I}_antisense_reads.bam

gffread sample_${I}.gtf -T -o sample_${I}.gff3
gffread -E sample_${I}.gff3 -T -o sample_${I}.bed

cd antisense_assembly_${I}
gffread transcripts.gtf -T -o transcripts.gff3
gffread -E transcripts.gff3 -T -o transcripts.bed
echo "sample_${I}:" $(realpath transcripts.gtf) >> ../../../../scripts/transcripts_c${CONDITION}_count.txt

cd ..

cd ../../../scripts
touch condition_count.txt
cd
cd $SAMPLE_DIRECTORY

NUM_LINES_TRANSCRIPTS=$(wc -l < ../../../scripts/transcripts_c${CONDITION}_count.txt)
if [ $NUM_LINES_TRANSCRIPTS -eq $NUMBER_FASTA_CONDITION ]
then
	echo "Transcripts for condition $CONDITION are processed"
	echo "Intersecting for condition $CONDITION"
	echo "Condition $CONDITION is completed" >> ../../../scripts/condition_count.txt
	cd
	cd $INSDIR/$EXP
	sbatch --job-name=int_${CONDITION} --output=scripts/out_int_${CONDITION} --error=scripts/er_int_${CONDITION} ../intersect_c${CONDITION}.sh  scripts/transcripts_c${CONDITION}_count.txt $CONDITION $INSDIR $EXP $NUMBER_CONDITIONS $NUMBER_SAMPLES $PATHS
	cd
	cd $SAMPLE_DIRECTORY
else
	echo "Still processing for other conditions"
fi
