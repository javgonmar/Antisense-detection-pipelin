#! /bin/bash

PARAMS=$1

if [ $# -gt 1 ]
then
        echo "The number of arguments is: $#"
        echo "Usage: antipipe.sh <params.file>"
        echo "An example of params.file can be found in params.txt"
        exit
elif [ $# -lt 1 ]
then
        echo "===================================="
        echo "You need to add a params file (.txt)"
        echo "===================================="
        exit
fi

echo "============================"
echo "   Extracting parameters    "
echo "============================"

WD=$(grep working_directory $PARAMS | awk '{print($2)}')
echo "Working directory path is $WD"
INSDIR=$(grep installation_directory $PARAMS | awk '{print($2)}')
echo "Installation directory is $INSDIR"
EXP=$(grep experiment_name $PARAMS | awk '{print($2)}')
echo "Experiment name is $EXP"

GENOME=$(grep genome_path $PARAMS | awk '{print($2)}')
echo "Genome path is $GENOME"
ANNOT=$(grep annotation_path $PARAMS | awk '{print($2)}')
echo "Annotation path is $ANNOT"

NUMBER_CONDITIONS=$(grep number_conditions $PARAMS | awk '{print($2)}')
echo "Number of conditions is $NUMBER_CONDITIONS"
NUMBER_SAMPLES=$(grep number_samples $PARAMS | awk '{print($2)}')
echo "The number of samples is $NUMBER_SAMPLES"


echo "==========================="
echo "    Creating workspace"
echo "==========================="

cd $WD

mkdir $EXP
cd $EXP
mkdir samples genome annotation results scripts
cd samples

TOTAL_FASTA_FILES=0
for ((i=1; i <=${NUMBER_CONDITIONS}; i++));
do
	cd ../..
	eval NUMBER_FASTA_CONDITION_${i}=$(grep number_fasta_condition_${i} $INSDIR/$PARAMS | awk '{print($2)}')
	eval echo "Number fasta files for condition $i is \$NUMBER_FASTA_CONDITION_${i}"

	TOTAL_FASTA_FILES=$((TOTAL_FASTA_FILES + NUMBER_FASTA_CONDITION_${i}))
	echo "Number of total fasta in iteration $i is $TOTAL_FASTA_FILES"

	cd $EXP/samples
	mkdir condition_${i}
	cd
	cd $INSDIR
	cp sample_proc.sh sample_proc_condition${i}.sh
	cp intersect.sh intersect_c${i}.sh
	cd $EXP/samples
done

if [[ "$TOTAL_FASTA_FILES" -ne "$NUMBER_SAMPLES" ]]
then
	echo "Error: The number of fasta files does NOT fit with the number of samples expected"
	echo "Please revise the params file given"
	echo "Exiting"
	exit 1
fi

cd ..

cp $GENOME genome/genome.fa
cp $ANNOT annotation/annot.gtf

cd annotation
gffread annot.gtf -T -o annot.gff3
gffread -E annot.gff3 -T -o annot.bed
sortBed -i annot.bed > annot_sorted.bed
cd ..

echo "============================"
echo "   Creating genome index"
echo "============================"

GEN_PARAM=$(grep genome_parameter $INSDIR/$PARAMS | awk '{print($2)}')
echo "Genome length parameter is $GEN_PARAM"

cd genome
STAR --runMode genomeGenerate --genomeDir index --genomeFastaFiles genome.fa --sjdbGTFfile ../annotation/annot.gtf --genomeSAindexNbases $GEN_PARAM
cd ..

for ((a=1; a <=${NUMBER_CONDITIONS}; a++));
do
	NUMBER_FASTA_CONDITION=$(grep number_fasta_condition_${a} $INSDIR/$PARAMS | awk '{print($2)}')
	for ((i=1; i <=${NUMBER_FASTA_CONDITION}; i++));
	do
		echo "Sample $i for condition $a"
		SAMPLE_PATH=$(grep fasta_path_condition_${a}_${i} $INSDIR/$PARAMS | awk '{print($2)}') 
		cd samples/condition_${a}
		mkdir -p  sample_${i}
		cp ${SAMPLE_PATH} $WD/$EXP/samples/condition_${a}/sample_${i}/sample_${i}.fq.gz
		cd ../../..
		sbatch --job-name=c${a}_${i} --output=$EXP/scripts/out_c${a}_${i} --error=$EXP/scripts/er_c${a}_${i} sample_proc_condition${a}.sh $WD $INSDIR $SAMPLE_PATH $WD/$EXP/genome/index $WD/$EXP/samples/condition_${a}/sample_${i} ${i} $NUMBER_SAMPLES $EXP $NUMBER_FASTA_CONDITION $a $NUMBER_CONDITIONS $PARAMS
		cd $EXP
	done
done

