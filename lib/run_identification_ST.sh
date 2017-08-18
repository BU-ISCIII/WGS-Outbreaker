#!/bin/bash

set -e
set -u
set -x

#Variables
USE_SGE=$1
KMERFINDER=$2
OUTPUT_DIR=$3
BACT_DB_PATH=$4
KMERFINDER_PATH=$5
THREADS=$6
FASTQ_R1_list=$7
FASTQ_R2_list=$8
sample_number=$9
SAMPLE_NAMES=${10}
OUTPUT_CONCAT_NAMES=${11}
OUTPUT_KMERFINDER_NAMES=${12}

##create directories
mkdir -p $OUTPUT_DIR/kmerfinder

jobid_trimmomatic=$( cat $OUTPUT_DIR/logs/jobids.txt | grep -w "TRIMMOMATIC" | cut -d ':' -f2 )

if [ -d $OUTPUT_DIR/QC/trimmomatic ]; then
	DIR=$OUTPUT_DIR/QC/trimmomatic
	CONCAT_ARGS="${SGE_ARGS} -hold_jid ${jobid_trimmomatic}"
else
	DIR=$OUTPUT_DIR/raw
fi


CONCAT_CMD="$SCRIPTS_DIR/concat.sh $DIR $OUTPUT_DIR/kmerfinder $SAMPLE_NAMES $FASTQ_R1_list $FASTQ_R2_list $OUTPUT_CONCAT_NAMES "
KMERFINDER_CMD="$SCRIPTS_DIR/kmerfinder.sh $OUTPUT_DIR/kmerfinder $OUTPUT_DIR/kmerfinder $THREADS $SAMPLE_NAMES $KMERFINDER_PATH $BACT_DB_PATH $OUTPUT_CONCAT_NAMES $OUTPUT_KMERFINDER_NAMES "

if [ $KMERFINDER == "YES" ]; then
	if [ "$USE_SGE" = "1" ]; then
		CONCATFILES=$( qsub $CONCAT_ARGS -N $JOBNAME.CONCATFILES $CONCAT_CMD)
	jobid_concat=$( echo $CONCATFILES | cut -d ' ' -f3 | cut -d '.' -f1 )
	echo -e "CONCAT FILES:$jobid_concat\n" >> $OUTPUT_DIR/logs/jobsids.txt 
			
		KMERFINDER_ARG="{$SGE_ARGS -pe openmp $THREADS -hold_jib $jobid_concat}"
		KMERFINDER=$( qsub $KMERFINDER_ARG -t1-$sample_number -N $JOBNAME.KMERFINDER $KMERFINDER_CMD)
	jobid_kmerfinder=$( echo $KMERFINDER | cut -d ' ' -f3 | cut -d '.' -f1 )
	echo -e "KMERFINDER:$jobif_kmerfinder\n" >> $OUTPUT_DIR/logs/jobsids.txt
	else
		for count in  `seq 1 $sample_number`; do
			echo "Running concat files on sample $count"
			CONCATFILES=$($CONCAT_CMD $count)
			echo "Running kmerfinder on sample $count"
			KMERFINDER=$($KMERFINDER_CMD $count)
		done
	fi	
fi



