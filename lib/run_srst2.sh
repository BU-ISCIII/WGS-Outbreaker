#!/bin/bash

#usage




# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x

USE_SGE=$1
SRST2=$2
OUTPUT_DIR=$3
SRST2_DB_PATH=$4
THREADS=$5
FASTQ_R1_LIST=$6
FASTQ_R2_LIST=$7
sample_number=$8
sample_names=$9
resistance_list=${10}
plasmid_list=${11}
mlst_list=${12}

#create directories
mkdir -p $OUTPUT_DIR/srst2

jobid_trimmomatic=$( cat $OUTPUT_DIR/logs/jobids.txt | grep -w "TRIMMOMATIC" | cut -d ':' -f2 )

if [ -d $OUTPUT_DIR/QC/trimmomatic ]; then
        DIR=$OUTPUT_DIR/QC/trimmomatic
        SRST2_ARGS="${SGE_ARGS} -hold_jid ${jobid_trimmomatic}"
else
        DIR=$OUTPUT_DIR/raw
fi

RESISTANCE_CMD="$SCRIPTS_DIR/resistance.sh $DIR $OUTPUT_DIR/srst2 $SRST2_DB_PATH $sample_names $FASTQ_R1_LIST $FASTQ_R2_LIST $resistance_list"
#PLASMID_CMD="$SCRIPTS_DIR/plasmid.sh $DIR $OUTPUT_DIR/srst2 $SRST2_DB_PATH $sample_names $FASTQ_R1_list $FASTQ_R2_list $plasmid_list"
#MLST_CMD="$SCRIPTS_DIR/mlst.sh $DIR $OUTPUT_DIR/srst2 $SRST2_DB_PATH $sample_names $FASTQ_R1_list $FASTQ_R2_list $mlst_list"


if [ $SRST2 == "YES" ]:then
	if [ "$USE_SG" = "1" ]; then
		RESISTANCE=$( qsub $SRST2_ARG -N $JOBNAME.RESISTANCE $RESISTANCE_CMD )
		jobid_resistance=$( echo $RESISTANCE | cut -d ' ' -f3 | cut -d '.' -f1 )
		echo -e "RESISTANCE FILES:$jobid_resistance\n" >> $OUTPUT_DIR/logs/jobids.txt

	#	PLASMIDS=$( qsub $SRST2_ARG -N $JOBNAME.PLASMIDS $PLASMID_CMD )
        #       jobid_plasmid=$( echo $PLASMIDS | cut -d ' ' -f3 | cut -d '.' -f1 )
        #        echo -e "PLASMID FILES:$jobid_plasmid\n" >> $OUTPUT_DIR/logs/jobids.txt
	
	#	MLST=$( qsub $SRST2_ARG -N $JOBNAME.MLST $MLST_CMD )
        #       jobid_mlst=$( echo $MLST | cut -d ' ' -f3 | cut -d '.' -f1 )
        #        echo -e "MLST FILES:$jobid_mlst\n" >> $OUTPUT_DIR/logs/jobids.txt
		
	else
		for count in `seq 1 $sample_number`;do
		echo "Running resistance on sample $count"
		RESISTANCE=$($RESISTANCE_CMD $count)
	#	echo "Running plasmid on sample $count"
        #       PLASMIDS=$($PLASMID_CMD $count)
	#	echo "Running mlst on sample $count"
        #       MLST=$($MLST_CMD $count)
		done
	fi
fi

