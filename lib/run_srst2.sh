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
SRST2_DB_PATH_ARGannot=$4
SRST2_DB_PATH_PlasmidFinder=$5
SRST2_DB_PATH_mlst=$6
SRST2_DB_PATH_mlst_definitions=$7
THREADS=$8
FASTQ_compress_R1_list=$9
FASTQ_compress_R2_list=${10}
sample_number=${11}
sample_names=${12}
resistance_list=${13}
plasmid_list=${14}
mlst_list=${15}
#SRST2_DELIMITER=${16}
SUMMARY_RESULTS=${16}

SAMPLE=$( echo $sample_names | tr ":" "\n" | head -$sample_number | tail -1)


#create directories
mkdir -p $OUTPUT_DIR/srst2/$SAMPLE/summary


jobid_trimmomatic=$( cat $OUTPUT_DIR/logs/jobids.txt | grep -w "TRIMMOMATIC" | cut -d ':' -f2 )

if [ -d $OUTPUT_DIR/QC/trimmomatic ]; then
        DIR=$OUTPUT_DIR/QC/trimmomatic
        SRST2_ARGS="${SGE_ARGS} -hold_jid ${jobid_trimmomatic}"
else
        DIR=$OUTPUT_DIR/raw
fi

RESISTANCE_CMD="$SCRIPTS_DIR/resistance.sh $DIR $OUTPUT_DIR/srst2/$SAMPLE $SRST2_DB_PATH_ARGannot $sample_names $FASTQ_compress_R1_list $FASTQ_compress_R2_list $resistance_list"
PLASMID_CMD="$SCRIPTS_DIR/plasmid.sh $DIR $OUTPUT_DIR/srst2/$SAMPLE $SRST2_DB_PATH_PlasmidFinder $sample_names $FASTQ_compress_R1_list $FASTQ_compress_R2_list $plasmid_list"
MLST_CMD="$SCRIPTS_DIR/mlst.sh $DIR $OUTPUT_DIR/srst2/$SAMPLE $SRST2_DB_PATH_mlst $SRST2_DB_PATH_mlst_definitions $sample_names $FASTQ_compress_R1_list $FASTQ_compress_R2_list $mlst_list"
SUMMARY_CMD="$SCRIPTS_DIR/summary_srst2.sh $OUTPUT_DIR/srst2/$SAMPLE $OUTPUT_DIR/srst2/summary $SUMMARY_RESULTS $sample_names" 

if [ $SRST2 == "YES" ];then
	if [ "$USE_SGE" = "1" ]; then
		RESISTANCE=$( qsub $SRST2_ARG -N $JOBNAME.RESISTANCE $RESISTANCE_CMD )
		jobid_resistance=$( echo $RESISTANCE | cut -d ' ' -f3 | cut -d '.' -f1 )
		echo -e "RESISTANCE FILES:$jobid_resistance\n" >> $OUTPUT_DIR/logs/jobids.txt

		PLASMIDS=$( qsub $SRST2_ARG -N $JOBNAME.PLASMIDS $PLASMID_CMD )
                jobid_plasmid=$( echo $PLASMIDS | cut -d ' ' -f3 | cut -d '.' -f1 )
                echo -e "PLASMID FILES:$jobid_plasmid\n" >> $OUTPUT_DIR/logs/jobids.txt
	
		MLST=$( qsub $SRST2_ARG -N $JOBNAME.MLST $MLST_CMD )
                jobid_mlst=$( echo $MLST | cut -d ' ' -f3 | cut -d '.' -f1 )
                echo -e "MLST FILES:$jobid_mlst\n" >> $OUTPUT_DIR/logs/jobids.txt
		
		SUMMARY=$( qsub $SRST2_ARG -N $JOBNAME.SUMMARY $SUMMARY_CMD )
                jobid_mlst=$( echo $SUMMARY | cut -d ' ' -f3 | cut -d '.' -f1 )
                echo -e "SUMAMRY FILES:$jobid_summary\n" >> $OUTPUT_DIR/logs/jobids.txt

		
	else
		for count in `seq 1 $sample_number`; do
			echo "Running resistance on sample $count"
			RESISTANCE=$($RESISTANCE_CMD $count)
			echo "Running plasmid on sample $count"
               		PLASMIDS=$($PLASMID_CMD $count)
			echo "Running mlst on sample $count"
        		MLST=$($MLST_CMD $count)
			echo "Running summary on sample $count"
			SUMMARY=$($SUMMARY_CMD $count)
		done
	fi
fi
