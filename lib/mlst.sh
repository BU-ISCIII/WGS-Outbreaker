#!/bin/bash

#usage




# Test whether the script is being executed with sge or not.

if [ -z $SGE_TASK_ID ]; then
        use_sge=0
else
        use_sge=1
fi

# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x

INPUT_DIR=$1
OUTPUT_DIR=$2
SRST2_DB_PATH_mlst=$3
SRST2_DB_PATH_definitions=$4
SAMPLE_NAMES=$5
FASTQ_R1_LIST=$6
FASTQ_R2_LIST=$7
mlst_list=$8
SRST2_DELIMITER=$9

if [ "$use_sge" = "1" ]; then
        sample_number=$SGE_TASK_ID
else
        sample_number=${10}
fi

SAMPLE=$( echo $SAMPLE_NAMES | tr ":" "\n" | head -$sample_number | tail -1)
FASTQ_R1_NAME=$( echo $FASTQ_R1_LIST | tr ":" "\n" | head -$sample_number | tail -1)
FASTQ_R2_NAME=$( echo $FASTQ_R2_LIST | tr ":" "\n" | head -$sample_number | tail -1)
OUTPUT_NAME=$( echo $mlst_list | tr ":" "\n" | head -$sample_number | tail -1)

echo -e "Running mlst.sh for $SAMPLE"

srst2 --input_pe $INPUT_DIR/$SAMPLE/$FASTQ_R1_NAME $INPUT_DIR/$SAMPLE/$FASTQ_R2_NAME --forward trimmed_R1 --reverse trimmed_R2 --log --output $OUTPUT_DIR/$OUTPUT_NAME --mlst_db $SRST2_DB_PATH_mlst --mlst_definitions $SRST2_DB_PATH_definitions --mlst_delimiter $SRST2_DELIMITER

echo -e "mlst.sf for $SAMPLE finished"



