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
SRST2_DB_PATH_ARGannot=$3
sample_names=$4
FASTQ_compress_R1_list=$5
FASTQ_compress_R2_list=$6
OUTPUT_FILE=$7


if [ "$use_sge" = "1" ]; then
        sample_number=$SGE_TASK_ID
else
        sample_number=$8
fi

SAMPLE=$( echo $sample_names | tr ":" "\n" | head -$sample_number | tail -1)
FASTQ_R1_COMPRESS_NAME=$( echo $FASTQ_compress_R1_list | tr ":" "\n" | head -$sample_number | tail -1)
FASTQ_R2_COMPRESS_NAME=$( echo $FASTQ_compress_R2_list | tr ":" "\n" | head -$sample_number | tail -1)
OUTPUT_NAME=$( echo $OUTPUT_FILE | tr ":" "\n" | head -$sample_number | tail -1)

echo -e "Running resistance.sh for $SAMPLE \n\n"

srst2 --input_pe $INPUT_DIR/$SAMPLE/$FASTQ_R1_COMPRESS_NAME $INPUT_DIR/$SAMPLE/$FASTQ_R2_COMPRESS_NAME --forward trimmed_R1 --reverse trimmed_R2 --log --output $OUTPUT_DIR/$SAMPLE/$OUTPUT_NAME --gene_db $SRST2_DB_PATH_ARGannot

echo -e "Resistance.sh for $SAMPLE finished \n\n"
