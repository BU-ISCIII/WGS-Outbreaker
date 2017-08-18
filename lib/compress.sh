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
SAMPLE_NAMES=$3
FASTQ_PAIRED_R1=$4
FASTQ_PAIRED_R2=$5
FASTQ_UNPAIRED_R1=$6
FASTQ_UNPAIRED_R2=$7


if [ "$use_sge" = "1" ]; then
        sample_number=$SGE_TASK_ID
else
        sample_number=$8
fi

SAMPLE=$( echo $SAMPLE_NAMES | tr ":" "\n" | head -$sample_number | tail -1)
FASTQ_PA_R1_NAME=$( echo $FASTQ_PAIRED_R1 | tr ":" "\n" | head -$sample_number | tail -1)
FASTQ_PA_R2_NAME=$( echo $FASTQ_PAIRED_R2 | tr ":" "\n" | head -$sample_number | tail -1)
FASTQ_UNPA_R1_NAME=$( echo $FASTQ_UNPAIRED_R1 | tr ":" "\n" | head -$sample_number | tail -1)
FASTQ_UNPA_R2_NAME=$( echo $FASTQ_UNPAIRED_R2 | tr ":" "\n" | head -$sample_number | tail -1)

echo -e "compressing $SAMPLE \n"

gzip $INPUT_DIR/$SAMPLE/$FASTQ_PA_R1_NAME
gzip $INPUT_DIR/$SAMPLE/$FASTQ_PA_R2_NAME
gzip $INPUT_DIR/$SAMPLE/$FASTQ_UNPA_R1_NAME
gzip $INPUT_DIR/$SAMPLE/$FASTQ_UNPA_R2_NAME

echo -e "compressging $SAMPLE finished \n"

