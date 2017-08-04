#!/bin/bash

#usage: concat.sh

#Test whether the script is being used with sge or not.

if [-z $SGE_TASK_ID]; then
	use_sge=0

else
	use_sge=1

fi

# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status 
set -e  
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u 

## Usage



#Print a trace of simple commands, for commands, case commands, select commands, and arithmetic for commands and their arguments or associated word lists after they are expanded and before they are executed
set -x
echo `date`

# Variables
INPUT_DIR=$1
OUTPUT_DIR=$2
SAMPLE_NAMES=$3
FASTQ_FILES_R1=$4
FASTQ_FILES_R2=$5
OUTPUT_CONCAT_NAMES=$6

if [" use_sge" = "1"]; then
	sample_number=$SGE_TASK_ID
else
 	sample_number=$7
fi

SAMPLE=$( echo $SAMPLE_NAMES | tr ":" "\n" | head -$sample_number | tail -1)
FASTQ_R1=$( echo $FASTQ_FILES_R1 | tr ":" "\n" | head -$sample_number | tail -1)
FASTQ_R2=$( echo $FASTQ_FILES_R2 | tr ":" "\n" | head -$sample_number | tail -1)

echo -e "concat files for $SAMPLE \n\n"

zcat $INPUT_DIR/$SAMPLE/FASTQ_R* -o $OUTPUT_CONCAT_NAMES

echo -e "Concat $SAMPLE finished \n\n"
