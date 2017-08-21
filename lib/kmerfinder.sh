#!/bin/bash

# Test whether the script is being executed with sge or not.
if [ -z $SGE_TASK_ID ]; then 
	use_sge=0
else 
	use_sge=1
fi

# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status 
#set -e  
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u 

## Usage




set -x
echo `date`


# Variables
INPUT_DIR=$1
OUTPUT_DIR=$2
THREADS=$3
SAMPLE_NAMES=$4
KMERFINDER_PATH=$5
BACT_DB_PATH=$6
OUTPUT_CONCAT_NAMES=$7
OUTPUT_KMERFINDER_NAMES=$8 



if [ "$use_sge" = "1" ]; then
	sample_number=$SGE_TASK_ID
else
	sample_number=$9
fi

SAMPLE=$( echo $SAMPLE_NAMES | tr ":" "\n" | head -$sample_number | tail -1)
INPUT_NAME=$( echo $OUTPUT_CONCAT_NAMES | tr ":" "\n" | head -$sample_number | tail -1)
OUTPUT_NAME=$( echo $OUTPUT_KMERFINDER_NAMES | tr ":" "\n" | head -$sample_number | tail -1)


echo -e "Running kmerfinder for $SAMPLE \n"

python $KMERFINDER_PATH/findtemplate.py \
	-i $OUTPUT_DIR/$INPUT_NAME \
	-t $BACT_DB_PATH \
	-o $OUTPUT_DIR/$OUTPUT_NAME \
	-x ATGAC \
	-w

echo -e "kmerfinder for $SAMPLE finished \n"

echo -e "Remove concat file for $SAMPLE \n"

rm $OUTPUT_DIR/$INPUT_NAME

echo -e "Concat file for $SAMPLE deleted"
