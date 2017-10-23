#!/bin/bash
## Author: A.Hernandez
## Usage: compress.sh


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

# VARIABLES

input_dir=$1
output_dir=$2
samples=$3
trimmedFastqArray_paired_R1_list=$4
trimmedFastqArray_paired_R2_list=$5
trimmedFastqArray_unpaired_R1_list=$6
trimmedFastqArray_unpaired_R2_list=$7


if [ "$use_sge" = "1" ]; then
        sample_count=$SGE_TASK_ID
else
        sample_count=$8
fi

sample=$( echo $samples | tr ":" "\n" | head -$sample_count | tail -1)
trimmedFastqArray_paired_R1=$( echo $trimmedFastqArray_paired_R1_list | tr ":" "\n" | head -$sample_count | tail -1)
trimmedFastqArray_paired_R2=$( echo $trimmedFastqArray_paired_R2_list | tr ":" "\n" | head -$sample_count | tail -1)
trimmedFastqArray_unpaired_R1=$( echo $trimmedFastqArray_unpaired_R1_list | tr ":" "\n" | head -$sample_count | tail -1)
trimmedFastqArray_unpaired_R2=$( echo $trimmedFastqArray_unpaired_R2_list | tr ":" "\n" | head -$sample_count | tail -1)

echo -e "compressing $sample \n"

gzip $input_dir/$sample/$trimmedFastqArray_paired_R1
gzip $input_dir/$sample/$trimmedFastqArray_paired_R2
gzip $input_dir/$sample/$trimmedFastqArray_unpaired_R1
gzip $input_dir/$sample/$trimmedFastqArray_unpaired_R2

echo -e "compressging $sample finished \n"

