#!/bin/bash
##Author:A.Hernandez
## Usage :kmerfinder.sh ...


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
#Print commands and their arguments as they are executed.
set -x

echo `date`


# VARIABLES

threads=$1
dir=$2
output_dir=$3
samples=$4
concatFastq_list=$5
kmerfinderST_list=$6
kmerfinder_path=$7
bact_db_path=$8


if [ "$use_sge" = "1" ]; then
	sample_count=$sge_task_id
else
	sample_count=$9
fi

sample=$( echo $samples | tr ":" "\n" | head -$sample_count | tail -1)
concatFastq=$( echo $concatFastq_list | tr ":" "\n" | head -$sample_count | tail -1)
kmerfinderST=$( echo $kmerfinderST_list | tr ":" "\n" | head -$sample_count | tail -1)


echo -e "Running kmerfinder for $sample \n"

python $kmerfinder_path/findtemplate.py \
	-i $dir/$concatFastq \
	-t $bact_db_path \
	-o $output_dir/$kmerfinderST \
	-x ATGAC \
	-w

echo -e "kmerfinder for $sample finished \n"

echo -e "Remove concat file for $sample \n"

rm $output_dir/$concatFastq

echo -e "Concat file for $sample deleted"
