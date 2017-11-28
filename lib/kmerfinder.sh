#!/bin/bash
## Author:A.Hernandez
## version v2.0.

if [ $# -eq 0 ];then
        echo -e "\nScript to run kmerfinder\n"
        echo -e "Usage: kmerfinder.sh input_dir output_dir samples_list concat_files kmerfinderST_list kmerfinder_path bacterial_db_path"
        exit
fi

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

# VARIABLES

dir=$1
output_dir=$2
samples=$3
concatFastq_list=$4
kmerfinderST_list=$5
kmerfinder_path=$6
bact_db_path=$7

if [ "$use_sge" = "1" ]; then
	sample_count=$SGE_TASK_ID
else
	sample_count=$8
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
