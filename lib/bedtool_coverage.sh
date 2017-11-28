#!/bin/bash
## Author A.Hernandez
## version v2.0

if [ $# -eq 0 ]; then
        echo -e "\nScrip to run bedtools genomecov\n"
        echo -e "Usage: bedtool_coverage.sh input_dir output_dir samples_list duplicateBAM_list coverage_list coverage_graph_list reference_path"
        exit
fi

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

#VARIABLES

dir=$1
output_dir=$2
samples=$3
duplicateBamArray_list=$4
coverage_list=$5
coverage_graph_list=$6
ref_path=$7

if [ "$use_sge" = "1" ]; then
        sample_count=$SGE_TASK_ID
else
        sample_count=$8
fi

sample=$( echo $samples | tr ":" "\n" | head -$sample_count | tail -1)
duplicateBamArray=$( echo $duplicateBamArray_list | tr ":" "\n" | head -$sample_count | tail -1)
coverage_sample=$( echo $coverage_list | tr ":" "\n" | head -$sample_count | tail -1)
coverage_graph_sample=$( echo $coverage_graph_list | tr ":" "\n" | head -$sample_count | tail -1)

bedtools genomecov -ibam $dir/$sample/$duplicateBamArray -g $ref_path > $output_dir/$coverage_sample
