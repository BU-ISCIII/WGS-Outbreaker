#!/bin/bash
## Author S. Monzon
## version v2.0
## usage bamutil.sh

if [ $# -eq 0 ]; then
	echo -e "\nScrip to run bam stats\n"
	echo -e "Usage: bamutils.sh output_dir samples_list bamRG_list bamStat_list exome_enrichement"
	exit
fi

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
output_dir=$1
samples=$2
mappingArray_rg_list=$3
bamstatArray_pre_list=$4
exome_enrichement=$5

if [ "$use_sge" = "1" ]; then
  	sample_count=$SGE_TASK_ID
else
  	sample_count=$6
fi

sample=$( echo $samples | tr ":" "\n" | head -$sample_count | tail -1)
mappingArray_rg=$( echo $mappingArray_rg_list | tr ":" "\n" | head -$sample_count | tail -1)
bamstatArray_pre=$( echo $bamstatArray_pre_list | tr ":" "\n" | head -$sample_count | tail -1)

echo "- BAMUtil Analysis: $sample"
if [ $exome_enrichement == "NO" ];then
	bam stats --in $output_dir/$sample/$mappingArray_rg --basic --baseSum 2> $output_dir/$sample/$bamstatArray_pre
else
	bam stats --in $output_dir/$sample/$mappingArray_rg --regionList $exome_enrichement --basic --baseSum 2> $output_dir/$sample/$bamstatArray_pre
fi
