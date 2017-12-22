#!/bin/bash
## Author S. Monzon & A. Hernandez
## version v2.0

if [ $# -eq 0 ];then
        echo -e "\nScript to mark duplicates with picardtools\n"
        echo -e "Usage: picard_duplicates.sh input_dir samples_list BAM_rg_list BAM_duplicate_list picard_path"
        exit
fi

# Test whether the script is being executed with sge or not
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
samples=$2
mappingArray_rg_list=$3
duplicateBamArray_list=$4
picard_path=$5                                                                                                                                                                                                                   
if [ "$use_sge" = "1" ]; then
	sample_count=$SGE_TASK_ID                                                                     
else                                                                                                        
   	sample_count=$6
fi                                              

sample=$( echo $samples | tr ":" "\n" | head -$sample_count | tail -1)
mappingArray_rg=$( echo $mappingArray_rg_list | tr ":" "\n" | head -$sample_count | tail -1) 
duplicateBamArray=$( echo $duplicateBamArray_list | tr ":" "\n" | head -$sample_count | tail -1) 

echo -e "- Duplicate Filter: $sample"

java $JAVA_RAM -jar $picard_path/picard.jar MarkDuplicates ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false INPUT="$input_dir/$sample/$mappingArray_rg" OUTPUT="$input_dir/$sample/$duplicateBamArray" METRICS_FILE="$input_dir/$sample/$sample.duplicates.stats" TMP_DIR="$TEMP"

#Index BAM file
samtools index $input_dir/$sample/$duplicateBamArray
