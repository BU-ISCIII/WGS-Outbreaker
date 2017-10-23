#!/bin/bash

##Author: A. Hernandez
#help
#Usage: cfsan_readgroup.sh

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

#VARIABLES

threads=$1
dir=$2
samples=$3
dedup_sam_list=$4
dedup_bam_list=$5
platform=$6
picard_path=$7

if [ "$use_sge" = "1" ]; then                                                                                                                                                                                               
   	sample_count=$SGE_TASK_ID                                                                    
else                                                                                                        
   	sample_count=$8                                                                                
fi

sample=$( echo $samples | tr ":" "\n" | head -$sample_count | tail -1)
dedup_sam_file=$( echo $dedup_sam_list | tr ":" "\n" | head -$sample_count | tail -1)
dedup_bam_file=$( echo $dedup_bam_list | tr ":" "\n" | head -$sample_count | tail -1)

java $JAVA_RAM -jar $picard_path/picard.jar AddOrReplaceReadGroups I="$dir/$sample/$dedup_sam_file" O="$dir/$sample/$dedup_bam_file" SO=coordinate ID=S$sample LB=L$sample PL=$platform PU=$sample SM=$sample TMP_DIR="$TEMP"
