#!/bin/bash

##Author: A.Hernandez
#help
#usage: cfsan_sort_sam.sh


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
sort_sam_list=$4
picard_path=$5

if [ "$use_sge" = "1" ]; then                                                                                                                                                                                               
   	sample_count=$SGE_TASK_ID                                                                      
else                                                                                                        
   	sample_count=$6                                                                                   
fi

sample=$( echo $samples | tr ":" "\n" | head -$sample_count | tail -1)
sort_sam_file=$( echo $sort_sam_list | tr ":" "\n" | head -$sample_count | tail -1)

java $JAVA_RAM -jar $picard_path/picard.jar SortSam INPUT="$dir/$sample/reads.sam" OUTPUT="$dir/$sample/$sort_sam_file" SO=coordinate TMP_DIR="$TEMP"
