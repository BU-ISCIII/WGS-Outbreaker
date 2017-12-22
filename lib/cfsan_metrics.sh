#!/bin/bash
## Author: A. Hernandez
## version v2.0

if [ $# -eq 0 ];then
        echo -e "\nScript to run cfsan collect_metrics\n"
        echo -e "Usage: cfsan_metrics.sh input_dir samples_list reference_path"
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

#VARIABLES

dir=$1
samples=$2
cfsan_ref_path=$3

if [ "$use_sge" = "1" ]; then
	sample_count=$SGE_TASK_ID                                                                   
else                                                                                                        
   	sample_count=$4                                                                               
fi

sample=$( echo $samples | tr ":" "\n" | head -$sample_count | tail -1)

cfsan_snp_pipeline collect_metrics -o $dir/samples/$sample/metrics $dir/samples/$sample $cfsan_ref_path
