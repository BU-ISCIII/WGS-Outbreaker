#!/bin/bash
## Author:A.Hernandez
## version v2.0

if [ $# -eq 0 ];then
        echo -e "\nScript to create directories for each sample for cfsan pipeline\n"
        echo -e "Usage: cfsan_create_sample_directories.sh samples_list input_dir sample_count"
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

samples=$1
dir=$2
sample_count=$3

touch $dir/sampleDirectories.txt
sample=$( echo $samples | tr ":" "\n" | head -$sample_count | tail -1)
echo $dir/samples/$sample >> $dir/sampleDirectories.txt
