#!/bin/bash
#author:A. Hernandez
#help
#usage: cfsan_filter_regions.sh


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


#VARIABLES






if [ "$use_sge" = "1" ]; then
        sample_count=$sge_task_id
else
        sample_count=$8
fi


cfsan_snp_pipeline filter_regions -n  sample


