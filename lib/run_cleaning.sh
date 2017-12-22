#!/bin/bash
## Author: A.Hernandez
## version v2.0

if [ $# -eq 0  ]; then
        echo -e "\nExecute cleaning\u"
        echo "Usage: run_cleaning.sh <config.file>"
        exit
fi

#Execute processing_config.sh

CONFIG_FILE=$1

# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x

clean_cmd="$SCRIPTS_DIR/clean.sh\
	$output_dir \
	$mapping \
	$duplicate_filter \
	$cfsan \
	$srst2" 

#In HPC
if [ "$use_sge" = "1" ]; then
	clean_args=${SGE_ARGS} -pe openmp $threads}
	clean=$( qsub clean_args -N $JOBNAME.CLEAN $clean_cmd)
	else
	echo "Running clean mapping"
	clean=$( $clean_cmd)
fi
