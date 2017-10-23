#!/bin/bash

##Author: A. Hernandez
#help
#Usage: cfsan_combine_metrics.sh


# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x

#VARAIBLES

dir=$1

cfsan_snp_pipeline combine_metrics -n metrics -o $dir/metrics.tsv $dir/sampleDirectories.txt
