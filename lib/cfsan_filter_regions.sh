#!/bin/bash
#author:A. Hernandez
#help
#usage: cfsan_filter_regions.sh


# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x


#VARIABLES

dir=$1
cfsan_ref_path=$2

cfsan_snp_pipeline filter_regions -n var.flt.vcf $dir/sampleDirectories.txt $cfsan_ref_path


