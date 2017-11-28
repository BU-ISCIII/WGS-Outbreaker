#!/bin/bash
## Author:A. Hernandez
## version v2.0

if [ $# -eq 0 ];then
        echo -e "\nScript to run cfsan filter_regions\n"
        echo -e "Usage: cfsan_filter_regions.sh input_dir reference_path edge_length max_snp window_size"
        exit
fi


# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x

#VARIABLES

dir=$1
cfsan_ref_path=$2
edge_length=$3
max_snp=$4
window_size=$5

cfsan_snp_pipeline filter_regions --edge_length $edge_length --window_size $window_size --max_snp $max_snp -n var.flt.vcf $dir/sampleDirectories.txt $cfsan_ref_path
