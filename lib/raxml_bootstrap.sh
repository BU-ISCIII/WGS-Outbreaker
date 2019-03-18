#!/bin/bash
## Author: A. Hernandez
## version v2.0

if [ $# -eq 0 ];then
        echo -e "\nScript to run RAxML bootstrap\n"
        echo -e "Usage: raxml_bootstrap.sh input_dir output_dir snp_matrix model_raxml"
        exit
fi

if [ -z $JOB_ID ]; then
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
output_dir=$2
snp_msa=$3
model=$4
boots=$5

if [ "$use_sge" = "1" ]; then
	raxmlHPC-MPI-AVX -m $model -V -b 12345 -w $output_dir -n RAXML_TREE_BOOTSTRAP -p 12345 -s $dir/$snp_msa -N $boots
else
	raxmlHPC-AVX -m $model -V -b 12345 -w $output_dir -n RAXML_TREE_BOOTSTRAP -p 12345 -s $dir/$snp_msa -N $boots
fi
