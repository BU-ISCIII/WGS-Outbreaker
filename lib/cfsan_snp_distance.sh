#!/bin/bash
## Author: A. Hernandez
## version v2.0

if [ $# -eq 0 ];then
        echo -e "\nScript to run cfsan distance\n"
        echo -e "Usage: cfsan_snp_distance.sh input_dir"
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

cfsan_snp_pipeline distance -p $dir/snp_distance_pairwise.tsv -m $dir/snp_distance_matrix.tsv $dir/snpma.fasta
cfsan_snp_pipeline distance -p $dir/snp_distance_pairwise_preserved.tsv -m $dir/snp_distance_matrix_preserved.tsv $dir/snpma_preserved.fasta