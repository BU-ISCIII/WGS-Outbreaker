#!/bin/bash
## Author: A. Hernandez
# Help
# usage: raxml.sh ....
#


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

raxmlHPC-MPI-AVX -m $model -V -n $output_dir/RAXML_TREE_INFERENCE -p 12345 -s $dir/$snp_msa -N 100

raxmlHPC-MPI-AVX -m $model -V -b 12345 -n $output_dir/RAXML_TREE_BOOTSTRAP -p 12345 -s $dir/$snp_msa -N 100

raxmlHPC-AVX -f b -p 12345  -m $model -t $output_dir/RAxML_bestTree.RAXML_TREE_INFERENCE -z RAxML_bootstrap.RAXML_TREE_BOOTSTRAP -n RAXML_TREE_ANNOT


