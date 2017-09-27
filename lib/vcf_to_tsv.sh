#!/bin/bash

#Author: A.Hernandez

# Help 
# usage: run_cfsan.sh ....
#

# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x

#VARIABLES

threads=$1
dir=$2
vcfsnpsfilArray_list=$3
tsv_file=$4

bcftools query -f '%CHROM\t%POS\t%REF\t[%TGT\t]\n' $dir/$vcfsnpsfilArray_list --print-header -o $dir/$tsv_file
