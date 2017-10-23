#!/bin/bash

#Author: A.Hernandez

# Help 
# usage: vcf_to_tsv.sh ....
#

# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x

#VARIABLES

dir=$1
vcfsnps_list=$2
tsv_file=$3

bcftools query -f '%CHROM\t%POS\t%REF\t[%TGT\t]\n' $dir/$vcfsnps_list --print-header -o $dir/$tsv_file
