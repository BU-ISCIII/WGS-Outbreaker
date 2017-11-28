#!/bin/bash
## Author A. Hernández
## version v2.0


if [ $# -eq 0 ];then
        echo -e "\nScript to crete references for bowtie2\n"
        echo -e "Usage: reference_bowtie.sh reference_path genome_name"
        exit
fi

# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x

#VARIABLES

ref_path=$1
genome_name=$2

bowtie2-build $ref_path $genome_name
