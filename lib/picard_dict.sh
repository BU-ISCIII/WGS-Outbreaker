#!/bin/bash
## Author A. Hernández
## version v2.0


if [ $# -eq 0 ];then
        echo -e "\nScript to create dictionary with picartools\n"
        echo -e "Usage: picard_dict.sh reference_path genome_name JAVA_RAM picard_path"
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
JAVA_RAM=$3
picard_path=$4

java $JAVA_RAM -jar $picard_path/picard.jar CreateSequenceDictionary R= $ref_path O= $genome_name.dict
