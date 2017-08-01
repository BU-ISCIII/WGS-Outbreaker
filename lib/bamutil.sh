#!/bin/bash
## Author S. Monzon
## version v2.0


# Test whether the script is being executed with sge or not.
if [ -z $SGE_TASK_ID ]; then
  	use_sge=0
else
  	use_sge=1
fi

# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u

## Usage

if [ $# != 5 -a "$use_sge" == "1" ]; then
  	echo "usage: ............"
  	exit
elif [ $# != 6 -a "$use_sge" == "0" ]; then
  	echo "usage: ............"
   	exit
fi

#Print a trace of simple commands, for commands, case commands, select commands, and arithmetic for commands and their arguments or associated word lists after they are expanded and before they are executed
set -x
echo `date`

# Variables

DIR_BAM=$1
SAMPLE_NAMES=$2
OUTPUT_BAM_SORTED_NAMES=$3
OUTPUT_BAMSTAT_NAMES=$4
EXOME_ENRICHMENT=$5

if [ "$use_sge" = "1" ]; then
  	sample_number=$SGE_TASK_ID
else
  	sample_number=$6
fi

SAMPLE=$( echo $SAMPLE_NAMES | tr ":" "\n" | head -$sample_number | tail -1)
OUTPUT_BAM_SORTED_NAME=$( echo $OUTPUT_BAM_SORTED_NAMES | tr ":" "\n" | head -$sample_number | tail -1)
OUTPUT_BAMSTAT_NAME=$( echo $OUTPUT_BAMSTAT_NAMES | tr ":" "\n" | head -$sample_number | tail -1)

echo "- BAMUtil Analysis: $SAMPLE"
if [ $EXOME_ENRICHMENT == "NO" ];then
	bam stats --in $DIR_BAM/$SAMPLE/$OUTPUT_BAM_SORTED_NAME --basic --baseSum 2> $DIR_BAM/$SAMPLE/$OUTPUT_BAMSTAT_NAME
else
	bam stats --in $DIR_BAM/$SAMPLE/$OUTPUT_BAM_SORTED_NAME --regionList $EXOME_ENRICHMENT --basic --baseSum 2> $DIR_BAM/$SAMPLE/$OUTPUT_BAMSTAT_NAMES
fi
