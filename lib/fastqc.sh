#!/bin/bash
## Author S. Monzon
## version v2.0

####################
###   Commands  ####
####################

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

if [ $# != 6 -a "$use_sge" == "1" ]; then
	echo "usage: <input dir> <output dir> <sample names (s1:s2:sn)> <fastq files R1 (f1:f2:fn)> <fastq files R2 (f1:f2:fn)> <threads>"
	exit
elif [ $# != 7 -a "$use_sge" == "0" ]; then                                                                         
	echo "usage: <input dir> <output dir> <sample names (s1:s2:sn)> <fastq files R1 (f1:f2:fn)> <fastq files R2 (f1:f2:fn)> <threads> <sample_number>" 
 	exit
fi                                                                                                          

#Print a trace of simple commands, for commands, case commands, select commands, and arithmetic for commands and their arguments or associated word lists after they are expanded and before they are executed
set -x
echo `date`

## Variables
DIR=$1
OUTPUT_DIR=$2
THREADS=$6      
SAMPLE_NAMES=$3
FASTQ_FILES_R1=$4
FASTQ_FILES_R2=$5


if [ "$use_sge" = "1" ]; then
	sample_number=$SGE_TASK_ID
else
	sample_number=$7
fi


SAMPLE=$( echo $SAMPLE_NAMES | tr ":" "\n" | head -$sample_number | tail -1)
FASTQ_R1=$( echo $FASTQ_FILES_R1 | tr ":" "\n" | head -$sample_number | tail -1) 
FASTQ_R2=$( echo $FASTQ_FILES_R2 | tr ":" "\n" | head -$sample_number | tail -1)  

echo -e "Running FastQC for $SAMPLE....\n"

# Results folder per sample creation
mkdir -p $OUTPUT_DIR/QC/fastqc/$SAMPLE

fastqc --noextract -o $OUTPUT_DIR/QC/fastqc/$SAMPLE -t $THREADS $DIR/$SAMPLE/$FASTQ_R1 $DIR/$SAMPLE/$FASTQ_R2	 

echo -e "Preprocessing for $SAMPLE finished \n\n"
