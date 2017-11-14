#!/bin/bash
## Author A. Hernández
## version v2.0

# Help 
# usage: run_references.sh ....
#

CONFIG_FILE=$1

source $SCRIPTS_DIR/processing_config.sh --"$CONFIG_FILE"

# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
#set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x



bwa_cmd="$SCRIPTS_DIR/reference_bwa.sh \
	$ref_path"

bowtie_cmd="$SCRIPTS_DIR/reference_bowtie.sh \
	$ref_path
	$genome_name"

picard_cmd="$SCRIPTS_DIR/picard_dict.sh \
	$ref_path
	$genome_name
	$JAVA_RAM
	$picard_path"

samtools_cmd="$SCRIPTS_DIR/samtools_fai.sh \
	$ref_path"

if [ "$use_sge" = "1" ]; then

	if [ ! -f $ref_path.amb ]; then

		bwa_index=$( qsub $SGE_ARGS -N BWA_INDEX $bwa_cmd)
	fi

	if [ ! -f $genome_name.1.bt2 ]; then

        	bowtie_index=$( qsub $SGE_ARGS -N BOWTIE_INDEX $bowtie_cmd)
	fi
 
	if [ ! -f $genome_name.dict ]; then
	
        	picard_dict=$( qsub $SGE_ARGS -pe openmp $threads -l h_vmem=$vmem -N PICARD_DICT $picard_cmd)
	fi

	if [ ! -f $ref_path.fai ]; then

       		samtools_fai=$( qsub $SGE_ARGS -N SAMTOOLS_FAI $samtools_cmd)
	fi

else
	if [ ! -f $ref_path.amb ]; then
		bwa_index=$($bwa_cmd)
	fi 

	if [ ! -f $ref_path.1.bt2 ]; then
		bowtie_index=$($bowtie_cmd)
	fi

	if [ ! -f $ref_path.dict ]; then
		picard_dict=$($picard_cmd)
	fi
	if [ ! -f $ref_path.fai ]; then
		samtools_fai=$($samtools_cmd)
	fi
fi
