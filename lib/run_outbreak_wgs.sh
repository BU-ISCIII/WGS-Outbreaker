#!/bin/bash
## Author A.Hernandez
## version v2.0
## usage: run_outbreak_wgs.sh <config.file>

###############
## VARIABLES ##
###############

# Exit immediately if a pipeline, which may consist of a single simple command, a list, or a compound command returns a non-zero status
set -e
# Treat unset variables and parameters other than the special parameters ‘@’ or ‘*’ as an error when performing parameter expansion. An error message will be written to the standard error, and a non-interactive shell will exit
set -u
#Print commands and their arguments as they are executed.
set -x

# Configuration
CONFIG_FILE=$1

#Global VARIABLES
export JOBNAME="exome_pipeline_v2.0"
export SCRIPTS_DIR=$( cat $CONFIG_FILE | grep -w 'SCRIPTS_DIR' | cut -d '=' -f2 )
export TEMP=$( cat $CONFIG_FILE | grep -w 'TEMP_DIR' | cut -d '=' -f2 )
export JAVA_RAM=$( cat $CONFIG_FILE | grep -w 'JAVA_RAM' | cut -d '=' -f2 )

source $SCRIPTS_DIR/processing_config.sh

## SGE args
if [ "$use_sge"="1" ]; then
 mkdir -p $output_dir/logs
  export SGE_ARGS="-V -j y -b y -wd $output_dir/logs -m a -M $email"
fi


# Execute preprocessing
$SCRIPTS_DIR/run_preprocessing.sh

# Execute mapping
if [ $mapping == "YES" ]; then
	$SCRIPTS_DIR/run_mapping.sh
fi

# Execute kmerfinder
if [ $kmerfinder == "YES" ];then
	$SCRIPTS_DIR/run_identification_ST.sh
fi

# Execure srst2
if [ $srst2 == "YES" ]; then
	$SCRIPTS_DIR/run_srst2.sh
fi

#Execute CFSAN
if [ $cfsan == "YES" ]; then
	$SCRIPTS_DIR/run_cfsan.sh
fi

# Execute variant Calling
if [ $variant_calling == "YES" ]; then
        $SCRIPTS_DIR/run_variantCalling_haploid.sh
fi

