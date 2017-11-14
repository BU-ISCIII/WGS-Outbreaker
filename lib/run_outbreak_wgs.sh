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

export SCRIPTS_DIR=$( cat $CONFIG_FILE | grep -w 'SCRIPTS_DIR' | cut -d '=' -f2 )

source $SCRIPTS_DIR/processing_config.sh --"$CONFIG_FILE"

# Execute preprocessing
$SCRIPTS_DIR/run_preprocessing.sh $CONFIG_FILE

#Execute references

$SCRIPTS_DIR/run_references.sh $CONFIG_FILE


# Execute mapping
if [ $mapping == "YES" ]; then
	$SCRIPTS_DIR/run_mapping.sh $CONFIG_FILE
fi

# Execute kmerfinder
if [ $kmerfinder == "YES" ];then
	$SCRIPTS_DIR/run_identification_ST.sh $CONFIG_FILE

fi

# Execure srst2
if [ $srst2 == "YES" ]; then
	$SCRIPTS_DIR/run_srst2.sh $CONFIG_FILE

fi

#Execute CFSAN
if [ $cfsan == "YES" ]; then
	$SCRIPTS_DIR/run_cfsan.sh $CONFIG_FILE
fi

# Execute variant Calling
if [ $variant_calling == "YES" ]; then
        $SCRIPTS_DIR/run_variantCalling_haploid.sh $CONFIG_FILE
fi

#Execute RAxML
if [ $raxml == "YES" ]; then
	$SCRIPTS_DIR/run_raxml.sh $CONFIG_FILE

fi

#Execute stats

if [ $stats == "YES" ]; then
	 $SCRIPTS_DIR/run_stats.sh $CONFIG_FILE

fi
