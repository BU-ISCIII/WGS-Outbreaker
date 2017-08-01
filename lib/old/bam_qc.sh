#!/bin/bash
## Author S. Monzon
## version v1.3

source Config.file 

####################
###   Commands  ####
####################
mkdir $RESULTS_DIR/QC/bam_qc
mkdir $RESULTS_DIR/QC/bam_qc/preduplicates
mkdir $RESULTS_DIR/QC/bam_qc/wo_duplicates

## BAM stats using teqc pre duplicate filter
R --vanilla --slave --args $RESULTS_DIR/Alignment/ -align.bam$ $EXOME_ENRICHMENT $RESULTS_DIR/QC/bam_qc/preduplicates < $SCRIPTS_PATH/teqc.R 
## BAM stats using teqc post duplicate filter
R --vanilla --slave --args $RESULTS_DIR/Alignment/ -wodup.bam$ $EXOME_ENRICHMENT $RESULTS_DIR/QC/bam_qc/wo_duplicates < $SCRIPTS_PATH/teqc.R 
