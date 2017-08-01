#!/bin/bash
## Author S. Monzon
## version v1.3

source Config.file 

####################
###   Commands  ####
####################

echo -e "Prepare reference: bwa and samtools indexing"
$BWA_PATH/bwa index -a bwtsw $REF_PATH
samtools faidx $REF_PATH