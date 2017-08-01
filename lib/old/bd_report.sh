#!/bin/bash
## Author S. Monzon
## version v1.3

source Config.file 

####################
###   Commands  ####
####################
echo -e "Generating BD and Variants Report"

dat=`date +"%Y%m%d_%H%M"`

mkdir $RESULTS_DIR/"$dat""_Report_Variants"
cd $RESULTS_DIR/"$dat""_Report_Variants"

perl $installDir/variants_database/generateDB.pl -db "$dat""_mutations.sqlite"
perl $installDir/variants_database/includetier.pl -d $RESULTS_DIR - t $THREADS
perl $installDir/variants_database/parserDB.pl -data $RESULTS_DIR -db "$dat""_mutations.sqlite"
perl $installDir/variants_database/generateReport.pl -db "$dat""_mutations.sqlite"




