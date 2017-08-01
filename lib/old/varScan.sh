#!/bin/bash
## Author S. Monzon
## version v1.3

source Config.file 

mkdir $RESULTS_DIR/Variants/
mkdir $RESULTS_DIR/ControlVsTumor


for f in $RESULTS_DIR/Alignment/*.bam
do

filename=`basename $f`
echo -e "3) Variant Calling: $filename"

refname=`basename $REF_PATH`
refdir=`dirname $REF_PATH`
$SAMTOOLS_PATH/samtools mpileup -d 500 -u -l $EXOME_ENRICHMENT -Q 30 -q 40 -C50 -f $REF_PATH $f > $RESULTS_DIR/Variants/variants_$filename.pileup
$SCRIPTS_PATH/pileline/cmd/pileline-rfilter.sh -A $RESULTS_DIR/Variants/variants_$filename.pileup ./ -b $EXOME_ENRICHMENT > $RESULTS_DIR/Variants/variants_$filename.pileupOntarget 
$SCRIPTS_PATH/samtools.pl varFilter -D100 $RESULTS_DIR/Variants/variants_$filename.pileupOntarget | awk '($3=="*"&&$6>=50)||($3!="*"&&$6>=20)' > $RESULTS_DIR/Variants/variants_$filename.pileupOntarget_filt 
wait

done


for f in `ls $RESULTS_DIR/Variants/ | grep "^variants_.*$Control.*.pileupOntarget_filt"` 
do
	echo -e "4) Variant Analysis"

	filename=`basename $f`
	sample=$(perl $SCRIPTS_PATH/processII.pl $filename $Control)

	vcontrol="variants_$sample$Control-wodup.bam.pileupOntarget_filt"
	vcase="variants_$sample$Case-wodup.bam.pileupOntarget_filt"

	java -jar $SCRIPTS_PATH/VarScan.v2.2.11.jar somatic $RESULTS_DIR/Variants/$vcontrol $RESULTS_DIR/Variants/$vcase $RESULTS_DIR/$filename-varscan

done