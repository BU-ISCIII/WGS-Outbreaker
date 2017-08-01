#!/bin/bash
## Author S. Monzon
## version v1.3

source Config.file

FOLDER=$1
FOLDER_OUT=`basename $FOLDER`
mkdir -p $RESULTS_DIR/OUTPUT
mkdir -p $RESULTS_DIR/OUTPUT/$FOLDER_OUT

## Combine indels and snvs in one vcf file.
 vcf-concat $FOLDER/results/passed.somatic.indels.vcf $FOLDER/results/passed.somatic.snvs.filter.vcf > $FOLDER/results/passed.somatic.indels.snvs.vcf \
## Reheader
grep "##" $FOLDER/results/passed.somatic.indels.vcf >> $FOLDER/results/header.txt
grep "##" $FOLDER/results/passed.somatic.snvs.vcf >> $FOLDER/results/header.txt
grep "#CHROM" $FOLDER/results/passed.somatic.snvs.vcf >> $FOLDER/results/header.txt
bcftools reheader -h $FOLDER/results/header.txt $FOLDER/results/passed.somatic.indels.snvs.vcf > $FOLDER/results/passed.somatic.indels.snvs.reheader.vcf

$SCRIPTS_PATH/variant_effect_predictor/variant_effect_predictor.pl -i $FOLDER/results/passed.somatic.indels.snvs.reheader.vcf --format vcf --output_file $FOLDER/results/passed.somatic.indels.snvs.reheader.effect.vcf --everything -cache -dir /home/smonzon/.vep --offline --vcf --force_overwrite
$SCRIPTS_PATH/parserExcell/parserExcell.pl -d $RESULTS_DIR/Strelka -s strelka
	
## Annotate
java -jar $SCRIPTS_PATH/kggseq/kggseq.jar \
	--buildver hg19 \
	--no-gty-vcf-file $FOLDER/results/passed.somatic.indels.snvs.reheader.vcf \
	--indiv-pheno NORMAL:1,TUMOR:1 \
	--db-gene refgene \
	--db-score dbnsfp \
	--genome-annot \
	--db-filter ESP5400,dbsnp141,1kg201305,exac \
	--rare-allele-freq 1 \
	--cosmic-annot \
	--out $RESULTS_DIR/OUTPUT/$FOLDER_OUT/results/passed.somatic.indels.snv.filter.annotated.txt

R --vanilla --slave --args $FOLDER/results/passed.somatic.indels.snvs.reheader.effect.vcf.txt $RESULTS_DIR/OUTPUT/$FOLDER_OUT/results/passed.somatic.indels.snv.filter.annotated.txt.flt.txt $RESULTS_DIR/OUTPUT/$FOLDER_OUT/results < $installDir/lib/merge_output_filtering.R