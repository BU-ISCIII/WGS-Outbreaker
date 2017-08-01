#!/bin/bash
## Author S. Monzon
## version v1.3

source Config.file 

## Workflow Strelka
date
echo -e "Running Strelka Workflow...\n"
mkdir $RESULTS_DIR/Strelka/
cd $RESULTS_DIR/Strelka/

# Step 2. Copy configuration ini file from default template set to a
#         local copy, possibly edit settings in local copy of file:
#
cp $STRELKA_INSTALL_DIR/strelka/etc/strelka_config_bwa_default.ini config.ini

for f in `ls $RESULTS_DIR/Alignment/ | grep "^.*$Control-wodup.bam$"` 
do
	echo -e "4) Variant Analysis"

	filename=`basename $f`
	sample=$(perl $installDir/requisites/processII.pl $filename $Control)

	bamcontrol="$sample$Control-wodup.bam"
	echo $bamcontrol
	bamcase="$sample$Case-wodup.bam"
	echo $bamcase

	# Step 3. Configure:
	
	$STRELKA_INSTALL_DIR/configureStrelkaWorkflow.pl --normal=$RESULTS_DIR/Alignment/$bamcontrol --tumor=$RESULTS_DIR/Alignment/$bamcase --ref=$REF_PATH --config=config.ini --output-dir=$RESULTS_DIR/Strelka/$filename-Analysis

	# Step 4. Run Analysis	
	
	cd ./$filename-Analysis
	make -j $THREADS
	cd ./results

	## Filtrar por enriquecimiento de exoma
	# vcftools --vcf passed.somatic.snvs.vcf --out $filename-Strelka_ExomeFilter --bed $EXOME_ENRICHMENT --recode --recode-INFO-all
	# vcftools --vcf passed.somatic.indels.vcf --out $filename-Strelka_indel_ExomeFilter --bed $EXOME_ENRICHMENT --recode --recode-INFO-all

	grep -vP "^chr.*_gl.*(_random)*" passed.somatic.snvs.vcf > passed.somatic.snvs.filter.vcf

	# $installDir/requisites/variant_effect_predictor/variant_effect_predictor.pl -i passed.somatic.snvs.filter.vcf --format vcf --output_file $filename-Strelka_snv_Effect.vcf --everything -cache -dir $cacheDir --offline --vcf --force_overwrite
	# $installDir/requisites/variant_effect_predictor/variant_effect_predictor.pl -i passed.somatic.indels.vcf --format vcf --output_file $filename-Strelka_indel_Effect.vcf --everything -cache -dir $cacheDir --offline --vcf --force_overwrite 

	# perl $installDir/requisites/parserExcell.pl -d $filename-Strelka_snv_Effect.vcf -s strelka
	# perl $installDir/requisites/parserExcell.pl -d $filename-Strelka_indel_Effect.vcf -s strelkaindel

	cd ..
	cd ..

done


