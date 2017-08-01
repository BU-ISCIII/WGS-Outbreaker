#!/bin/bash
## Author S. Monzon
## version v1.3

source Config.file 

## Workflow Strelka
date
echo -e "Running TrioCaller...\n"
mkdir $RESULTS_DIR/triocaller/
mkdir $RESULTS_DIR/triocaller/samtools_variants/
mkdir $RESULTS_DIR/triocaller/triocaller_corrected_variants/
cp $installDir/etc/triocaller.ped $RESULTS_DIR/triocaller

for f in `ls $RESULTS_DIR/Alignment/*wodup.bam`
do
	file=`basename $RESULTS_DIR/Alignment/$f`
	#java -Djava.io.tmpdir=$TEMP -Xmx24G -jar $SCRIPTS_PATH/picard/AddOrReplaceReadGroups.jar I=$f O=$RESULTS_DIR/Alignment/$file-RG.bam RGID=HWI-ST731_6 RGLB=undefined RGPL=Illumina RGPU=undefined RGSM=$file VALIDATION_STRINGENCY=LENIENT  
	#samtools index $RESULTS_DIR/Alignment/$file-RG.bam
done

#samtools mpileup -Iuf $REF_PATH $RESULTS_DIR/Alignment/*RG.bam | bcftools view -bvcg - > $RESULTS_DIR/triocaller/samtools_variants/variants.bcf

#bcftools view $RESULTS_DIR/triocaller/samtools_variants/variants.bcf  > $RESULTS_DIR/triocaller/samtools_variants/variants.vcf

#### Refining genotypes using LD information and trio constraints ####
#awk -f $SCRIPTS_PATH/separate_indels_snps_vcf.awk $RESULTS_DIR/triocaller/samtools_variants/variants.vcf
$SCRIPTS_PATH/TrioCaller-06262012/TrioCaller --vcf $RESULTS_DIR/triocaller/samtools_variants/variants_snv.vcf --pedfile $RESULTS_DIR/triocaller/triocaller.ped --rounds 10 --prefix $RESULTS_DIR/triocaller/triocaller_corrected_variants/variants.triocaller.vcf

echo -e "TrioCaller finished.\n"