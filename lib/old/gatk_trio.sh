#!/bin/bash
## Author S. Monzon
## version v1.3

source Config.file 

#Creación de carpetas de resultados
echo -e "Starting Variant Calling workflow"
date
echo -e "\n"
mkdir $RESULTS_DIR/gatk_trio
mkdir $RESULTS_DIR/gatk_trio/Realignment/
mkdir $RESULTS_DIR/gatk_trio/Recalibration/
mkdir $RESULTS_DIR/gatk_trio/Variants/

for file in `ls $RESULTS_DIR/Alignment/*-wodup.bam`
do
   name=`basename $file`
	echo -e "1) Local Realignment: $file"	

   # Generar información de read en el fichero BAM, es solicitado por GATK.
   java -Djava.io.tmpdir=$TEMP $java_ram -jar $SCRIPTS_PATH/picard/AddOrReplaceReadGroups.jar I=$file O=$file-RG.bam RGID=HWI-ST699 RGLB=undefined RGPL=Illumina RGPU=undefined RGSM=$name VALIDATION_STRINGENCY=LENIENT  

   # Indexar el nuevo fichero BAM generado.
   $SAMTOOLS_PATH/samtools index $file-RG.bam

   ## Realineamiento local - preparación.
   java -Djava.io.tmpdir=$TEMP $java_ram -jar $SCRIPTS_PATH/gatk/GenomeAnalysisTK.jar \
   		-T RealignerTargetCreator \
   		-I $file-RG.bam \
      -known $INDELS_1000G \
   		-R $REF_PATH \
   		-o $RESULTS_DIR/gatk_trio/$name-IndelRealigner.intervals \
   		-nt $THREADS \
   		-S LENIENT \
   		-log $RESULTS_DIR/gatk_trio/$name-targetCreator.log

   ## Realineamiento alrededor de indels.
	java -Djava.io.tmpdir=$TEMP $java_ram -jar $SCRIPTS_PATH/gatk/GenomeAnalysisTK.jar \
   		-T IndelRealigner \
   		-I $file-RG.bam \
      -known $INDELS_1000G \
   		-R $REF_PATH \
   		-targetIntervals $RESULTS_DIR/gatk_trio/$name-IndelRealigner.intervals \
   		-o $RESULTS_DIR/gatk_trio/Realignment/$name-realigned.bam \
   		-S LENIENT \
   		-log $RESULTS_DIR/gatk_trio/$name-realigner.log

   echo -e "2) Base Quality Recalibration"

   ## Recalibración de la calidad de las bases.
   java -Djava.io.tmpdir=$TEMP $java_ram -jar $SCRIPTS_PATH/gatk/GenomeAnalysisTK.jar \
	    -T BaseRecalibrator \
	    -R $REF_PATH \
   	  -I $RESULTS_DIR/gatk_trio/Realignment/$name-realigned.bam \
      -knownSites $SNPS_hg19 \
      -knownSites $HAPMAP_33 \
   	  -cov ReadGroupCovariate \
   	  -cov QualityScoreCovariate \
   	  -cov CycleCovariate \
   	  -cov ContextCovariate \
      -o $RESULTS_DIR/gatk_trio/Recalibration/$name-recal1_data.grp \
   	  -nct $THREADS \
   	  -S LENIENT \
   	  -log $RESULTS_DIR/gatk_trio/$name-recal.log

    java -Djava.io.tmpdir=$TEMP $java_ram -jar $SCRIPTS_PATH/gatk/GenomeAnalysisTK.jar \
      -T BaseRecalibrator \
      -BQSR $RESULTS_DIR/gatk_trio/Recalibration/$name-recal1_data.grp \
      -R $REF_PATH \
      -I $RESULTS_DIR/gatk_trio/Realignment/$name-realigned.bam \
      -knownSites $SNPS_hg19 \
      -knownSites $HAPMAP_33 \
      -cov ReadGroupCovariate \
      -cov QualityScoreCovariate \
      -cov CycleCovariate \
      -cov ContextCovariate \
      -o $RESULTS_DIR/gatk_trio/Recalibration/$name-recal2_data.grp \
      -nct $THREADS \
      -S LENIENT \
      -log $RESULTS_DIR/gatk_trio/$name-recal.log

    java -Djava.io.tmpdir=$TEMP $java_ram -jar $SCRIPTS_PATH/gatk/GenomeAnalysisTK.jar \
      -T AnalyzeCovariates \
      -R $REF_PATH \
      -before $RESULTS_DIR/gatk_trio/Recalibration/$name-recal1_data.grp \
      -after $RESULTS_DIR/gatk_trio/Recalibration/$name-recal2_data.grp \
      -csv $RESULTS_DIR/gatk_trio/Recalibration/$name-BQSR.csv \
      -plots $RESULTS_DIR/gatk_trio/Recalibration/$name-BQSR.pdf \
      -log $RESULTS_DIR/gatk_trio/$name-analyzecovariates.log

    java -Djava.io.tmpdir=$TEMP $java_ram -jar $SCRIPTS_PATH/gatk/GenomeAnalysisTK.jar \
      -T PrintReads \
      -R $REF_PATH \
      -I $RESULTS_DIR/gatk_trio/Realignment/$name-realigned.bam \
      -BQSR $RESULTS_DIR/gatk_trio/Recalibration/$name-recal1_data.grp \
      -o $RESULTS_DIR/gatk_trio/Recalibration/$name-recalibrated.bam \
      -nct $THREADS \
      -S LENIENT \
      -log $RESULTS_DIR/gatk_trio/$name-print_recal.log
done

echo -e "Variant calling with HaplotypeCaller"

find -P $RESULTS_DIR/gatk_trio/Recalibration -name "*recalibrated.bam" > $RESULTS_DIR/gatk_trio/Recalibration/bam.list

# Variant calling. Sólo variantes.
java -Djava.io.tmpdir=$TEMP $java_ram -jar $SCRIPTS_PATH/gatk/GenomeAnalysisTK.jar \
     -T HaplotypeCaller \
     -R $REF_PATH \
     -I $RESULTS_DIR/gatk_trio/Recalibration/bam.list \
     --dbsnp $SNPS_hg19 \
     -o $RESULTS_DIR/gatk_trio/Variants/trio.snps.indels.raw.vcf \
     -stand_call_conf 30.0 \
     -stand_emit_conf 10.0 \
     -S LENIENT \
     -log $RESULTS_DIR/gatk_trio/$name-variants.log

java -Djava.io.tmpdir=$TEMP $java_ram -jar $SCRIPTS_PATH/gatk/GenomeAnalysisTK.jar \
     -R $REF_PATH \
     -T SelectVariants \
     -V $RESULTS_DIR/gatk_trio/Variants/trio.snps.indels.raw.vcf \
     -o $RESULTS_DIR/gatk_trio/Variants/trio.snps.raw.vcf \
     -selectType SNP \
     -nt $THREADS \
     -S LENIENT \
     -log $RESULTS_DIR/gatk_trio/$name-selectSNP.log

# Se filtran según los parámetros por defecto que recomienda GATK. El filtro por cobertura es muy bajo
# y el verdadero filtro se realizará después de realizar la comparación tumor frente a control. 
java -Djava.io.tmpdir=$TEMP $java_ram -jar $SCRIPTS_PATH/gatk/GenomeAnalysisTK.jar \
     -R $REF_PATH \
     -T VariantFiltration \
     -V $RESULTS_DIR/gatk_trio/Variants/trio.snps.raw.vcf \
     -o $RESULTS_DIR/gatk_trio/Variants/trio.snps.filtered.vcf \
     --clusterWindowSize 10 \
     --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
     --filterName "HARD_TO_VALIDATE" \
     --filterExpression "DP <5 " \
     --filterName "LowCoverage" \
     --filterExpression "QUAL <30.0 " \
     --filterName "VeryLowQual" \
     --filterExpression "QUAL >30.0 && QUAL <50.0 " \
     --filterName "LowQual" \
     --filterExpression "QD <1.5 " \
     --filterName "LowQD" \
     --filterExpression "SB >-10.0 " \
     --filterName "StrandBias" \
     --filterExpression "FS >60.0 " \
     --filterName "p-value StrandBias" \
     --filterExpression "HaplotypeScore >13.0 " \
     --filterName "HaplotypeScore" \
     -S LENIENT \
     -log $RESULTS_DIR/gatk_trio/$name-filterSNPs.log

   ## Realizamos el mismo procedimiento con los indels.
echo -e "Select and Filter Indels"
java -Djava.io.tmpdir=$TEMP $java_ram -jar $SCRIPTS_PATH/gatk/GenomeAnalysisTK.jar \
     -R $REF_PATH \
     -T SelectVariants \
     -V $RESULTS_DIR/gatk_trio/Variants/trio.snps.indels.raw.vcf \
     -o $RESULTS_DIR/gatk_trio/Variants/trio.indels.raw.vcf \
     -selectType INDEL \
     -nt $THREADS \
     -S LENIENT \
     -log $RESULTS_DIR/gatk_trio/trio.selectIndels.log

java -Djava.io.tmpdir=$TEMP $java_ram -jar $SCRIPTS_PATH/gatk/GenomeAnalysisTK.jar \
     -T VariantFiltration \
     -R $REF_PATH \
     -V $RESULTS_DIR/gatk_trio/Variants/trio.indels.raw.vcf \
     -o $RESULTS_DIR/gatk_trio/Variants/trio.indels.filtered.vcf \
     --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
     --filterName "LowQD|StrandBias|ReadPosRankSum" \
     --filterExpression "QD < 2.0 || MQ0 > 50" \
     --filterName "Nov09filters" \
     -S LENIENT \
     -log $RESULTS_DIR/gatk_trio/trio.filterIndels.log

echo -e "Between Sample phasing"  
java -Djava.io.tmpdir=$TEMP $java_ram -jar $SCRIPTS_PATH/gatk/GenomeAnalysisTK.jar \
     -R $REF_PATH \
     -T PhaseByTransmission \
     -V $RESULTS_DIR/gatk_trio/Variants/trio.snps.filtered.vcf \
     -ped $installDir/etc/triocaller.ped \
     -o $RESULTS_DIR/gatk_trio/Variants/trio.snps.filtered.phased.vcf \
     -log $RESULTS_DIR/gatk_trio/trio.snps.PhaseByTransmission.log

java -Djava.io.tmpdir=$TEMP $java_ram -jar $SCRIPTS_PATH/gatk/GenomeAnalysisTK.jar \
     -R $REF_PATH \
     -T PhaseByTransmission \
     -V $RESULTS_DIR/gatk_trio/Variants/trio.indels.filtered.vcf \
     -ped $installDir/etc/triocaller.ped \
     -o $RESULTS_DIR/gatk_trio/Variants/trio.indels.filtered.phased.vcf \
     -log $RESULTS_DIR/gatk_trio/trio.indels.PhaseByTransmission.log


echo -e "Across Sample phasing"   
java -Djava.io.tmpdir=$TEMP $java_ram -jar $SCRIPTS_PATH/gatk/GenomeAnalysisTK.jar \
     -T ReadBackedPhasing \
     -R $REF_PATH \
     -I $RESULTS_DIR/gatk_trio/Recalibration/bam.list \
     --variant $RESULTS_DIR/gatk_trio/Variants/trio.snps.filtered.phased.vcf \
     -L $RESULTS_DIR/gatk_trio/Variants/trio.snps.filtered.phased.vcf \
     -o $RESULTS_DIR/gatk_trio/Variants/trio.snps.rephased.vcf \
     --phaseQualityThresh 20.0 \
     -log $RESULTS_DIR/gatk_trio/trio.snps.ReadBackedPhasing.log

java -Djava.io.tmpdir=$TEMP $java_ram -jar $SCRIPTS_PATH/gatk/GenomeAnalysisTK.jar \
     -T ReadBackedPhasing \
     -R $REF_PATH \
     -I $RESULTS_DIR/gatk_trio/Recalibration/bam.list \
     --variant $RESULTS_DIR/gatk_trio/Variants/trio.indels.filtered.phased.vcf \
     -L $RESULTS_DIR/gatk_trio/Variants/trio.indels.filtered.phased.vcf \
     -o $RESULTS_DIR/gatk_trio/Variants/trio.indels.rephased.vcf \
     --phaseQualityThresh 20.0 \
     -log $RESULTS_DIR/gatk_trio/trio.indels.ReadBackedPhasing.log
