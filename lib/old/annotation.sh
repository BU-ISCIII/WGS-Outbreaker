#!/bin/bash
## Author S. Monzon
## version v1.3

source Config.file

mkdir $RESULTS_DIR/KGGSEQ_ANNOTATION

## Annotation script. It uses kggseq software and annotates genomics and proteomics information from known databases. Also it filters mutations according to several disease models given a genotype.

## SNPS ##
## De novo model. Annotations: refgene,dbsnp,effect predictions(SIFT,LRT,Mutation Assesor,etc.),uniprot,filter per allele freq, pubmed and omim.
java -jar /opt/kggseq/kggseq.jar -Xmx24G --buildver hg19 --vcf-file $RESULTS_DIR/gatk_trio/Variants/trio.snps.indels.filtered.vcf \
									  --no-lib-check \
									  --no-resource-check \
									  --ped-file $PED_FILE \
									  --genotype-filter 7 \
									  --db-gene refgene \
									  --db-score dbnsfp \
									  --genome-annot \
									  --rare-allele-freq 0.005 \
									  --omim-annot \
									  --out $RESULTS_DIR/KGGSEQ_ANNOTATION/kggseq_snp_annotation_denovo

java -jar /opt/kggseq/kggseq.jar -Xmx24G --buildver hg19 --vcf-file $RESULTS_DIR/gatk_trio/Variants/trio.snps.indels.rephased.vcf \
									  --no-lib-check \
									  --no-resource-check \
									  --ped-file $PED_FILE \
									  --genotype-filter 7 \
									  --db-gene refgene \
									  --db-score dbnsfp \
									  --genome-annot \
									  --rare-allele-freq 0.005 \
									  --omim-annot \
									  --out $RESULTS_DIR/KGGSEQ_ANNOTATION/kggseq_snp_phase_annotation_denovo

## Double-hit-gene-trio-filter
java -jar /opt/kggseq/kggseq.jar -Xmx24G --buildver hg19 --vcf-file $RESULTS_DIR/gatk_trio/Variants/trio.snps.indels.filtered.vcf \
									  --no-lib-check \
									  --no-resource-check \
									  --ped-file $PED_FILE \
									  --double-hit-gene-trio-filter \
									  --db-gene refgene \
									  --db-score dbnsfp \
									  --genome-annot \
									  --rare-allele-freq 0.005 \
									  --omim-annot \
									  --out $RESULTS_DIR/KGGSEQ_ANNOTATION/kggseq_snp_annotation_doublehit

java -jar /opt/kggseq/kggseq.jar -Xmx24G --buildver hg19 --vcf-file $RESULTS_DIR/gatk_trio/Variants/trio.snps.indels.rephased.vcf \
									  --no-lib-check \
									  --no-resource-check \
									  --ped-file $PED_FILE \
									  --double-hit-gene-trio-filter \
									  --db-gene refgene \
									  --db-score dbnsfp \
									  --genome-annot \
									  --rare-allele-freq 0.005 \
									  --omim-annot \
									  --out $RESULTS_DIR/KGGSEQ_ANNOTATION/kggseq_snp_phase_annotation_doublehit

## Recessive
# java -jar /opt/kggseq/kggseq.jar -Xmx24G --buildver hg19 --vcf-file $RESULTS_DIR/gatk_trio/Variants/trio.snps.indels.filtered.vcf \
# 									  --no-lib-check \
# 									  --no-resource-check \
# 									  --ped-file $PED_FILE \
# 									  --genotype-filter 1 \
# 									  --db-gene refgene \
# 									  --db-score dbnsfp \
# 									  --genome-annot \
# 									  --rare-allele-freq 0.005 \
# 									  --omim-annot \
# 									  --out $RESULTS_DIR/KGGSEQ_ANNOTATION/kggseq_snp_annotation_recessive

# java -jar /opt/kggseq/kggseq.jar -Xmx24G --buildver hg19 --vcf-file $RESULTS_DIR/gatk_trio/Variants/trio.snps.indels.rephased.vcf \
# 									  --no-lib-check \
# 									  --no-resource-check \
# 									  --ped-file $PED_FILE \
# 									  --genotype-filter 1 \
# 									  --db-gene refgene \
# 									  --db-score dbnsfp \
# 									  --genome-annot \
# 									  --rare-allele-freq 0.005 \
# 									  --omim-annot \
# 									  --out $RESULTS_DIR/KGGSEQ_ANNOTATION/kggseq_snp_phase_annotation_recessive

## Compound-heterozygosity
# java -jar /opt/kggseq/kggseq.jar -Xmx24G --buildver hg19 --vcf-file $RESULTS_DIR/gatk_trio/Variants/trio.snps.indels.filtered.vcf \
# 									  --no-lib-check \
# 									  --no-resource-check \
# 									  --ped-file $PED_FILE \
# 									  --genotype-filter 2 \
# 									  --db-gene refgene \
# 									  --db-score dbnsfp \
# 									  --genome-annot \
# 									  --rare-allele-freq 0.005 \
# 									  --omim-annot \
# 									  --out $RESULTS_DIR/KGGSEQ_ANNOTATION/kggseq_snp_annotation_compoundhet

# java -jar /opt/kggseq/kggseq.jar -Xmx24G --buildver hg19 --vcf-file $RESULTS_DIR/gatk_trio/Variants/trio.snps.indels.rephased.vcf \
# 									  --no-lib-check \
# 									  --no-resource-check \
# 									  --ped-file $PED_FILE \
# 									  --genotype-filter 2 \
# 									  --db-gene refgene \
# 									  --db-score dbnsfp \
# 									  --genome-annot \
# 									  --rare-allele-freq 0.005 \
# 									  --omim-annot \
# 									  --out $RESULTS_DIR/KGGSEQ_ANNOTATION/kggseq_snp_phase_annotation_compoundhet

## Dominant without consanguineous mating
# java -jar /opt/kggseq/kggseq.jar -Xmx24G --buildver hg19 --vcf-file $RESULTS_DIR/gatk_trio/Variants/trio.snps.indels.filtered.vcf \
# 									  --no-lib-check \
# 									  --no-resource-check \
# 									  --ped-file $PED_FILE \
# 									  --genotype-filter 3 \
# 									  --db-gene refgene \
# 									  --db-score dbnsfp \
# 									  --genome-annot \
# 									  --rare-allele-freq 0.005 \
# 									  --omim-annot \
# 									  --out $RESULTS_DIR/KGGSEQ_ANNOTATION/kggseq_snp_annotation_dominantwoconsang

# java -jar /opt/kggseq/kggseq.jar -Xmx24G --buildver hg19 --vcf-file $RESULTS_DIR/gatk_trio/Variants/trio.snps.indels.rephased.vcf \
# 									  --no-lib-check \
# 									  --no-resource-check \
# 									  --ped-file $PED_FILE \
# 									  --genotype-filter 3 \
# 									  --db-gene refgene \
# 									  --db-score dbnsfp \
# 									  --genome-annot \
# 									  --rare-allele-freq 0.005 \
# 									  --omim-annot \
# 									  --out $RESULTS_DIR/KGGSEQ_ANNOTATION/kggseq_snp_phase_annotation_dominantwoconsang

## Dominant with full penetrance causal mutation(s) in the sample
# java -jar /opt/kggseq/kggseq.jar -Xmx24G --buildver hg19 --vcf-file $RESULTS_DIR/gatk_trio/Variants/trio.snps.indels.filtered.vcf \
# 									  --no-lib-check \
# 									  --no-resource-check \
# 									  --ped-file $PED_FILE \
# 									  --genotype-filter 4 \
# 									  --db-gene refgene \
# 									  --db-score dbnsfp \
# 									  --genome-annot \
# 									  --rare-allele-freq 0.005 \
# 									  --omim-annot \
# 									  --out $RESULTS_DIR/KGGSEQ_ANNOTATION/kggseq_snp_annotation_dominantfullpen

# java -jar /opt/kggseq/kggseq.jar -Xmx24G --buildver hg19 --vcf-file $RESULTS_DIR/gatk_trio/Variants/trio.snps.indels.rephased.vcf \
# 									  --no-lib-check \
# 									  --no-resource-check \
# 									  --ped-file $PED_FILE \
# 									  --genotype-filter 4 \
# 									  --db-gene refgene \
# 									  --db-score dbnsfp \
# 									  --genome-annot \
# 									  --rare-allele-freq 0.005 \
# 									  --omim-annot \
# 									  --out $RESULTS_DIR/KGGSEQ_ANNOTATION/kggseq_snp_phase_annotation_dominantfullpen

## Dominant without heterogeneity
# java -jar /opt/kggseq/kggseq.jar -Xmx24G --buildver hg19 --vcf-file $RESULTS_DIR/gatk_trio/Variants/trio.snps.indels.filtered.vcf \
# 									  --no-lib-check \
# 									  --no-resource-check \
# 									  --ped-file $PED_FILE \
# 									  --genotype-filter 5 \
# 									  --db-gene refgene \
# 									  --db-score dbnsfp \
# 									  --genome-annot \
# 									  --rare-allele-freq 0.005 \
# 									  --omim-annot \
# 									  --out $RESULTS_DIR/KGGSEQ_ANNOTATION/kggseq_snp_annotation_dominantwohet

# java -jar /opt/kggseq/kggseq.jar -Xmx24G --buildver hg19 --vcf-file $RESULTS_DIR/gatk_trio/Variants/trio.snps.indels.rephased.vcf \
# 									  --no-lib-check \
# 									  --no-resource-check \
# 									  --ped-file $PED_FILE \
# 									  --genotype-filter 5 \
# 									  --db-gene refgene \
# 									  --db-score dbnsfp \
# 									  --genome-annot \
# 									  --rare-allele-freq 0.005 \
# 									  --omim-annot \
# 									  --out $RESULTS_DIR/KGGSEQ_ANNOTATION/kggseq_snp_phase_annotation_dominantfwohet

## Dominant without heterogeneity
# java -jar /opt/kggseq/kggseq.jar -Xmx24G --buildver hg19 --vcf-file $RESULTS_DIR/gatk_trio/Variants/trio.snps.indels.filtered.vcf \
# 									  --no-lib-check \
# 									  --no-resource-check \
# 									  --ped-file $PED_FILE \
# 									  --genotype-filter 6 \
# 									  --db-gene refgene \
# 									  --db-score dbnsfp \
# 									  --genome-annot \
# 									  --rare-allele-freq 0.005 \
# 									  --omim-annot \
# 									  --out $RESULTS_DIR/KGGSEQ_ANNOTATION/kggseq_snp_annotation_fullpen

# java -jar /opt/kggseq/kggseq.jar -Xmx24G --buildver hg19 --vcf-file $RESULTS_DIR/gatk_trio/Variants/trio.snps.indels.rephased.vcf \
# 									  --no-lib-check \
# 									  --no-resource-check \
# 									  --ped-file $PED_FILE \
# 									  --genotype-filter 6 \
# 									  --db-gene refgene \
# 									  --db-score dbnsfp \
# 									  --genome-annot \
# 									  --rare-allele-freq 0.005 \
# 									  --omim-annot \
# 									  --out $RESULTS_DIR/KGGSEQ_ANNOTATION/kggseq_snp_phase_annotation_fullpen

#############
### Parseo y anotación sin modelo de herencia
#############

java -jar /opt/kggseq/kggseq.jar -Xmx24G --buildver hg19 --vcf-file $RESULTS_DIR/gatk_trio/Variants/trio.snps.rephased.vcf \
									  --no-lib-check \
									  --no-resource-check \
									  --ped-file $PED_FILE \
									  --db-gene refgene \
									  --db-score dbnsfp \
									  --genome-annot \
									  --rare-allele-freq 0.005 \
									  --omim-annot \
									  --out $RESULTS_DIR/KGGSEQ_ANNOTATION/kggseq_snp_phase_annotation_nomodel