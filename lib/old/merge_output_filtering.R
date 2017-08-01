#!/bin/bash
## Author S. Monzon
## version v1.3

args <- commandArgs(trailingOnly = TRUE)
variants_file <- args[1]
variants_annot_file <- args[2]
results_dir <- args[3]

##Merge and filtering variants
variants <- read.table(variants_file,header=T,sep="\t")
variants_annot <- read.table(variants_annot_file,header=T,sep="\t")

variants$merged <- paste(variants$CHROM,variants$POS,sep="_")
variants_annot$merged <- paste("chr",variants_annot$Chromosome,"_",variants_annot$StartPosition,sep="")

table_comp <- merge(variants,variants_annot,by="merged",all.x=T,all.y=T)
write.table(table_comp,file=paste(results_dir,"/passed.somatic.indes.snvs.annotated.all.txt",sep=""),sep="\t",row.names=F)

