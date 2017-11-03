#!/usr/bin/env Rscript

library(ape)

args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
output <-args[2]
fasta_file <- args[3]
distance <-args[4]
distance_pair <- args[5]

setwd(dir)

#MATRIX SNP DISMATCHES

myfasta <- read.FASTA(fasta_file)

setwd(output)

dis_matrix <- as.data.frame(dist.dna(myfasta, model = "N", as.matrix = TRUE))

dis_matrix$names <- row.names(dis_matrix)

dis_matrix <- data.frame(dis_matrix$names,dis_matrix[,2:ncol(dis_matrix)-1])

write.table(dis_matrix, file= distance, col.names= TRUE, row.names = FALSE, sep = "\t")

#TABLE PAIRWISE SNP DISMATCHES

dis_pair <- dist.dna(myfasta, model = "N", as.matrix = TRUE)

matrix_pair <- as.data.frame(as.table(dis_pair))

write.table(matrix_pair, file= distance_pair, col.names= TRUE, row.names = TRUE, sep = "\t")
