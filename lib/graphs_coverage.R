#!/usr/bin/env Rscript

library(ggplot2)
library(plyr)

args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]

setwd(dir)

# get file names
print(files <- list.files(pattern="coverage.csv$"))

# load and parse the files

cov_graph <- NULL
pdf(file="coverage_graph.pdf",width=15)
maxCov=200

for (f in files){
 	df=read.table(f, sep="\t")
	colnames(df) <- c("chr","covThreshold","fractionAtThisCoverage","ChrLength","diffFracBelowThreshold")
 	cov <- ddply(df,.(chr),summarize,covThreshold=covThreshold,fracAboveThreshold=1-cumsum(diffFracBelowThreshold))

	cov_graph <- rbind(cov_graph,cbind(cov, sample= gsub("^(AA-[0-9a-zA-Z]+).*$", "\\1", f)))


	p2<-ggplot(subset(cov, covThreshold<maxCov & chr == "genome"), aes(x= covThreshold, y=100* fracAboveThreshold)) +
geom_line() +
ylim(0, 100) +
theme_bw() +
theme(axis.text.x = element_text(size = 10.5,angle=75, vjust=0.5), strip.text.x = element_text(size=6.5)) +
labs(title=paste("Genome Coverage",f,sep=" "), x="Depth of coverage", y="Percentage of coverage")

	plot(p2)
}

dev.off()

cov_genome <- cov_graph[cov_graph$chr == "genome",]
cov_table <- by(cov_genome,cov_genome[,"sample"],function(x) x$fracAboveThreshold[x$covThreshold == 10 | x$covThreshold==20 | x$covThreshold== 30 | x$covThreshold == 50| x$covThreshold == 70| x$covThreshold == 100])
cov_table <- do.call(rbind,cov_table)
colnames(cov_table) = c("10x","20x","30x","50x","70x","100x")
cov_table <- as.data.frame(cov_table)
cov_table$names <- row.names(cov_table)
cov_table <- data.frame(cov_table$names,cov_table[,2:ncol(cov_table)-1])
write.table(cov_table, file= "coverage_table.csv", col.names= TRUE, row.names = FALSE, sep = "\t")
