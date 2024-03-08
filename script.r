#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
setwd("/home/alejandro/PostDoc/heackaton/Charlotte/~/out_dir/")

file1=args[1]
file2=args[2]
DT1=read.table(file1,sep="\t",header=T)
DT2=read.table(file2,sep="\t",header=T)



output <- apply(cbind(DT1[,7], DT2[,7]),MARGIN = 1, sum, na.rm = TRUE)
sum_df <- data.frame(genes=DT1$Geneid, sample=output)


write.table(sum_df,"feature_counts_latest.csv", row.names = FALSE)
