args = commandArgs(TRUE)

nc= args[3]
data<-read.table(nc, header=T, sep='\t', quote=NULL, comment='')
over_fdr<-p.adjust(data$pvalue, method="BH")
a<-cbind(data, over_fdr)
na=paste(args[1],'/',args[2],'_FETresults_FDR.csv', sep='')
write.csv(a, file=na, row.names=F, quote = F)