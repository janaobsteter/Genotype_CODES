#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("One argument must be supplied (input file)", call.=FALSE)
}

x<-read.table(paste(args[1],".lmiss",sep=""),header=T)
try(x[is.na(x$F_MISS),]$F_MISS <- 0, silent=TRUE)
try(x[x$F_MISS<=0.0001,]$F_MISS <- 0.00011)
ylabels=c("0","5K","10K","15K","20K")
xlabels=c("0.0001","0.001","0.01","0.1","1")
#par(mfrow=c(1,1))
pdf(paste(args[1], "_lmiss.pdf", sep=""), width=13/2.54, height=13/2.54, pointsize=12)
#hist(log10(x$F_MISS),axes=F,xlim=c(-4,0),col="RED",ylab="Number of SNPs",xlab="Fraction of missing data",main="All SNPs",ylim=c(0,20000))
hist(log10(x$F_MISS),xaxt="n",xlim=c(-4,0),col="RED",ylab="Number of SNPs",xlab="Fraction of missing data",main=NULL,cex.lab=1.5)
axis(side=2,labels=F)
#mtext(ylabels,side=2,las=2, at=c(0,5000,10000,15000,20000),line=1)
axis(side=1,labels=F)
mtext(xlabels,side=1,at=c(-4,-3,-2,-1,0),line=1)
abline(v=log10(0.10),lty=2)
dev.off()
