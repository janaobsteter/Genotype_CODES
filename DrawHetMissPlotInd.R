#!/usr/bin/env Rscript

library("geneplotter")

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("One argument must be supplied (input file)", call.=FALSE)
}

imiss=read.table(paste(args[1],".imiss",sep=""),h=T)
imiss$logF_MISS = log10(imiss[,6])
het=read.table(paste(args[1],".het",sep=""),h=T)
het$meanHet = (het$N.NM. - het$O.HOM.)/het$N.NM.
colors  <- densCols(imiss$logF_MISS,het$meanHet,colramp = colorRampPalette(c("grey", "black")))
pdf(paste(args[1], "_imiss-vs-het.pdf", sep=""), width=13/2.54, height=13/2.54, pointsize=12)
plot(imiss$logF_MISS,het$meanHet, col=colors, xlim=c(-5,0),ylim=c(0,1),pch=20, xlab="Proportion of missing genotypes", ylab="Heterozygosity rate",axes=F,cex.lab=1.5)
axis(2,at=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),tick=T)
axis(1,at=c(-5,-4,-3,-2,-1,0),labels=c("0.00001","0.0001","0.001","0.01","0.1","1"))
abline(h=mean(het$meanHet, na.rm=TRUE)-(6*sd(het$meanHet, na.rm=TRUE)),col="BLACK",lty=2)
abline(h=mean(het$meanHet, na.rm=TRUE)+(6*sd(het$meanHet, na.rm=TRUE)),col="BLACK",lty=2)
#abline(v=-1.522879, col="BLUE", lty=2)
abline(v=log10(0.10), col="BLACK", lty=2)
dev.off()

MISS <- data.frame(imiss[imiss$F_MISS > 0.10, c("FID","IID")])
HET  <- data.frame(het[het$meanHet > (mean(het$meanHet, na.rm=TRUE)+(6*sd(het$meanHet, na.rm=TRUE))) | het$meanHet < (mean(het$meanHet)-(6*sd(het$meanHet))), c("FID","IID")])

names(MISS) <- c("ID")
names(HET)  <- c("ID")

EXCLUDE <- unique(rbind(MISS,HET))

write.table(EXCLUDE, file="IndividualsToExlcudeMissHet.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)
