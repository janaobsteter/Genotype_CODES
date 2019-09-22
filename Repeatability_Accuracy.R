#kg mleka
varG <- 8.2
varPE <- 3.9
varTE <- 9

#kg beljakovin
varG <- 0.008
varPE <- 0.004
varTE <- 0.011

varP <- varG + varPE + varTE


t <- (varG + varPE) / varP
t

#t <- 0.49
#t <- 0.35

t = 0.35

r <- function (n, h2, t)  {
  sqrt((n * h2)  / (1 + (n-1)*t) )
}

acc <- data.frame(N = NA, h2 = NA, r = NA)

for (her in c(0.03, 0.1, 0.25, 0.35, 0.5)) {
  for (num in c(1:11)) {
    acc <- rbind(acc, c(num, her, r(num, her, t )))
  }
}

acc <- acc[-1,]
acc
acc$h2 <- as.factor(acc$h2)
acc$N <- as.factor(acc$N)
ggplot(data=acc, aes(x=N, y=r, group=h2, colour=h2)) + geom_line() + ggtitle(paste0("Repeatability: ", round(t, 3))) + xlab("Number of records") + 
  ylab("Accuracy") + scale_y_continuous(breaks=c(seq(0, 1.1, by=0.1), 1))
ggplot(data=acc[acc$h2== 0.25,], aes(x=N, y=r, group=h2, colour=h2)) + geom_line() + ggtitle(paste0("Repeatability: ", round(t, 2))) + xlab("Number of records") + 
  ylab("Accuracy") + scale_y_continuous(breaks=c(seq(0, 0.8, by=0.05), 1))

acc[acc$h2 == 0.25,]

t <- read.csv("~/Documents/PhD/Projects/inProgress/Amount_of_phenotypisation/Heritabilies_repeatabilities_BSWSlo.csv", header=TRUE, skip=1)
#how does heritability change with repeatability
ggplot(data = t, aes(x=h2.1, y=t.1)) + geom_point() + xlab("h2") + ylab("t")
#how does pe2 change with repeatability
ggplot(data = t, aes(x=h2.1, y=pe2.1)) + geom_point() + xlab("h2") + ylab("pe2")
#how do h2 and pe2 change with t
d <- t[,c("h2.1", "pe2.1", "e2.1", "t.1")]
library(reshape)
D <- melt(d, id.vars = "t.1"  )
ggplot(data = D, aes(x=t.1, y=value, group=variable, colour=variable)) + geom_point() + xlab("t") + ylab("Proportion of variance explained") + 
  scale_colour_manual(values = c("red2", "royalblue", "forestgreen"), labels=c("h2", "pe2", "e2" ))

D1 <- melt(d, id.vars = "h2.1" , measure.vars = c("pe2.1", "e2.1", "t.1") )
ggplot(data = D1, aes(x=h2.1, y=value, group=variable, colour=variable)) + geom_point() + xlab("h2") + ylab("Proportion of variance explained") + 
  scale_colour_manual(values = c("red2", "royalblue", "forestgreen"), labels=c("pe2", "e2", "t" ))


#sprememba točnost s številom ponovitev
change <- function (t, n) {sqrt(1 / (t + ((1-t) / n)))}

accChange <- data.frame(t = NA, n=NA, change=NA)
for (t in seq(0.1, 0.8, by=0.1)) {
  for (n in 1:10) {
    accChange <- rbind(accChange, c(t, n, change(t, n)))
  }
}

accChange <- accChange[-1,]
accChange$t <- as.factor(accChange$t)
accChange$n <- as.factor(accChange$n)
ggplot(data=accChange, aes(x=n, y=change, group=t, colour=t)) + geom_line() + ylab("Percentage change in r")
ggplot(data=accChange, aes(x=t, y=change, group=n, colour=n)) + geom_line() + ylab("Percentage change in r")


acc <- read.csv("Accuracy_CatAge_permEnv.csv")
accA <- summarySE(data=acc, measurevar = "COR", groupvars = c("strategy", "scenario", "AgeCat"))

accSel <- accA[accA$AgeCat %in% c("genTest1", "cak5", "vhlevljeni1", "mladi2", "potomciNP0", "telF1"),]
write.csv(accSel, "Accuracy_aggregate_permEnv.csv", quote=FALSE, row.names=FALSE)

#toćnost potomk - telF ob odbiri
#sqrt(.25R^2_oce + .25R^2_mama))

k <- accA[accA$AgeCat %in% c(paste0("k", 3:6)),]
kA <- summarySE(data=k, measurevar = "COR", groupvars = c("strategy", "scenario"))
kA

pb <- accA[accA$AgeCat %in% c(paste0("pb", 5:10)) & accA$scenario=="Class",]
pbA <- summarySE(data=pb, measurevar = "COR", groupvars = c("strategy", "scenario"))

gpb <- accA[accA$AgeCat %in% c(paste0("gpb", 1:5)) & accA$scenario=="Gen",]
gpbA <- summarySE(data=gpb, measurevar = "COR", groupvars = c("strategy", "scenario"))

kC <- kA[kA$scenario == "Class",]
kG <- kA[kA$scenario == "Gen",]

colnames(kC)[4] <- "accK"
colnames(pbA)[4] <- "accPB"
C <- merge(kC[,c(1,2,4)], pbA[,c(1,2,4)], by=c("strategy", "scenario"))
colnames(kG)[4] <- "accK"
colnames(gpbA)[4] <- "accGPB"
G <- merge(kG[,c(1,2,4)], gpbA[,c(1,2,4)], by=c("strategy", "scenario"))


C$accOFF <- sqrt(C$accK^2 + C$accPB^2) / 2
telC <- accA[accA$AgeCat == "telF1" & accA$scenario == "Class",]
C <- C[order(C$strategy),]
C$realAccOff <- telC$COR

G$accOFF <- sqrt(G$accK^2 + G$accGPB^2) / 2
telG <- accA[accA$AgeCat == "telF1" & accA$scenario == "Gen",]
G <- G[order(G$strategy),]
G$realAccOff <- telG$COR


