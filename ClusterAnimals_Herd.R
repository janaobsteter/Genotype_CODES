ped <- read.table("/home/jana/Documents/PhD/CompBio/TestingGBLUP/PedigreeAndGeneticValues_cat.txt", header=TRUE)
#ped$herd <- kmeans(ped[,c("Mother", "gvNormUnres1", "phenoNormUnres1")], 100)

pedK <- ped[ped$cat=="k",c("Indiv", "Mother", "gvNormUnres1", "phenoNormUnres1")]

pedK$cluster <- NA
pedK$cluster <- kmeans(pedK[,2:4], 100)$cluster
write.table(pedK, "/home/jana/Documents/PhD/CompBio/TestingGBLUP/PedCows_HERDS.txt", quote=FALSE, row.names=FALSE)
herdNo <- as.data.frame(table(pedK$cluster))
colnames(herdNo) <- c("Herd", "NoAnim")
write.table(herdNo, "/home/jana/Documents/PhD/CompBio/HerdNo.txt", quote=FALSE, row.names=FALSE)

gen <- read.csv("/home/jana/Documents/PhD/CompBio/TestingGBLUP/IndForGeno.txt", header=FALSE)
length(intersect(gen$V1, pedK$Indiv))

chosen  <- read.csv("/home/jana/Documents/PhD/CompBio/TestingGBLUP/ChosenInd.txt", header=FALSE, sep=" ")
length(intersect(chosen$V15, pedK$Indiv))



####################################
#Newstart
setwd("/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenSLO/")
setwd("/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenSLO/TE")
sol <- read.table("solutions", header=FALSE, skip=1)
colnames(sol) [3] <- "renID"
colnames(sol) [4] <- "EBV"
ped <- read.table("SimulatedData/PedigreeAndGeneticValues_cat.txt", header=TRUE)
pedB <- read.table("renadd01.ped")
length(which(pedB$V1 != pedB$V10))
colnames(pedB)[1] <- "renID"
colnames(pedB)[10] <- "Indiv"
pedBlup <- read.table("Blupf90.ped")


ped$EBV <- renSol$V3
cor(ped$EBV, ped$gvNormUnres1)


#najprej združi rešitev in renadd
pedB <- merge(pedB, sol, by="renID")
ped <- merge(ped, pedB, by="Indiv")
cor(ped$gvNormUnres1, ped$EBV)
cor(ped$EBV, ped$gvNormUn)

#v zadnji generaciji
ped58 <- ped[ped$Generation==59,]
cor(ped58$gvNormUnres1, ped58$EBV) # to je normalno - parent average

#v geotipiziranih živalih
genInd <- read.csv("IndForGeno.txt", header=FALSE)
pedG <- ped[ped$Indiv %in% genInd$V1, ]
cor(pedG$EBV, pedG$gvNormUnres1)

pedNew <- pedG[pedG$Generation ==60,]
cor(pedNew$gvNormUnres1, pedNew$EBV)


#Use blup.dat from before
dat <- read.table("Blupf90.dat")
colnames(dat)[1] <- "Indiv"
datP <- merge(dat, ped, by="Indiv", all.x=T)



#################3
#preveri v sel file
setwd("/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenSLO/")
ped <- read.table("SimulatedData/PedigreeAndGeneticValues_cat.txt", header=TRUE)
sol <- read.table("TESTBLUP/KravepoGen/renumbered_Solutions")
sol <- read.table("TESTBLUP//renumbered_Solutions")
colnames(sol) <- c("renID", "Indiv", "EBV")
SOL <- merge(sol, ped, by="Indiv")
cor(SOL$gvNormUnres1, SOL$EBV)
solR <- read.table("TESTBLUP/KravepoGen/renumbered_Solutions")
solO <- read.table("renumbered_Solutions_60")
solbPB <- read.table("./TESTBLUP/brezPB/renumbered_Solutions")
ped$EBV <- sol$V3
ped$EBVR <- solR$V3
ped$EBVO <- solO$V3
ped$EBVF <- solbPB$V3
ped$EBV2 <- sol$V3

cor(ped$EBV, ped$gvNormUnres1) #aktivne krave + pb
cor(ped$EBVR, ped$gvNormUnres1) #razporejene
cor(ped$EBVO, ped$gvNormUnres1) #original iz simulacije
cor(ped$EBVF, ped$gvNormUnres1) #no progeny tested sires, only females
cor(ped$EBV2, ped$gvNormUnres1) #no progeny tested sires, only females

#v zadnji generaciji
ped58 <- ped[ped$Generation==60,]
SOL58 <- SOL[SOL$Generation==60,]
cor(ped58$gvNormUnres1, ped58$EBV) # to je normalno - parent average
cor(ped58$gvNormUnres1, ped58$EBVR) # to je normalno - parent average
cor(ped58$gvNormUnres1, ped58$EBVO) # to je normalno - parent average
cor(ped58$gvNormUnres1, ped58$EBVF) # to je normalno - parent average
cor(ped58$gvNormUnres1, ped58$EBV2) # to je normalno - parent average
cor(SOL58$gvNormUnres1, SOL58$EBV) # to je normalno - parent average

#v geotipiziranih živalih
genIndR <- read.csv("./TESTBLUP/KravepoGen/IndForGeno.txt", header=FALSE)
genInd <- read.csv("./TESTBLUP/KraveAktivne/IndForGeno.txt", header=FALSE)
genIndF <- read.csv("./TESTBLUP/brezPB/IndForGeno.txt", header=FALSE)
genIndO <- read.csv("./TESTBLUP/IndForGeno.txt", header=FALSE)
pedGR <- ped58[ped58$Indiv %in% genIndR$V1, ] #razporejene
pedGO <- ped58[ped58$Indiv %in% genIndO$V1, ] #original
pedGF <- ped58[ped58$Indiv %in% genIndF$V1, ] #original
pedG <- ped58[ped58$Indiv %in% genInd$V1, ] #aktivne
SOLG <- SOL58[SOL58$Indiv %in% genIndO$V1, ] #aktivne
cor(pedGR$EBVR, pedGR$gvNormUnres1)
cor(pedGO$EBVO, pedGO$gvNormUnres1)
cor(pedG$EBV, pedG$gvNormUnres1)
cor(pedGF$EBVF, pedGF$gvNormUnres1)
cor(SOLG$EBV, SOLG$gvNormUnres1)

library(plyr)

func <- function(ped)
{
  return(data.frame(COR = cor(ped$EBV, ped$gvNormUnres1)))
}

genCor <- ddply(ped, .(ped$Generation), func)
cor <- ggplot(data=genCor, aes(x=genCor$`ped$Generation`, y=genCor$COR)) + geom_point() +ylim(-0.1,1)

func <- function(ped)
{
  return(data.frame(COR = cor(ped$EBVR, ped$gvNormUnres1)))
}

genCorR <- ddply(ped, .(ped$Generation), func)
corR <- ggplot(data=genCorR, aes(x=genCorR$`ped$Generation`, y=genCorR$COR)) + geom_point() + ylim(-0.1,1)

func <- function(ped)
{
  return(data.frame(COR = cor(ped$EBVO, ped$gvNormUnres1)))
}

genCorO <- ddply(ped, .(ped$Generation), func)
corO <- ggplot(data=genCorO, aes(x=genCorO$`ped$Generation`, y=genCorO$COR)) + geom_point()+ ylim(-0.1,1)

func <- function(ped)
{
  return(data.frame(COR = cor(ped$EBVF, ped$gvNormUnres1)))
}

genCorF <- ddply(ped, .(ped$Generation), func)
corF <- ggplot(data=genCorF, aes(x=genCorF$`ped$Generation`, y=genCorF$COR)) + geom_point()+ ylim(-0.1,1)
library(Rmisc)
multiplot(corO, corR, corF, cor, cols=4)


#Use blup.dat from before
dat <- read.table("Blupf90.dat")
colnames(dat)[1] <- "Indiv"
datP <- merge(dat, ped, by="Indiv", all.x=T)



#to so genotipizirane živali po kategorijah
ind <- read.table("./TESTBLUP/KraveAktivne//IndForGeno.txt")
colnames(ind)[1] <- "Indiv"
ind <- merge(ind, ped, by="Indiv")
ind$cat <- as.factor(ind$cat)
cat <- ggplot(ind, aes(x=cat, fill=cat)) + geom_bar()
gen <- ggplot(ind[ind$Generation %in% 45:60,], aes(x=Generation, fill=Generation)) + geom_bar() + xlim(45, 60)
table(ind$cat)
table(ind$Generation)
sum(table(ind$Generation))
nrow(ind)

length(intersect(dat$Indiv, ind$Indiv))
nrow(dat[dat$Indiv %in% ind$Indiv,])
datI <- datP[datP$Indiv %in% ind$Indiv,]
catP <- ggplot(datI, aes(x=cat, fill=cat)) + geom_bar()+ ylim(0, 25000)
genP <- ggplot(datI[datI$Generation %in% 45:60,], aes(x=Generation, fill=Generation)) + geom_bar() + xlim(45, 60)+ ylim(0, 7000)

#to so genotipizirane živali po kategorijah
ind <- read.table("./TESTBLUP/KravepoGen/IndForGeno.txt")
colnames(ind)[1] <- "Indiv"
ind <- merge(ind, ped, by="Indiv")
ind$cat <- as.factor(ind$cat)
catR <- ggplot(ind, aes(x=cat, fill=cat)) + geom_bar()
genR <- ggplot(ind[ind$Generation %in% 45:60,], aes(x=Generation, fill=Generation)) + geom_bar() + xlim(45, 60)
table(ind$cat)
table(ind$Generation)
nrow(ind)

length(intersect(dat$Indiv, ind$Indiv))
nrow(dat[dat$Indiv %in% ind$Indiv,])
datI <- datP[datP$Indiv %in% ind$Indiv,]
catPR <- ggplot(datI, aes(x=cat, fill=cat)) + geom_bar()+ ylim(0, 25000)
genPR <- ggplot(datI[datI$Generation %in% 45:60,], aes(x=Generation, fill=Generation)) + geom_bar() + xlim(45, 60)+ ylim(0, 7000)


###


ind <- read.table("./TESTBLUP/brezPB//IndForGeno.txt")
colnames(ind)[1] <- "Indiv"
ind <- merge(ind, ped, by="Indiv")
ind$cat <- as.factor(ind$cat)
catF <- ggplot(ind, aes(x=cat, fill=cat)) + geom_bar()
genF <- ggplot(ind[ind$Generation %in% 45:60,], aes(x=Generation, fill=Generation)) + geom_bar() + xlim(45, 60)
table(ind$cat)
table(ind$Generation)
nrow(ind)

length(intersect(dat$Indiv, ind$Indiv))
nrow(dat[dat$Indiv %in% ind$Indiv,])
datI <- datP[datP$Indiv %in% ind$Indiv,]
catPF <- ggplot(datI, aes(x=cat, fill=cat)) + geom_bar()+ ylim(0, 25000)
genPF <- ggplot(datI[datI$Generation %in% 45:60,], aes(x=Generation, fill=Generation)) + geom_bar() + xlim(45, 60)+ ylim(0, 7000)


###



ind <- read.table("IndForGeno.txt")
colnames(ind)[1] <- "Indiv"
ind <- merge(ind, ped, by="Indiv")
ind$cat <- as.factor(ind$cat)
catO <- ggplot(ind, aes(x=cat, fill=cat)) + geom_bar()
genO <- ggplot(ind[ind$Generation %in% 45:60,], aes(x=Generation, fill=Generation)) + geom_bar() + xlim(45, 60)
table(ind$cat)
table(ind$Generation)
nrow(ind)

length(intersect(dat$Indiv, ind$Indiv))
nrow(dat[dat$Indiv %in% ind$Indiv,])
datI <- datP[datP$Indiv %in% ind$Indiv,]
catPO <- ggplot(datI, aes(x=cat, fill=cat)) + geom_bar() + ylim(0, 25000)
genPO <- ggplot(datI[datI$Generation %in% 45:60,], aes(x=Generation, fill=Generation)) + geom_bar() + xlim(45, 60) + ylim(0, 7000)

###


multiplot(catO, catR, catF,cat, cols=4)
multiplot(genO, genR, genF, gen, cols=4)
multiplot(catPO, catPR, catPF,catP, cols=4)
multiplot(genPO, genPR, genPF, genP, cols=4)
###

#tukaj pa poglej še, koliko imaš fenotipskih podatkov po kategorijah in generaijcah
##########


setwd("~/Documents/PhD/CompBio/TestingGBLUP/")
dat <- read.table("Blupf90.dat")
ped <- read.table("PedigreeAndGeneticValues_cat.txt", header=TRUE)
colnames(dat) <- c("Indiv", "Ph", "sex")

DAT <- merge(dat, ped, by="Indiv")
