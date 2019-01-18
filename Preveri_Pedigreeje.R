ped <- read.csv("~/Class_20Rep.txt",sep=" ")
ped <- read.csv("~/Class_1Year.txt",sep=" ")
ped <- read.csv("~/GenSLO5_1Pb.txt",sep=" ")
ped <- read.csv("~/GenSLO15_1Pb.txt",sep=" ")
ped <- read.csv("~/Class_1Year.txt",sep=" ")
ped <- read.csv("~/Class_1Year.txt",sep=" ")
ped <- read.csv("~/PEDGI.txt",sep=" ")



pedL <- ped[ped$Generation %in% 40:60,]

#age of sires of cakajoci biki
Cak <- ped[ped$cat=="cak", c("Indiv", "Generation", "Father")]
colnames(Cak) <- c("Male", "MaleBirth", "Sire")

CakSire <- ped[ped$Indiv %in% Cak$Sire, c("Indiv", "Generation")]
colnames(CakSire) <- c("Sire", "SireBirth")

CAKSIRE <- merge(Cak, CakSire, by="Sire")
CAKSIRE$SireAge <- CAKSIRE$MaleBirth - CAKSIRE$SireBirth
CAKSIRE[order(CAKSIRE$MaleBirth),]

#najdi stare očete
ped[ped$Indiv %in% c(397358, 397394),]
#najdi mlade očete čakajočih iste generacije
ped[ped$Indiv ==431974,]
#vsi očetje te generacije
ped[ped$cat=="pb" & ped$Generation %in% c(46:50),]


table(pedL$Father, pedL$Generation)

ped50 <- ped[ped$Generation==50,]
table(ped$Generation[ped$Indiv %in% ped50$Father])
table(ped50$Father)

FatherUse <- data.frame()
for (father in unique(pedL$Father)) {
  range <- max(ped$Generation[ped$Father == father]) - min(ped$Generation[ped$Father == father])
  FatherUse <- rbind(FatherUse, c(father, range))
}
colnames(FatherUse) <- c("Father", "NoOffspring")

FatherUse$NoOffspring <- as.numeric(FatherUse$NoOffspring)
hist(FatherUse$NoOffspring)
table(FatherUse$NoOffspring)


ped[ped$Indiv %in% FatherUse$Father[FatherUse$NoOffspring > 5],]
table(ped$cat[ped$Indiv %in% FatherUse$Father[FatherUse$NoOffspring > 5]])
table(ped$Generation[ped$Indiv %in% FatherUse$Father[FatherUse$NoOffspring > 5]])

ped[ped$Indiv %in% FatherUse$Father[FatherUse$NoOffspring == 0],]
table(ped$cat[ped$Indiv %in% FatherUse$Father[FatherUse$NoOffspring == 0]])
table(ped$Generation[ped$Indiv %in% FatherUse$Father[FatherUse$NoOffspring == 0]])


#
#tukaj pa rabiš preverit, koliko očetov je odbranih na eno generacijo
#še vedno
ped <- read.csv("~/Class_1Pb.txt",sep=" ")
ped <- read.csv("~/Gen_1Pb.txt",sep=" ")
ped <- read.csv("~/GenSLO_1Pb.txt",sep=" ")
ped <- read.csv("~/BmGen_1Pb.txt",sep=" ")

table(ped$Generation[ped$cat == 'pb'])
table(ped$Generation[ped$cat == 'mladi'])


#koliko časa so očetje v uporabi
#kdaj jih odbereš
table(ped$Generation[ped$cat=="pb"])
table(ped$Generation[ped$cat=="cak"])
table(ped$Generation[ped$cat=="mladi"])
table(ped$Generation[ped$cat=="vhlevljeni"])
table(ped$Generation[ped$cat=="potomciNP"])

ped <- ped[ped$Generation %in% 40:60,]
for (father in unique(ped$Indiv[ped$cat=="pb"])) {
  offspring <- unique(ped$Generation[ped$Father==father])
  diff <- max(offspring) - min(offspring)
  print(c(father, diff))
}


ped[ped$Father==406080,]
use <- unique(ped$Generation[ped$Father==406080])
for (year in use) {
  off <- nrow(ped[ped$Generation==year & ped$Father==406080,])
  print(c(year, off))
}


#preveri točnost - GI
pedC <- read.table("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/PedigreeAndSolutions.txt", header=TRUE)
pedGS <- read.table("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/Ped_Gen1.txt", header=TRUE)
pedSLO <- read.table("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/Ped_GenSLO1.txt", header=TRUE)
pedGSC <- read.table("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/Ped_OtherCowsGen1.txt", header=TRUE)
pedGSBD <- read.table("~/Documents/PhD/Projects/inProgress/GenomicStrategies_SireUse/Results/Ped_BmGen1.txt", header=TRUE)
#age of sires of cakajoci biki
Cak <- ped[ped$cat=="cak", c("Indiv", "Generation", "Father")]
colnames(Cak) <- c("Male", "MaleBirth", "Sire")

CakSire <- ped[ped$Indiv %in% Cak$Sire, c("Indiv", "Generation")]
colnames(CakSire) <- c("Sire", "SireBirth")

CAKSIRE <- merge(Cak, CakSire, by="Sire")
CAKSIRE$SireAge <- CAKSIRE$MaleBirth - CAKSIRE$SireBirth
CAKSIRE[order(CAKSIRE$MaleBirth),]

#najdi stare očete
ped[ped$Indiv %in% c(388722, 388722),]
#najdi mlade očete čakajočih iste generacije
ped[ped$Indiv ==431974,]
#vsi očetje te generacije
ped[ped$cat=="pb" & ped$Generation %in% c(46:50),]

#očetje generacije 55
sires <- CAKSIRE$Sire[CAKSIRE$MaleBirth==55]
pedS <- ped[ped$Indiv %in% sires,]
pedS[,c("Indiv", "gvNormUnres1", "X54")]

#primerjaj najboljšega v njegovem letu
ped[ped$Generation==45 & ped$cat=="pb",]
ped[ped$Generation==45 & ped$cat=="pb",]
ped[ped$Generation==45 & ped$cat=="pb",]


agg <- aggregate(ped$gvNormUnres1 ~ ped$Generation + ped$cat, FUN="mean")
aggS <- agg[agg$`ped$cat`=="pb",]
colnames(aggS) <- c("Gen", "Cat", "TGV")
ggplot(data=aggS, aes(x=Gen, y=TGV)) + geom_path()
pedC$Indiv <- as.factor(pedC$Indiv)
pedC$Generation <- as.factor(pedC$Generation)
pedGS$Indiv <- as.factor(pedGS$Indiv)
pedGS$Generation <- as.factor(pedGS$Generation)
pedSLO$Indiv <- as.factor(pedSLO$Indiv)
pedSLO$Generation <- as.factor(pedSLO$Generation)
pedGSC$Indiv <- as.factor(pedGSC$Indiv)
pedGSC$Generation <- as.factor(pedGSC$Generation)
pedGSBD$Indiv <- as.factor(pedGSBD$Indiv)
pedGSBD$Generation <- as.factor(pedGSBD$Generation)
PT <- ggplot(data=pedC[pedC$cat=="pb"  & pedC$Generation %in% 35:60,], aes(x=Indiv, y=gvNormUnres1, group=Generation, fill=Generation)) + geom_bar(stat="identity") + ggtitle("PT") + ylim(0, 20) + theme(legend.position = "none")
GS <- ggplot(data=pedGS[pedGS$cat=="gpb" & pedGS$Generation %in% 35:60,], aes(x=Indiv, y=gvNormUnres1, group=Generation, fill=Generation)) + geom_bar(stat="identity") + ggtitle("GS") +theme(legend.position = "none")
ggplot(data=pedGS[pedGS$cat=="pb",], aes(x=Indiv, y=gvNormUnres1, group=Generation, fill=Generation)) + geom_bar(stat="identity")
SLO <- ggplot(data=pedSLO[pedSLO$cat=="pb" &  pedSLO$Generation %in% 35:60,], aes(x=Indiv, y=gvNormUnres1, group=Generation, fill=Generation)) + geom_bar(stat="identity") + ggtitle("GS-PS") + ylim(0,20) + theme(legend.position = "none")
GSC <- ggplot(data=pedGSC[(pedGSC$cat=="pb" | pedGSC$cat=="gpb")  &  pedGSC$Generation %in% 35:60,], aes(x=Indiv, y=gvNormUnres1, group=Generation, fill=Generation)) + geom_bar(stat="identity") + ggtitle("GS-C") + ylim(0,20)+ theme(legend.position = "none")
GSBD <- ggplot(data=pedGSBD[pedGSBD$cat=="pb" &  pedGSBD$Generation %in% 35:60,], aes(x=Indiv, y=gvNormUnres1, group=Generation, fill=Generation)) + geom_bar(stat="identity") + ggtitle("GS-BD") + ylim(0,20) + theme(legend.position = "none")

library(Rmisc)
multiplot(PT, SLO, GSC, GSBD, GS)

meansd <- data.frame("Scenario" = NA, "MeanSD" = NA)
sd <- data.frame("Scenario" = NA, "SD" = NA)

pedc_PB <- pedC[pedC$cat=="pb"  & pedC$Generation %in% 35:60,]
meansd <- rbind(meansd, c("PT", mean(aggregate(pedc_PB$gvNormUnres1 ~ pedc_PB$Generation, FUN="sd")[,2])))
sd <- rbind(sd, c("PT", sd(pedc_PB$gvNormUnres1)))

pedslo_PB <- pedSLO[pedSLO$cat=="pb"  & pedSLO$Generation %in% 35:60,]
meansd <- rbind(meansd, c("GS-PS", mean(aggregate(pedslo_PB$gvNormUnres1 ~ pedslo_PB$Generation, FUN="sd")[,2])))
sd <- rbind(sd, c("GS-PS", sd(pedslo_PB$gvNormUnres1)))


pedgsc_PB <- pedGSC[(pedGSC$cat=="pb" | pedGSC$cat=="gpb")  & pedGSC$Generation %in% 35:60,]
meansd <- rbind(meansd, c("GS-C", mean(aggregate(pedgsc_PB$gvNormUnres1 ~ pedgsc_PB$Generation, FUN="sd")[,2])))
sd <- rbind(sd, c("GS-C", sd(pedgsc_PB$gvNormUnres1)))

pedgsbd_PB <- pedGSBD[pedGSBD$cat=="pb"  & pedGSBD$Generation %in% 35:60,]
meansd <- rbind(meansd, c("GS-BD", mean(aggregate(pedgsbd_PB$gvNormUnres1 ~ pedgsbd_PB$Generation, FUN="sd")[,2])))
sd <- rbind(sd, c("GS-BD", sd(pedgsbd_PB$gvNormUnres1)))

pedgs_PB <- pedGS[(pedGS$cat=="pb" | pedGS$cat=="gpb")  & pedGS$Generation %in% 35:60,]
pedgs_PB <- pedGS[pedGS$cat=="gpb"  & pedGS$Generation %in% 35:60,]
meansd <- rbind(meansd, c("GS", mean(aggregate(pedgs_PB$gvNormUnres1 ~ pedgs_PB$Generation, FUN="sd")[,2])))
sd <- rbind(sd, c("GS", sd(pedgs_PB$gvNormUnres1)))

sd <- sd[-1,]
meansd <-  meansd[-1,]

#točnost za bika 388722 in 406064 (star in mlad oče v generaciji 55)
star <- ped[ped$Indiv == 388722,]
starSol <- as.data.frame(t(star[,c(paste0("X", 40:60))]))
starSol$Gen <- 40:60
starSol$Age <- "star"
colnames(starSol)[1] <- "EBV"


mlad <- ped[ped$Indiv == 406064,]
table(ped$Generation[ped$Father==406064], ped$sex[ped$Father==406064])
table(ped$Generation[ped$Father==388722], ped$sex[ped$Father==388722])
mladSol <- as.data.frame(t(mlad[,c(paste0("X", 40:60))]))
mladSol$Gen <- 40:60
mladSol$Age <- "mlad"
colnames(mladSol)[1] <- "EBV"

S <- rbind(mladSol[,c("Gen", "Age", "EBV")], starSol[,c("Gen", "Age", "EBV")])
ggplot(data=S, aes(x=Gen, y=EBV, group=Age, colour=Age)) + geom_line() + geom_hline(yintercept = mlad$gvNormUnres1, color="red")+ geom_hline(yintercept = star$gvNormUnres1, color="blue")
