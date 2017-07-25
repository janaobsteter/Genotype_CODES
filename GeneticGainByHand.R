#check genotypes individuals
genind <- read.table("~/bin/AlphaSim1.05Linux//IndForGeno.txt")
colnames(genind) <- "Indiv"
ped = read.table('~/bin/AlphaSim1.05Linux///SimulatedData/PedigreeAndGeneticValues_cat.txt', header=T)

genG <- merge(genind, ped, by="Indiv", all.x=T)
table(genG$Generation)
table(genG$Generation, genG$cat)

###############
#GenGain - zračunaj na roke
GIClass <- read.csv("/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Class/GenInts.txt", sep=" ")
GIClass <- GIClass[GIClass$Gen %in% 45:60,]
GIClass$Path <- paste0(GIClass$line, GIClass$sex)
GIC_mean <- aggregate(GIClass$genInt ~ GIClass$Path, FUN=mean)

GICAcc <- read.csv("/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Class/AccuraciesEBVPA.csv", sep=",")
GIC_Accmean <- aggregate(GICAcc$corEBV ~ GICAcc$Cat, FUN=mean)

GIC_mean$p <- c(0.9, 0.015,0.15, 1)
GIC_mean$i <- dnorm(qnorm(GIC_mean$p, lower.tail=FALSE)) / GIC_mean$p
GIC_mean$a <- c(0.42, 0.42, 0.90, 0.7)
sdG <- 1
GainC <- sum(GIC_mean$i * GIC_mean$a) * sdG / sum(GIC_mean$`GIClass$genInt`)

GainS <- read.table("/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Class/SimulatedData/PedigreeAndGeneticValues.txt", header=T)
GainY <- aggregate(GainS$gvNormUnres1 ~ GainS$Generation, FUN=mean)
DiffC <- (GainY$`GainS$gvNormUnres1`[60] - GainY$`GainS$gvNormUnres1`[41]) / 20

#GenGain - zračunaj na roke - genomska
GIClass <- read.csv("/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenSplosnaPop/GenInts.txt", sep=" ")
GIClass <- GIClass[GIClass$Gen %in% 45:60,]
GIClass$Path <- paste0(GIClass$line, GIClass$sex)
GIC_mean <- aggregate(GIClass$genInt ~ GIClass$Path, FUN=mean)

GICAcc <- read.csv("/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenSplosnaPop/AccuraciesEBVPA.csv", sep=",")
GIC_Accmean <- aggregate(GICAcc$corEBV ~ GICAcc$Cat, FUN=mean)

GIC_mean$p <- c(0.9, 0.015,0.15, 0.85)
GIC_mean$i <- dnorm(qnorm(GIC_mean$p, lower.tail=FALSE)) / GIC_mean$p
GIC_mean$a <- c(0.42, 0.42, 0.90, 0.6)
sdG <- 1
GainC <- sum(GIC_mean$i * GIC_mean$a) * sdG / sum(GIC_mean$`GIClass$genInt`)

GainS <- read.table("/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenSplosnaPop/SimulatedData/PedigreeAndGeneticValues.txt", header=T)
GainY <- aggregate(GainS$gvNormUnres1 ~ GainS$Generation, FUN=mean)
DiffGS <- (GainY$`GainS$gvNormUnres1`[60] - GainY$`GainS$gvNormUnres1`[41]) / 20

#GenGain - zračunaj na roke - genomska
GIClass <- read.csv("/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenSLO_BmGen//GenInts.txt", sep=" ")
GIClass <- GIClass[GIClass$Gen %in% 45:60,]
GIClass$Path <- paste0(GIClass$line, GIClass$sex)
GIC_mean <- aggregate(GIClass$genInt ~ GIClass$Path, FUN=mean)

GICAcc <- read.csv("/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenSLO_BmGen//AccuraciesEBVPA.csv", sep=",")
GIC_Accmean <- aggregate(GICAcc$corEBV ~ GICAcc$Cat, FUN=mean)

GIC_mean$p <- c(0.9, 0.015,0.15, 0.85)
GIC_mean$i <- dnorm(qnorm(GIC_mean$p, lower.tail=FALSE)) / GIC_mean$p
GIC_mean$a <- c(0.42, 0.42, 0.60, 0.90)
sdG <- 1
GainC <- sum(GIC_mean$i * GIC_mean$a) * sdG / sum(GIC_mean$`GIClass$genInt`)

GainS <- read.table("/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenSLO_BmGen//SimulatedData/PedigreeAndGeneticValues.txt", header=T)
GainY <- aggregate(GainS$gvNormUnres1 ~ GainS$Generation, FUN=mean)
DiffGen1 <- (GainY$`GainS$gvNormUnres1`[60] - GainY$`GainS$gvNormUnres1`[41]) / 20
