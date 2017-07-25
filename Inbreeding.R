ped <- read.table('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Gen/SimulatedData/PedigreeAndGeneticValues.txt', header=T)[2:4]
ped$Father[pedO$Father == 0] <- NA
ped$Mother[pedO$Mother == 0] <- NA

library(pedigreemm)

pedE <- editPed(sire = ped$Father, dam = ped$Mother, label=ped$Indiv)
pede<- with(pedE, pedigree(label=label, sire=sire, dam=dam))

pedO <- read.table('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Gen/SimulatedData/PedigreeAndGeneticValues.txt', header=T)[1:4]
pedO$Inbreeding <- inbreeding(pede)
genF <- aggregate(pedO$Inbreeding ~pedO$Generation, FUN=mean)
plot(genF$`pedO$Inbreeding` ~ genF$`pedO$Generation`, type='line')
######################################
pedC <- read.table('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Class//SimulatedData/PedigreeAndGeneticValues.txt', header=T)[2:4]
pedC$Father[pedC$Father == 0] <- NA
pedC$Mother[pedC$Mother == 0] <- NA

pedCe <- editPed(sire = pedC$Father, dam = pedC$Mother, label=pedC$Indiv)
pedeC<- with(pedCe, pedigree(label=label, sire=sire, dam=dam))

pedOC <- read.table('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Class/SimulatedData/PedigreeAndGeneticValues.txt', header=T)[1:4]
pedOC$Inbreeding <- inbreeding(pedeC)
genF <- aggregate(pedOC$Inbreeding ~pedOC$Generation, FUN=mean)
lines(genF$`pedOC$Inbreeding` ~ genF$`pedOC$Generation`, col="red")


