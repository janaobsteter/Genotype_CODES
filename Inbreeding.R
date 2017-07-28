
library(pedigreemm)
F <- data.frame(Generation=41:60)
#################################################
scenario = "Class"
dir = paste0('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_', scenario, '/SimulatedData/')
ped <- read.table(paste0(dir,'/PedigreeAndGeneticValues_cat.txt'), header=T)[2:4]
ped$Father[ped$Father == 0] <- NA
ped$Mother[ped$Mother == 0] <- NA


pedE <- editPed(sire = ped$Father, dam = ped$Mother, label=ped$Indiv)
pede<- with(pedE, pedigree(label=label, sire=sire, dam=dam))

ped$Inbreeding <- inbreeding(pede)
genF <- aggregate(ped$Inbreeding ~ped$Generation, FUN=mean)
colnames(genF) <- c("Generation", paste0("Inbreeding_", scenario))
F <- merge(F, genF, by="Generation")

write.table(ped, paste0(dir, '/PedigreeAndGeneticValues_CatInbreeding.txt'), quote=F, row.names=F)
#plot(genF$`pedO$Inbreeding` ~ genF$`pedO$Generation`, type='line')

######################################
scenario = "GenSLO"
dir = paste0('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_', scenario, '/SimulatedData/')
ped <- read.table(paste0(dir,'/PedigreeAndGeneticValues_cat.txt'), header=T)[2:4]
ped$Father[ped$Father == 0] <- NA
ped$Mother[ped$Mother == 0] <- NA


pedE <- editPed(sire = ped$Father, dam = ped$Mother, label=ped$Indiv)
pede<- with(pedE, pedigree(label=label, sire=sire, dam=dam))

ped$Inbreeding <- inbreeding(pede)
genF <- aggregate(ped$Inbreeding ~ped$Generation, FUN=mean)
colnames(genF) <- c("Generation", paste0("Inbreeding_", scenario))
F <- merge(F, genF, by="Generation")

write.table(ped, paste0(dir, '/PedigreeAndGeneticValues_CatInbreeding.txt'), quote=F, row.names=F)
#plot(genF$`pedO$Inbreeding` ~ genF$`pedO$Generation`, type='line')
######################################
scenario = "Gen_SplosnaPop"
dir = paste0('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_', scenario, '/SimulatedData/')
ped <- read.table(paste0(dir,'/PedigreeAndGeneticValues_cat.txt'), header=T)[2:4]
ped$Father[ped$Father == 0] <- NA
ped$Mother[ped$Mother == 0] <- NA


pedE <- editPed(sire = ped$Father, dam = ped$Mother, label=ped$Indiv)
pede<- with(pedE, pedigree(label=label, sire=sire, dam=dam))

ped$Inbreeding <- inbreeding(pede)
genF <- aggregate(ped$Inbreeding ~ped$Generation, FUN=mean)
colnames(genF) <- c("Generation", paste0("Inbreeding_", scenario))
F <- merge(F, genF, by="Generation")

write.table(ped, paste0(dir, '/PedigreeAndGeneticValues_CatInbreeding.txt'), quote=F, row.names=F)
#plot(genF$`pedO$Inbreeding` ~ genF$`pedO$Generation`, type='line')
######################################
scenario = "GenSLO_BmGen"
dir = paste0('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_', scenario, '/SimulatedData/')
ped <- read.table(paste0(dir,'/PedigreeAndGeneticValues_cat.txt'), header=T)[2:4]
ped$Father[ped$Father == 0] <- NA
ped$Mother[ped$Mother == 0] <- NA


pedE <- editPed(sire = ped$Father, dam = ped$Mother, label=ped$Indiv)
pede<- with(pedE, pedigree(label=label, sire=sire, dam=dam))

ped$Inbreeding <- inbreeding(pede)
genF <- aggregate(ped$Inbreeding ~ped$Generation, FUN=mean)
colnames(genF) <- c("Generation", paste0("Inbreeding_", scenario))
F <- merge(F, genF, by="Generation")

write.table(ped, paste0(dir, '/PedigreeAndGeneticValues_CatInbreeding.txt'), quote=F, row.names=F)
#plot(genF$`pedO$Inbreeding` ~ genF$`pedO$Generation`, type='line')
######################################
scenario = "Gen"
dir = paste0('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_', scenario, '/SimulatedData/')
dir = paste0('/home/jana/bin/AlphaSim1.05Linux/SimulatedData/')
ped <- read.table(paste0(dir,'/PedigreeAndGeneticValues_cat.txt'), header=T)[2:4]
ped$Father[ped$Father == 0] <- NA
ped$Mother[ped$Mother == 0] <- NA


pedE <- editPed(sire = ped$Father, dam = ped$Mother, label=ped$Indiv)
pede<- with(pedE, pedigree(label=label, sire=sire, dam=dam))

ped$Inbreeding <- inbreeding(pede)
genF <- aggregate(ped$Inbreeding ~ped$Generation, FUN=mean)
colnames(genF) <- c("Generation", paste0("Inbreeding_", scenario))
F <- merge(F, genF, by="Generation")

write.table(ped, paste0(dir, '/PedigreeAndGeneticValues_CatInbreeding.txt'), quote=F, row.names=F)
#plot(genF$`pedO$Inbreeding` ~ genF$`pedO$Generation`, type='line')

write.csv(F, "/home/jana/bin/AlphaSim1.05Linux/Inbreeding_Scenarios.csv", quote=F, row.names=F)


