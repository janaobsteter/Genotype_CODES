
library(pedigreemm)
F <- data.frame(Generation=41:60)
#################################################
for (scenario in c("Class", "GenSLO", "GenSplosnaPop", "GenSLO_BmGen", "Gen")) {
    dir = paste0('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_', scenario, '/SimulatedData/')
    pedO <- read.table(paste0(dir,'/PedigreeAndGeneticValues_cat.txt'), header=T)
    ped <- pedO[2:4,]
    ped$Father[ped$Father == 0] <- NA
    ped$Mother[ped$Mother == 0] <- NA
    
    
    pedE <- editPed(sire = ped$Father, dam = ped$Mother, label=ped$Indiv)
    pede<- with(pedE, pedigree(label=label, sire=sire, dam=dam))
    
    pedO$Inbreeding <- inbreeding(pede)
    genF <- aggregate(pedO$Inbreeding ~pedO$Generation, FUN=mean)
    colnames(genF) <- c("Generation", paste0("Inbreeding_", scenario))
    F <- merge(F, genF, by="Generation")
    
    write.table(pedO, paste0(dir, '/PedigreeAndGeneticValues_CatInbreeding.txt'), quote=F, row.names=F)
}

write.csv(F, "/home/jana/bin/AlphaSim1.05Linux/Inbreeding_Scenarios.csv", quote=F, row.names=F)

#plot(genF$`pedO$Inbreeding` ~ genF$`pedO$Generation`, type='line')
