library(pedigreemm)
F <- data.frame(Generation=1:60)
#################################################
for (scenario in c("Class1")) {
    dir = paste0('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_', scenario, '/SimulatedData/')
    pedO <- read.table(paste0(dir,'/PedigreeAndGeneticValues_cat.txt'), header=T)
    ped <- pedO[,2:4]
    ped$Father[ped$Father == 0] <- NA
    ped$Mother[ped$Mother == 0] <- NA
    
    
    pedE <- editPed(sire = ped$Father, dam = ped$Mother, label=ped$Indiv)
    pede<- with(pedE, pedigree(label=label, sire=sire, dam=dam))
    
    pedO$Inbreeding <- inbreeding(pede)
    genF <- aggregate(pedO$Inbreeding ~pedO$Generation, FUN=mean)
    colnames(genF) <- c("Generation", paste0("Inbreeding_", scenario))
    F <- merge(F, genF, by="Generation")
    
    write.table(pedO, paste0(dir, '/PedigreeAndGeneticValues_CatInbreeding.txt'))
}

write.csv(F, "/home/jana/bin/AlphaSim1.05Linux/Inbreeding_Scenarios.csv")

#plot(genF$`pedO$Inbreeding` ~ genF$`pedO$Generation`, type='line')
library(ggplot2)
library(reshape)

Ft <- read.csv("/home/jana/bin/AlphaSim1.05Linux/Inbreeding_Scenarios.csv")[,-1]
Ftmelt <- melt(Ft, id.vars = "Generation")
ggplot(data=Ftmelt, aes(x = Generation, y = value, fill = variable, colour=variable)) + geom_path() + ylab("Inbreeding") +
  scale_color_hue("Shema", labels=c("Conventional", "GenomicSLO", "GenBulls on Other Cows", "GenBulls on Bull Dams", "GenBulls on All Cows"))
