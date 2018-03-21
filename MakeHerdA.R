library(pedigree)
gen <- read.table('//home/jana/bin/AlphaSim1.05Linux/IndForGeno_5gen.txt') #tukaj so krave, pb in potomciNP
colnames(gen) <- "Indiv"

ped <- read.table('/home/jana/bin/AlphaSim1.05Linux/REALFillIn20BurnIn20/SimulatedData/PedigreeAndGeneticValues_cat.txt', sep=" ", header=TRUE)
pedCat <- merge(gen, ped, by="Indiv", all.x=TRUE)

herds <- read.table("/home/jana/Documents/PhD/CompBio/TestingGBLUP/PedCows_HERDS.txt", header=TRUE)
herds$Indiv <- as.character(herds$Indiv)

ped$Indiv <- as.character(ped$Indiv)

eva <- ped$Indiv[ped$cat =="potomciNP"]
pb <- ped$Indiv[ped$cat =="pb"]

herdNo <- unique(herds$cluster)

herdA <- data.frame("Herd1" = NA, "Herd2" = NA, "SumA" = NA)
#herdAinv <- data.frame("Herd1" = NA, "Herd2" = NA, "SumAinv" = NA)
herdsqA <- data.frame("Herd1" = NA, "Herd2" = NA, "SumAsq" = NA)
herdEvaA <- data.frame("Herd1" = NA,  "SumA" = NA)
#herdEvaAinv <- data.frame("Herd1" = NA,  "SumAinv" = NA)


library(pedigreemm)
#pedC <- pedigree(sire = ped$Father, dam = ped$Mother, label=ped$Indiv)
#getAInv(pedC)

for (herd1 in herdNo) {
  for (herd2 in herdNo) {
    #to je med dvema čredama
    herd1Ind <- herds$Indiv[herds$cluster %in% c(herd1)]
    herd2Ind <- herds$Indiv[herds$cluster %in% c(herd2)]
    herdInd <- herds$Indiv[herds$cluster %in% c(herd1, herd2, pb)]
    makeA(ped[,c(2,3,4)], which=c(ped$Indiv %in% herdInd)) #to so unique elementi plus vsak sam s sabo
    A <- read.table("A.txt")
    herdA <- rbind(herdA, c(herd1, herd2, sum(A$V3)))  
    
    #AinvH <- Ainv[(Ainv$V1 %in% c(herdInd,pb)) & (Ainv$V2 %in% c(herdInd, pb)),]
    #herdAinv <- rbind(herdAinv, c(herd1, herd2, sum(Ainv$V3)))  
    
    #to je med dvema čredama na kvadrat
    herdsqA <- rbind(herdsqA, c(herd1, herd2, sum(A$V3^2)))    
    
    #to je med čredno in napovedno populacijo
    refEvaInd <- herds$Indiv[herds$cluster %in% c(herd1, eva)]
    makeA(ped[,c(2,3,4)], which=c(ped$Indiv %in% refEvaInd))
    A <- read.table("A.txt")
    herdEvaA <- rbind(herdEvaA, c(herd1, sum(A$V3)))     
    }
}


write.table(herdA, "HerdA.txt", quote=FALSE)
write.table(herdsqA, "HerdsqaA.txt", quote=FALSE)
write.table(herdEvaA, "HerdEvaA.txt", quote=FALSE)