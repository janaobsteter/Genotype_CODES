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
herdAinv <- data.frame("Herd1" = NA, "Herd2" = NA, "SumAinv" = NA)
herdsqA <- data.frame("Herd1" = NA, "Herd2" = NA, "SumAsq" = NA)
herdEvaA <- data.frame("Herd1" = NA,  "SumA" = NA)
herdEvaAinv <- data.frame("Herd1" = NA,  "SumAinv" = NA)


library(pedigreemm)
pedC <- pedigree(sire = ped$Father, dam = ped$Mother, label=ped$Indiv)
getAInv(pedC)

makeAinv(ped[,c(2,3,4)])
Ai <- read.table("Ainv.txt")
nInd <- nrow(ped)
Ainv <- matrix(0,nrow = nInd,ncol = nInd)
# Ainv[as.matrix(Ai[,1:2])] <- Ai[,3]
dd <- diag(Ainv)
Ainv <- Ainv + t(Ainv)
diag(Ainv) <- dd

for (herd1 in herdNo) {
  for (herd2 in herdNo) {
    #to je med dvema čredama
    herd1Ind <- herds$Indiv[herds$cluster %in% c(herd1)]
    herd2Ind <- herds$Indiv[herds$cluster %in% c(herd2)]
    herdInd <- herds$Indiv[herds$cluster %in% c(herd1, herd2, pb)]
    makeA(ped[,c(2,3,4)], which=c(ped$Indiv %in% herdInd)) #to so unique elementi plus vsak sam s sabo
    A <- read.table("A.txt")
    herdA <- rbind(herdA, c(herd1, herd2, sum(A$V3)))  
    
    AinvH <- Ainv[(Ainv$V1 %in% c(herdInd,pb)) & (Ainv$V2 %in% c(herdInd, pb)),]
    herdAinv <- rbind(herdAinv, c(herd1, herd2, sum(Ainv$V3)))  
    
    #to je med dvema čredama na kvadrat
    herdsqA <- rbind(herdsqA, c(herd1, herd2, sum(A$V3^2)))    
    
    #to je med čredno in napovedno populacijo
    refEvaInd <- herds$Indiv[herds$cluster %in% c(herd1, eva)]
    makeA(pedT[,c(2,3,4)], which=c(pedT$Indiv %in% refEvaInd))
    A <- read.table("A.txt")
    herdEvaA <- rbind(herdEvaA, c(herd1, sum(A$V3)))     
    makeAinv(pedT[,c(2,3,4)], which=c(pedT$Indiv %in% refEvaInd))
    Ainv <- read.table("Ainv.txt")
    herdEvaAinv <- rbind(herdEvaAinv, c(herd1, sum(Ainv$V3)))    
    }
}







#pedC <- read.table("/home/jana/Documents/PhD/CompBio/TestingGBLUP/Ped", sep=" ", header=TRUE)
dataG <- ped$Indiv %in% gen$Indiv
dataGI <- ped[ped$cat %in% c("k", "potomciNP", "nr"), "Indiv"]
#trimPed - tukaj izubereš, koliko generacij nazaj
pedT <- ped[trimPed(ped[,c(2,3,4)], data=dataG, ngenback = 5),]
#write.table(pedT[,c(2,3,4)], "/home/jana/Documents/PhD/CompBio/TestingGBLUP/Pedigree.txt", col.names = FALSE, sep=",", row.names=FALSE, quote=FALSE)

makeA(pedT[,c(2,3,4)], which=c(pedT$Indiv %in% gen$Indiv))
A <- read.table("A.txt")

makeAinv(pedT[,c(2,3,4)])
Ai <- read.table("Ainv.txt")
nInd <- nrow(pedT)
Ainv <- matrix(0,nrow = nInd,ncol = nInd)
Ainv[as.matrix(Ai[,1:2])] <- Ai[,3]
dd <- diag(Ainv)
Ainv <- Ainv + t(Ainv)
diag(Ainv) <- dd

colnames(A) <- c("Indiv1", "Indiv2", "A")
colnames(Ai) <- c("Indiv1", "Indiv2", "Ai")

herds1 <- herds[,c("Indiv","cluster")]
colnames(herds1)[1] <- "Indiv1"
herds2 <- herds[,c("Indiv","cluster")]
colnames(herds2)[1] <- "Indiv2"

herdsA <- merge(A, herds1, by="Indiv1")
herdsA <- merge(herdsA, herds1, by="Indiv1")


b <- A[0:10,]

library(dplyr)
b$V1 <- as.factor(b$V1)
b$V2 <- as.factor(b$V2)
c <- combSummarise (b, var=c("V1", "V2"),  summarise="sum(V3)")




library(magrittr)
myData <- tbl_df(data.frame( var1 = rnorm(100), 
                             var2 = letters[1:3] %>%
                               sample(100, replace = TRUE) %>%
                               factor(), 
                             var3 = LETTERS[1:3] %>%
                               sample(100, replace = TRUE) %>%
                               factor(), 
                             var4 = month.abb[1:3] %>%
                               sample(100, replace = TRUE) %>%
                               factor()))

combSummarise (myData, var=c("var2", "var4"),  summarise=c("length(var1)"))
