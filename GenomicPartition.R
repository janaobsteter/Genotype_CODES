library(AlphaSimR)
nMales = 100
nFemales = 10000


nGenerationBurn = 20
nGenerationEval = 20

# ---- Base population genomes ----

founderPop = runMacs(nInd = nMales + nFemales,
                     nChr = 1,
                     segSites = 1000,
                     nThreads = 4,
                     # species = "GENERIC")
                     species = "CATTLE")

#save.image(file = "~/Documents/PhD/Projects/GenomicAlphaPart//FounderPop.RData")
#save(founderPop, file = "~/Documents/PhD/Projects/GenomicAlphaPart//FounderPopObject.RData")


load("~/Documents/PhD/Projects/GenomicAlphaPart/FounderPopObject.RData")
founderPop@nChr



setwd("/home/jana///Documents/PhD/Projects/GenomicAlphaPart/")
#setwd("/home/jana/Documents/KISdir/bin/AlphaSim1.05Linux/REAL20GenSel_Gen1//")
library(reshape)
library(AlphaSimR)
SP = SimParam$new(founderPop)
VarA = matrix(data = c(1.0, 0, 0, 1.0), nrow = 2); cov2cor(VarA)
VarE = matrix(data = c(3.0, 0, 0, 3.0), nrow = 2); cov2cor(VarE)
VarP = VarA + VarE; diag(VarA) / diag(VarP)
SP$addTraitA(nQtlPerChr = 10, mean = c(1, 1), var = diag(VarA), cor = cov2cor(VarA))
SP$setGender(gender = "yes_sys")

Base = newPop(founderPop)
Base = setPheno(Base, varE = VarE)

Dams = Base[Base@gender == "F"]
Sires = Base[Base@gender == "M"]

ped <- data.frame()
BVs <- data.frame()
geno <- data.frame()
for (gen in 1:5) {
  SelCand = randCross2(females = Dams, males = Sires, nCrosses = 100, nProgeny = 3)
  SelCand = setPheno(SelCand, varE = VarE)
  
  ped <- rbind(ped, data.frame(
    IId = SelCand@id,
    FId = SelCand@father,
    MId = SelCand@mother,
    Sex = SelCand@gender,
    Generation = gen)
  )
  
  BVs <- rbind(BVs, data.frame(IId = SelCand@id,
                               TGV1 = SelCand@gv[,1],
                               TGV2 = SelCand@gv[,2]))
  
  geno <- rbind(geno, as.data.frame(pullQtlGeno(SelCand)))
  
  Dams = selectInd(SelCand, 100, trait = function(x) rowMeans(scale(x)), use = "gv", gender = "F")
  Sires = selectInd(SelCand, 100, trait = function(x) rowMeans(scale(x)), use = "gv", gender = "M")
}


######### -- Inputs -- ############
nTraits = 2
colBV = c("TGV1", "TGV2")

# BVs <- data.frame(IId = Base@id,
#                   Pheno1 = Base@pheno[,1],
#                   Pheno2 = Base@pheno[,2],
#                   TGV1 = Base@gv[,1],
#                   TGV2 = Base@gv[,2]
#                   )
# 
# geno <- pullQtlGeno(Base)
geno$IId <- as.vector(rownames(geno))
 
addEff <- data.frame(SNP = 1:length(SP$traits[[1]]@addEff),
                    Trait1 = SP$traits[[1]]@addEff,
                    Trait2 = SP$traits[[2]]@addEff)

write.table(ped, "PedigreeSimulation.txt", row.names=FALSE, quote=FALSE)
write.table(BVs, "BreedinValuesSimulation.txt", row.names=FALSE, quote=FALSE)
write.table(geno, "GenotypesSimulation.txt", row.names=FALSE, quote=FALSE)
write.table(addEff, "AdditiveEffectsSimulation.txt", row.names=FALSE, quote=FALSE)


# #read in the pedigree
# #ped <- read.table("PedPart.txt")
# # ped <- read.table("Ped.txt")
# #Ped <- read.table("PedTotalPart.txt")
# Ped <- read.table("TotalPed.txt")
# colnames(ped) <- c("IId", "FId", "MId")
# colnames(Ped) <- c("IId", "FId", "MId", "Sex", "Generation", "Pheno", "TGV")
# 
# 
# #read in the genotype file
# #geno <- read.table("GenoPartIndInfo.txt")
# geno <- read.table("Geno.txt")
# #name the columns of the genotype file
# colnames(geno) <- c("IId", paste0("SNP", 1:(ncol(geno)-1) ))
# #melt the genotype file to get the genotypes in one column
# #genoM <- melt(geno, id.vars = "IId")
# #name the SNP column
# #colnames(genoM)[2] <- "SNP"

#read in additive effects of the QTLs
addEff <- read.table("AddEffects.txt")
#name the SNPs as SNP1, SNP2 ...
addEff$V1 <- paste0("SNP", addEff$V1)
#name the columns
colnames(addEff) <- c("SNP", "Eff")
#assign a value of the a and A allele (a in AlphaPeel is major, 0 in AlphaSimR is aa)
addEff$a <- -addEff$Eff / 2
addEff$A <- -addEff$a


################--- Read in Files ---#######################
setwd("/home/jana///Documents/PhD/Projects/GenomicAlphaPart/")
ped <- read.table("PedigreeSimulation.txt", header=TRUE)
BVs <- read.table( "BreedinValuesSimulation.txt", header=TRUE)
geno <- read.table( "GenotypesSimulation.txt", header=TRUE)
addEff <- read.table("AdditiveEffectsSimulation.txt", header=TRUE)



############-- Input parameters -- ##########
nTrait = 2
colBV = c("TGV1", "TGV2")

#KEEP ONLY the last generation individuals
#genoA <- geno[geno$IId %in% Ped$IId[Ped$Generation == 58],]
#only keep the individuals with genotyped parents
genoInd <- ped$IId[(ped$IId %in% geno$IId) & (ped$FId %in% geno$IId) & (ped$MId %in% geno$IId)]
genoA <- geno[geno$IId %in% genoInd,]
############################

###########################
#compute allele substituion effects
library(reshape)
addEff$SNP <- paste0("SNP", 1:nrow(addEff))
addEffm <- melt(addEff, id.vars = "SNP")
colnames(addEffm)[colnames(addEffm) == "variable"] <- "Trait"
colnames(addEffm)[colnames(addEffm) == "value"] <- "Eff"
addEffm$a <- -addEffm$Eff / 2
addEffm$A <- -addEffm$a

#GENOTYPES
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              

genoAInd <- as.data.frame(genoA$IId)
colnames(genoA) <- ifelse(colnames(genoA) == "IId", colnames(genoA), paste0("SNP", which(colnames(genoA) == colnames(genoA))))
colnames(geno) <- ifelse(colnames(geno) == "IId", colnames(geno), paste0("SNP", which(colnames(geno) == colnames(geno))))
locusMST <- data.frame(SNP = NA, Mean = NA, Sum = NA)




#set variables
nInd = nrow(genoA)
nLoci <- length(addEff$SNP)
mstArray <- array(dim      = c(nInd, nLoci, nTrait),
                  dimnames = list(c(genoA$IId),
                                  c(addEff$SNP),
                                  c(paste0("Trait", 1:nTrait))))
dim(mstArray)
dimnames(mstArray)




###########--- DO the partitioning --- ################
for (trait in 1:nTrait) {
  addEff_trait <- addEffm[addEffm$Trait == paste0("Trait", trait),]
  for (col in colnames(genoA)[colnames(genoA) !="IId"]) {
    print(paste0("This is ", col))
    genoM <- melt(geno[,c("IId", col)],  id.vars = "IId")
    colnames(genoM)[2] <- "SNP"
    
    genoMp <- merge(genoM, ped, by = "IId", all.x = TRUE)
    fathers <- genoMp[genoMp$IId %in% ped$FId,c("IId", "SNP", "value")]
    colnames(fathers) <- c("FId", "SNP", "genoF")
    mothers <- genoMp[genoMp$IId %in% ped$MId,c("IId", "SNP", "value")]
    colnames(mothers) <- c("MId", "SNP", "genoM")
    genoMp$FId <- as.factor(genoMp$FId)
    genoMp$MId <- as.factor(genoMp$MId)
    fathers$FId <- as.factor(fathers$FId)
    mothers$MId <- as.factor(mothers$MId)
    
    genoMp <- genoMp[genoMp$IId %in% genoInd,]
    genoMP <- merge(genoMp, fathers, by = c("FId", "SNP"), all.x=TRUE)
    genoMP <- merge(genoMP, mothers, by = c("MId", "SNP"), all.x=TRUE)
  
    #merge with additive effects
    genoMP <- merge(genoMP, addEff_trait, by="SNP", all.x=TRUE)
    #set the value of the individual genotype
    genoMP$effIId <- ifelse(genoMP$value == 1, 0, ifelse(genoMP$value == 0, -genoMP$Eff, genoMP$Eff))
    #find rthe difference between AlphaSim TGV and true TGV
    # trueTGV <- aggregate(genoMP$effIId ~ genoMP$IId, FUN="sum")
    # colnames(trueTGV) <- c("IId", "trueTGV")
    # #remove
    # PED <- merge(ped, trueTGV, by="IId")
    # PED$TGV <- PED$trueTGV
    #compute PA and MST
    # fat <- c[ped$IId %in% ped$FId, c("IId", "TGV")]
    # colnames(fat) <- c("FId", "TGV_F")
    # mot <- ped[ped$IId %in% ped$MId, c("IId", "TGV")]
    # colnames(mot) <- c("MId", "TGV_M")
    # #
    # PED <- merge(PED, fat, by="FId", all.x=TRUE)
    # PED <- merge(PED, mot, by="MId", all.x=TRUE)
    # PED$PA <- (PED$TGV_F + PED$TGV_M) / 2
    # PED$MST <- PED$TGV - PED$PA
    #compute average genotype of the parents
    genoMP$genoPA <- (genoMP$genoF + genoMP$genoM) / 2
    #compute average genotype value of the parents
    genoMP$effPA <- (genoMP$genoPA - 1) * genoMP$Eff
    #order
    genoMP <- genoMP[order(genoMP$IId),]
    #compute locus MST
    genoMP$locusMST <- genoMP$effIId - genoMP$effPA
    
    locusMST <- rbind(locusMST, 
                      c(unique(genoMP$SNP), mean(genoMP$locusMST), sum(genoMP$locusMST)))
    mstArray[,,paste0("Trait", trait)][,col] <- genoMP$locusMST    
  }
}

mstArray
#inspect trait 1, SNP8
hist(mstArray[,8,1 ])
hist(mstArray[,7,1 ])

#summarize the array
#by individuals
length(apply(mstArray, c(1,3), mean))
str(apply(mstArray, c(1,3), mean))
dim(apply(mstArray, c(1,3), mean))
indByTrait <- as.data.frame(apply(mstArray, c(1,3), mean))
head(indByTrait)
qplot(x = indByTrait[,1], geom="histogram", bins = 100)
ggplot(data = indByTrait, aes(x = rownames(indByTrait),y = Trait1)) + geom_bar(stat="identity") + 
  theme(axis.text.x = element_blank()) 
#summarise by loci
LociByTrait <- as.data.frame(apply(mstArray, c(2,3), mean))
dim(LociByTrait)
LociByTrait$SNP <- rownames(LociByTrait)
LociByTrait$SNP <- factor(LociByTrait$SNP, levels = paste0("SNP", 1:nrow(LociByTrait)))
qplot(x = LociByTrait[,1], geom="histogram", bins = 100)

LociByTraitm <- melt(LociByTrait, id.vars = "SNP")
ggplot(data = LociByTraitm, aes(x = SNP, y = value)) + geom_bar(stat="identity") + theme_bw() + theme(axis.title.x = element_blank()) + 
  facet_grid(rows = vars(variable))


write.table(locusMST, "LocusMST.txt", quote=FALSE, row.names=FALSE)

library(Rmisc)
plotInd <- function(ind) {
  Pedm <- melt(Ped[,c("IId", "TGV", "PA", "MST")], id.vars = c("IId"))
  Pedm$IId <- as.factor(Pedm$IId)
  mst <- ggplot(Pedm[Pedm$IId == ind,]) + 
    geom_bar(aes(x = IId, y = value, fill=variable), stat = "identity", position = "dodge")+ 
    scale_fill_manual("", values = c("black", "skyblue2", "red4")) +
    theme(legend.position = "top")
  lociMst <- ggplot(data = genoMP[genoMP$IId == ind,], aes(x = IId, y = locusMST, fill=SNP)) + 
    geom_bar(stat = "identity", position = "stack") + ylim(c(0, unique(genoMP$TGV[genoMP$IId == ind]))) + 
    theme(legend.position = "top")
  
  lay <- rbind(c(1,2,2),
               c(1,2,2),
               c(1,2,2))
  grid.arrange(mst, lociMst, layout_matrix = lay)
}

plotInd(25)






############################
############################
#HAPLOTIPI
#merge the genotypes with the additive effects
genoA <- merge(genoM, addEff, by = "SNP")
#compute the value of the genotype
genoA$snpA <- ifelse(genoA$value == 0, -genoA$Eff, ifelse(genoA$value == 1, 0, genoA$Eff))
#aggregate additive effects by individual
indGV <- aggregate(genoA$snpA ~genoA$IId, FUN="sum")
#rename the columns
colnames(indGV) <- c("IId", "snpA")
#merge values with the generation
indGV <- merge(indGV, Ped[,c("IId", "Generation")], by="IId")
#obtain the mean of the base population!!! (mean TGV to subtract)
subtractM <- mean(indGV$snpA[indGV$Generation == 1])



#read in the haplotypes
hap <- read.table("outputs/alphapeel.haps")
#rename the columns to IId and SNP1, SNP2 ...
colnames(hap) <- c("IId", paste0("SNP", 1:(ncol(hap)-1) ))

#seg <- read.table("outputs/alphapeel.seg")
#colnames(seg) <- c("IId", paste0("SNP", 1:(ncol(seg)-1) ))

#assign the haplotypes to the lines - Geno columns
hapProb <- c("aa", "aA", "Aa", "AA")
hap$Geno <- rep(hapProb, nrow(hap)/4)
#melt the haplotypes to get one SNP per line
hapM <- melt(hap, id.vars = c("IId", "Geno"))
#create a group IId-SNP
hapM$Group <- paste0(hapM$IId, hapM$variable)
#split for each individual and SNP and obtain the most likely haplotype
haplotyped <- do.call(rbind, lapply(split(hapM,hapM$Group), function(x) {return(x[which.max(x$value),])}))
#keep only the IId, SNP and haplotype
haplotyped <- haplotyped[,1:3]
#rename the SNP column
colnames(haplotyped)[3] <- "SNP"
#merge with additive effects
haplotyped <- merge(haplotyped, addEff, by="SNP")
#order by IId
haplotyped <- haplotyped[order(haplotyped$IId),]
#assign fathers contribution - if aa or aA, a came from fathers (otherwise A)
haplotyped$F_contr <- ifelse(haplotyped$Geno %in% c("aa", "aA"), haplotyped$a, haplotyped$A)
#assign mothers contribution - if aa or Aa, a came from fathers (otherwise A)
haplotyped$M_contr <- ifelse(haplotyped$Geno %in% c("aa", "Aa"), haplotyped$a, haplotyped$A)

#this is a check - if the father and mother contribution sum to the TGV
#sum the contributions
haplotyped$SNPCont <- haplotyped$F_contr + haplotyped$M_contr
#aggregate by IId to obtain the TGV
checkM <- aggregate(haplotyped$SNPCont ~ haplotyped$IId, FUN="sum")
#rename the columns
colnames(checkM) <- c("IId", "SNPContr")
#compare to indGV (computed from additive effects combined)
checkM$SNPContr[checkM$IId == 1] == indGV$snpA[indGV$IId == 1]


#ALI POTREBNO SPUSTIT ČEZ PEDIGRE???
##########################################################
#probabilistic
##########################################################

hap <- read.table("outputs/alphapeel.haps")
colnames(hap) <- c("IId", paste0("SNP", 1:(ncol(hap)-1) ))  
hapProb <- c("aa", "aA", "Aa", "AA")
hap$Geno <- rep(hapProb, nrow(hap)/4)

hapM <- melt(hap, id.vars = c("IId", "Geno"))
HAP <- spread(hapM, Geno, value)
colnames(HAP)[2] <- "SNP"

#read in the effects
addEff <- read.table("AddEffects.txt")
addEff$V1 <- paste0("SNP", addEff$V1)
colnames(addEff) <- c("SNP", "Eff")
addEff$a <- addEff$Eff / 2
addEff$A <- -addEff$a

HAP <- merge(HAP, addEff, by="SNP", all.x=TRUE)

#read in the pedigree
ped <- read.table("Ped.txt")
colnames(ped) <- c("IId", "FId", "MId")
HAP <- merge(HAP, ped, by="IId")


HAP$F_contr <- NA
HAP$M_contr <- NA

for (row in 1:nrow(HAP)) {
  HAP$F_contr[row] <- (HAP$aa[row] + HAP$aA[row]) * HAP$a[row] + (HAP$aA[row] + HAP$AA[row]) * HAP$A[row]
  HAP$M_contr[row] <- (HAP$aa[row] + HAP$Aa[row]) * HAP$a[row] + (HAP$aA[row] + HAP$AA[row]) * HAP$A[row]
}

#test
HAP$SNPEff <- HAP$F_contr + HAP$M_contr
M <- aggregate(HAP$SNPEff ~HAP$IId, FUN="mean")
gv <- read.table("TotalPed.txt")


####################
#this is the loci sum and mean from one generatio of simulated population
locusMST <- read.table("Documents/PhD/Projects/inProgress/GenomicAlphaPart/LocusMST.txt", header=TRUE)[-1,]
locusMST$SNP <- paste0("SNP", 1:nrow(locusMST))
locusMST$SNPnumber <- 1:nrow(locusMST)
head(locusMST)
locusMST$groupSum <- ifelse(abs(locusMST$Sum) > 5, "extreme", "normal")
summary(locusMST$Mean)
locusMST$groupMean <- ifelse(abs(locusMST$Mean) > 5e-04, "extreme", "normal")


ggplot(data = locusMST, aes(x = SNPnumber, y = Mean, colour= groupMean)) + geom_bar(stat = "identity")+
  scale_colour_manual("", values = c("red", "black")) + theme(legend.position = "none") + ylab("Mean MST")
ggplot(data = locusMST, aes(x = SNPnumber, y = Sum, colour = groupSum)) + geom_bar(stat = "identity") +
  scale_colour_manual("", values = c("red", "black")) + theme(legend.position = "none") + ylab("Sum MST")
summary(locusMST$Sum)


###############################
BiocManager::install("rhdf5")
library(rhdf5)

m1 <- matrix(rep(1:20000, each = 100), ncol = 20000, byrow = FALSE)
ex_file <- tempfile(fileext = ".h5")
h5write(m1, file = ex_file, name = "counts", level = 6)
h5save(m1, file = "counts.h5")
h5dump("counts.h5")

system.time(
  res1 <- h5read(file = ex_file, name = "counts", 
                 index = list(NULL, 1:10000))
)

index <- list(NULL, seq(from = 1, to = 20000, by = 2))
system.time(
  res2 <- h5read(file = ex_file, name = "counts", 
                 index = index)
)


start <- c(1,1)
stride <- c(1,2)
block <- c(100,1)
count <- c(1,10000)
system.time(
  res3 <- h5read(file = ex_file, name = "counts", start = start,
                 stride = stride, block = block, count = count)
)



h5createDataset(file = ex_file, dataset = "counts_chunked", 
                dims = dim(m1), storage.mode = "integer", 
                chunk = c(100,100), level = 6)
h5write(obj = m1, file = ex_file, name = "counts_chunked")


library(dplyr)
set.seed(1234)
columns <- sample(x = seq_len(20000), size = 10000, replace = FALSE) %>%
  sort()

f1 <- function(cols, name) { 
  h5read(file = ex_file, name = name, 
         index = list(NULL, cols))
}
system.time(res4 <- vapply(X = columns, FUN = f1, 
                           FUN.VALUE = integer(length = 100), 
                           name = 'counts'))


A = 1:7; B = 1:18; D = seq(0,1,by=0.1)
h5save(A, B, D, file="ex_save.h5")
h5dump("ex_save.h5")
