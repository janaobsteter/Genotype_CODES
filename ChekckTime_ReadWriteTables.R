library(AlphaSimR)
#do a test
#write genotypes to hdf5 file and loop
#loop through a regular table

load("~/Documents/PhD/Projects/inProgress/GenomicAlphaPart/FounderPopObject.RData")
founderPop@nChr


#geno <- pullQtlGeno(founderPop)

SP = SimParam$new(founderPop)
VarA = matrix(data = c(1.0, 0, 0, 1.0), nrow = 2); cov2cor(VarA)
VarE = matrix(data = c(3.0, 0, 0, 3.0), nrow = 2); cov2cor(VarE)
VarP = VarA + VarE; diag(VarA) / diag(VarP)
SP$addTraitA(nQtlPerChr = 1000, mean = c(1, 1), var = diag(VarA), cor = cov2cor(VarA))
SP$setGender(gender = "yes_sys")

Base = newPop(founderPop)
Base = setPheno(Base, varE = VarE)

Dams = Base[Base@gender == "F"]
Sires = Base[Base@gender == "M"]

ped <- data.frame()
BVs <- data.frame()
geno <- data.frame()
for (gen in 1:10) {
  SelCand = randCross2(females = Dams, males = Sires, nCrosses = 1000, nProgeny = 3)
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




bigGeno <- rbind(geno, geno)
bigGeno <- do.call(cbind, list(geno, geno, geno, geno))
bigGeno[1:10, 1:10]
colnames(bigGeno) <- paste0("QTL", 1:ncol(bigGeno))

dim(bigGeno)

library(rhdf5)
h5createFile("gPart.h5")
h5createGroup("gPart.h5", "genotypes")
h5delete("gPart.h5", "genotypes/geno_chunked")
h5createDataset(file = "gPart.h5", dataset = "genotypes/geno_chunked",
                dims = dim(geno), storage.mode = "integer",
                chunk = c(30, 100), level = 6)
h5save(as.matrix(geno), file = "gPart.h5", name = "genotypes/geno_chunked")
h5createDataset(file = "gPart.h5", dataset = "genotypes/geno_chunked1",
                dims = dim(geno), storage.mode = "integer",
                chunk = c(30, 100), level = 6)
h5save(as.matrix(geno), file = "gPart.h5", name = "genotypes/geno_chunked")
h5closeAll()
ls 


h5ls(file = "gPart.h5")
h5dump(file = "gPart.h5")

#test a for loop
mean1 <- ""
system.time(

  for (qtn in 1:1000) {
    mean1 <- c(mean1, mean(geno[,qtn]))
  }
)

mean2 <- ""
system.time(

  for (qtn in 1:1000) {
    mean2 <- c(mean2, mean(h5read(file = "gPart.h5", name = "genotypes/geno_chunked", 
                                index = list(NULL, qtn))))
  }
)

library(readr)
library(data.table)

processTable <- function() {
  write.table(bigGeno, "GenotypesTest.txt", quote=FALSE, row.names=FALSE)
  Geno <- read.table("GenotypesTest.txt", header=TRUE, sep= " ")
  (dim(Geno))
  mean1 <- c()
  genoPart <- bigGeno[, 700:1000]
  for (qtn in 1:300) {
    mean1 <- c(mean1, mean(genoPart[,qtn]))
  }
  system("rm GenotypesTest.txt")
  
}

processTable_readr <- function() {
  write_delim(bigGeno, "GenotypesTest.txt", delim=" ", col_names=TRUE)
  #Geno <- read_table2("GenotypesTest.txt", col_names=TRUE, col_types = cols(.default = "i"))
  Geno <- read_table2("GenotypesTest.txt", col_names=TRUE)
  (dim(Geno))
  mean1 <- c()
  genoPart <- bigGeno[, 700:1000]
  for (qtn in 1:300) {
    mean1 <- c(mean1, mean(genoPart[,qtn]))
  }
  system("rm GenotypesTest.txt")
  
}

processTable_data.table <- function() {
  fwrite(bigGeno, "GenotypesTest.txt", quote=FALSE, row.names=FALSE, sep=" ")
  Geno <- fread("GenotypesTest.txt", header=TRUE, sep=" ", select=700:1000)
  (dim(Geno))
  mean1 <- c()
  genoPart <- bigGeno[, 700:1000]
  for (qtn in 1:300) {
    mean1 <- c(mean1, mean(genoPart[[qtn]]))
  }
  system("rm GenotypesTest.txt")
  
}


processTable_hdf5 <- function() {
  h5createFile("gPart.h5")
  h5createGroup("gPart.h5", "genotypes")
  h5createDataset(file = "gPart.h5", dataset = "genotypes/geno_chunked",
                  dims = dim(bigGeno), storage.mode = "integer",
                  chunk = c(3000, 1000), level = 6)

  h5save(as.matrix(bigGeno), file = "gPart.h5", name = "genotypes/geno_chunked")
  h5closeAll()
  genoPart <- h5read(file = "gPart.h5", name = "genotypes/geno_chunked", 
                     index = list(NULL, 700:1000))
  print(dim(genoPart))  
  mean2 <- c()
  for (qtn in 1:300) {
    mean2 <- c(mean2, mean(genoPart[,qtn]))
  }
  system("rm gPart.h5")
}


processTable_feather <- function() {
  write_feather(bigGeno, "GenotypesTest.txt")
  Geno <- read_feather("GenotypesTest.txt", columns = 700:1000)
  (dim(Geno))
  mean1 <- c()
  for (qtn in 1:300) {
    mean1 <- c(mean1, mean(genoPart[[qtn]]))
  }
  system("rm GenotypesTest.txt")
}

print("write.table / read.table")
system.time(
  processTable()
)

print("readr - write_delim / read_table2")
system.time(
  processTable_readr()
)


print("data.table - fread / fwrite")
system.time(
  processTable_data.table()
)

print("hdf5")
system.time(
  processTable_hdf5()
)

print("father")
system.time(
  processTable_feather()
)




all(mean1 == mean2)
    
testG <- h5read(file = "gPart.h5", name = "genotypes/geno_chunked", 
       index = list(NULL, 1:100))
dim(testG)



## add BVs

processTable <- function() {
  write.table(bigGeno, "GenotypesTest.txt", quote=FALSE, row.names=FALSE)
  write.table(BVs, "BreedingValues.txt", quote=FALSE, row.names=FALSE)
  Geno <- read.table("GenotypesTest.txt", header=TRUE, sep= " ")
  BrVs <- read.table("BreedingValues.txt", header=TRUE, sep= " ")
  (dim(Geno))
  (dim(BVs))
  mean1 <- c()
  genoPart <- bigGeno[, 700:1000]
  for (qtn in 1:300) {
    mean1 <- c(mean1, mean(genoPart[,qtn]))
  }
  mean2 <- c()
  for (row in 1:nrow(BrVs)) {
    mean2 <- c(mean2, rowMeans(BrVs[row,2:3]))
  }
  system("rm GenotypesTest.txt")
  system("rm BreedingValues.txt")
  
}

processTable_readr <- function() {
  write_delim(bigGeno, "GenotypesTest.txt", delim=" ", col_names=TRUE)
  write_delim(BVs, "BreedingValues.txt", delim=" ", col_names=TRUE)
  #Geno <- read_table2("GenotypesTest.txt", col_names=TRUE, col_types = cols(.default = "i"))
  Geno <- read_table2("GenotypesTest.txt", col_names=TRUE)
  BrVs <- read_table2("BreedingValues.txt", col_names=TRUE)
  (dim(Geno))
  (dim(BrVs))
  mean1 <- c()
  genoPart <- bigGeno[, 700:1000]
  for (qtn in 1:300) {
    mean1 <- c(mean1, mean(genoPart[,qtn]))
  }
  mean2 <- c()
  for (row in 1:nrow(BrVs)) {
    mean2 <- c(mean2, rowMeans(BrVs[row,2:3]))
  }
  system("rm GenotypesTest.txt")
  system("rm BreedingValues.txt")
  
}

processTable_data.table <- function() {
  fwrite(bigGeno, "GenotypesTest.txt", quote=FALSE, row.names=FALSE, sep=" ")
  fwrite(bigGeno, "BreedingValues.txt", quote=FALSE, row.names=FALSE, sep=" ")
  Geno <- fread("GenotypesTest.txt", header=TRUE, sep=" ", select=700:1000)
  BrVs <- fread("BreedingValues.txt", header=TRUE, sep=" ")
  (dim(Geno))
  (dim(BrVs))
  mean1 <- c()
  genoPart <- bigGeno[, 700:1000]
  for (qtn in 1:300) {
    mean1 <- c(mean1, mean(genoPart[,qtn]))
  }
  mean2 <- c()
  for (row in 1:nrow(BrVs)) {
    mean2 <- c(mean2, rowMeans(BrVs[row,2:3]))
  }
  system("rm GenotypesTest.txt")
  system("rm BreedingValues.txt")
  
}


processTable_hdf5 <- function() {
  h5createFile("gPart.h5")
  h5createGroup("gPart.h5", "genotypes")
  h5createDataset(file = "gPart.h5", dataset = "genotypes/geno_chunked",
                  dims = dim(bigGeno), storage.mode = "integer",
                  chunk = c(3000, 1000), level = 6)
 
  h5save(as.matrix(bigGeno), file = "gPart.h5", name = "genotypes/geno_chunked")
  
  h5createGroup("gPart.h5", "breedingvalues")
  h5createDataset(file = "gPart.h5", dataset = "breedingvalues/bvs_chunked",
                  dims = dim(BVs), storage.mode = "double",
                  level = 6)
 
  h5save(as.matrix(BVs), file = "gPart.h5", name = "breedingvalues/bvs_chunked")
  h5closeAll()
  genoPart <- h5read(file = "gPart.h5", name = "genotypes/geno_chunked", 
                     index = list(NULL, 700:1000))
  BrVs <- h5read(file = "gPart.h5", name = "breedingvalues/bvs_chunked")
  print(dim(genoPart))
  print(dim(BrVs))
  mean1 <- c()
  for (qtn in 1:300) {
    mean1 <- c(mean1, mean(genoPart[,qtn]))
  }
  mean2 <- c()
  for (row in 1:nrow(BrVs)) {
    mean2 <- c(mean2, mean(BrVs[row,2:3]))
  }
  system("rm gPart.h5")
}


processTable_feather <- function() {
  write_feather(bigGeno, "GenotypesTest.txt")
  write_feather(BVs, "BreedingValues.txt")
  Geno <- read_feather("GenotypesTest.txt", columns = 700:1000)
  BrVs <- read_feather("BreedingValues.txt")
  (dim(Geno))
  (dim(BrVs))
  mean1 <- c()
  for (qtn in 1:300) {
    mean1 <- c(mean1, mean(genoPart[,qtn]))
  }
  mean2 <- c()
  for (row in 1:nrow(BrVs)) {
    mean2 <- c(mean2, rowMeans(BrVs[row,2:3]))
  }
  system("rm GenotypesTest.txt")
  system("rm BreedingValues.txt")
}

print("write.table / read.table")
system.time(
  processTable()
)

print("readr - write_delim / read_table2")
system.time(
  processTable_readr()
)


print("data.table - fread / fwrite")
system.time(
  processTable_data.table()
)

print("hdf5")
system.time(
  processTable_hdf5()
)

print("feather")
system.time(
  processTable_feather()
)
