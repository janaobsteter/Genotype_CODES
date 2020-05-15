#load the required packages
library(readr)
library(data.table)
library(hdf5)
library(feather)

#create sample genotype and breeding values tables
geno <- as.data.frame(matrix(sample(c(0:2), 1.2e+08, replace = TRUE), nrow = 30000, ncol = 4000))
BVs <- data.frame(ID = 1:30000, Trait1 = rnorm(30000), Trait2 = rnorm(30000))

'create function to :
  1) write the geno table
  2) read the geno table back in
  3) select 300 columns (do this on the read if possible)
  4) loop through the columns

Tested package: utils, readr, data.table, hdf5 and feather

# note: you can select columns with readr as well, but it is quite tedious with a large number of selected columns
# note: there is a lot of space to tweak the parameters with hdf5 package - especially the number of chunks, which also define the compression
'
 
#use utils write.table / read.table 
processTable <- function() {
  write.table(geno, "GenotypesTest.txt", quote=FALSE, row.names=FALSE)
  Geno <- read.table("GenotypesTest.txt", header=TRUE, sep= " ")
  mean1 <- c()
  genoPart <- geno[, 700:1000]
  for (qtn in 1:300) {
    mean1 <- c(mean1, mean(genoPart[,qtn]))
  }
  system("rm GenotypesTest.txt")
  
}

#use readr write_delim / read_table2
processTable_readr <- function() {
  write_delim(geno, "GenotypesTest.txt", delim=" ", col_names=TRUE)
  #Geno <- read_table2("GenotypesTest.txt", col_names=TRUE, col_types = cols(.default = "i"))
  Geno <- read_table2("GenotypesTest.txt", col_names=TRUE)
  mean1 <- c()
  genoPart <- geno[, 700:1000]
  for (qtn in 1:300) {
    mean1 <- c(mean1, mean(genoPart[,qtn]))
  }
  system("rm GenotypesTest.txt")
  
}

#use data.table fwrite / fread
processTable_data.table <- function() {
  fwrite(geno, "GenotypesTest.txt", quote=FALSE, row.names=FALSE, sep=" ")
  Geno <- fread("GenotypesTest.txt", header=TRUE, sep=" ", select=700:1000)
  mean1 <- c()
  genoPart <- geno[, 700:1000]
  for (qtn in 1:300) {
    mean1 <- c(mean1, mean(genoPart[[qtn]]))
  }
  system("rm GenotypesTest.txt")
  
}

#use hdf5 
processTable_hdf5 <- function() {
  h5createFile("gPart.h5")
  h5createGroup("gPart.h5", "genotypes")
  h5createDataset(file = "gPart.h5", dataset = "genotypes/geno_chunked",
                  dims = dim(geno), storage.mode = "integer",
                  chunk = c(3000, 1000), level = 6)
  
  h5save(as.matrix(geno), file = "gPart.h5", name = "genotypes/geno_chunked")
  h5closeAll()
  genoPart <- h5read(file = "gPart.h5", name = "genotypes/geno_chunked", 
                     index = list(NULL, 700:1000))
  mean2 <- c()
  for (qtn in 1:300) {
    mean2 <- c(mean2, mean(genoPart[,qtn]))
  }
  system("rm gPart.h5")
}

#use feather write_feather / read_feather
processTable_feather <- function() {
  write_feather(geno, "GenotypesTest.txt")
  Geno <- read_feather("GenotypesTest.txt", columns = 700:1000)
  mean1 <- c()
  for (qtn in 1:300) {
    mean1 <- c(mean1, mean(genoPart[[qtn]]))
  }
  system("rm GenotypesTest.txt")
}


###################################
# Run and time the functions
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


####################################################################
####################################################################
# Use the same package as before, but now write the BVs table as well
'create function to :
  1) write the geno table
  2) write the BVs table
  3) read the geno table back in
  4) read the BVs table back in
  5) select 300 columns of geno table (do this on the read if possible)
  6) loop through the columns of the geno table
  7) loop through the rows of the BVs table
'

processTable <- function() {
  write.table(geno, "GenotypesTest.txt", quote=FALSE, row.names=FALSE)
  write.table(BVs, "BreedingValues.txt", quote=FALSE, row.names=FALSE)
  Geno <- read.table("GenotypesTest.txt", header=TRUE, sep= " ")
  BrVs <- read.table("BreedingValues.txt", header=TRUE, sep= " ")
  mean1 <- c()
  genoPart <- geno[, 700:1000]
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
  write_delim(geno, "GenotypesTest.txt", delim=" ", col_names=TRUE)
  write_delim(BVs, "BreedingValues.txt", delim=" ", col_names=TRUE)
  Geno <- read_table2("GenotypesTest.txt", col_names=TRUE)
  BrVs <- read_table2("BreedingValues.txt", col_names=TRUE)
  mean1 <- c()
  genoPart <- geno[, 700:1000]
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
  fwrite(geno, "GenotypesTest.txt", quote=FALSE, row.names=FALSE, sep=" ")
  fwrite(geno, "BreedingValues.txt", quote=FALSE, row.names=FALSE, sep=" ")
  Geno <- fread("GenotypesTest.txt", header=TRUE, sep=" ", select=700:1000)
  BrVs <- fread("BreedingValues.txt", header=TRUE, sep=" ")
  mean1 <- c()
  genoPart <- geno[, 700:1000]
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
                  dims = dim(geno), storage.mode = "integer",
                  chunk = c(3000, 1000), level = 6)
  
  h5save(as.matrix(geno), file = "gPart.h5", name = "genotypes/geno_chunked")
  
  h5createGroup("gPart.h5", "breedingvalues")
  h5createDataset(file = "gPart.h5", dataset = "breedingvalues/bvs_chunked",
                  dims = dim(BVs), storage.mode = "double",
                  level = 6)
  
  h5save(as.matrix(BVs), file = "gPart.h5", name = "breedingvalues/bvs_chunked")
  h5closeAll()
  genoPart <- h5read(file = "gPart.h5", name = "genotypes/geno_chunked", 
                     index = list(NULL, 700:1000))
  BrVs <- h5read(file = "gPart.h5", name = "breedingvalues/bvs_chunked")
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
  write_feather(geno, "GenotypesTest.txt")
  write_feather(BVs, "BreedingValues.txt")
  Geno <- read_feather("GenotypesTest.txt", columns = 700:1000)
  BrVs <- read_feather("BreedingValues.txt")
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

#####################################
# Run and time the functions

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
