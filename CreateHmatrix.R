setwd("~/bin/AlphaMateLinux/OCSSloPop/CowSample/VseKrave/")
library(readr)

#1) a matrika za živali v optimizaciji + vse genotipizirane živali (v napovedi gPV)
a <- readr::read_table2("PedigreeNrm.txt", col_names=FALSE,  col_types = cols(X1 = col_character(), .default = col_double()))
#2) g matrika za vse genotipizirane živali (v napovedi gPV)
g <- readr::read_table2("GenotypeNrm.txt", col_names=FALSE,  col_types = cols(X1 = col_character(), .default = col_double()))

optInd <- read.table("IndOpt.txt")
genoInd <- read.table("IndGeno.txt")

#3) pridobi a12, a21 in a22
#del a matrike za sorodstvo negenotipiziranih z genotipiziranimi živali
a12 <- as.matrix(a[,-1][!(a$X1 %in% genoInd$V1), a$X1 %in% genoInd$V1])
#del a matrike za sorodstvo genotipiziranih z negenotipiziranimi
a21 <- as.matrix(a[,-1][(a$X1 %in% genoInd$V1), !(a$X1 %in% genoInd$V1)])

#del a matrike za genotipizirane živali
a22 <- as.matrix(a[,-1][a$X1 %in% genoInd$V1, a$X1 %in% genoInd$V1])
a22animals <- a$X1[a$X1 %in% genoInd$V1]
#inverza matrike
ia22 <- solve(a22)
#imena stolpcev = zivali
#colnames(ia22) <- a22animals




#4) izračunaj matriko H
#naredic matriko C
#zmnoži
a12_ia22 <- a12 %*% ia22
library(Matrix)
dim(a[-1])
dim(a12_ia22)
diffRow <- dim(a[,-1])[1] - dim(a12_ia22)[1]
diffCol <- dim(a[,-1])[2] - dim(a12_ia22)[2]

upperC <- cbind2(a12_ia22, matrix(0, dim(a12_ia22)[1], diffCol))
lowerC <- cbind2(matrix(0, diffRow, dim(a12_ia22)[1]), diag(diffRow))
C <- rbind2(upperC, lowerC)


#naredi identično matriko 1
I1 <- rbind2(diag(dim(a12_ia22)[1]), diag(diffRow))

#naredic matriko C1
#zmnoži
ia22_a12 <- ia22 %*% a21
dim(a[-1])
dim(ia22_a12)
diffRow <- dim(a[,-1])[1] - dim(ia22_a12)[1]
diffCol <- dim(a[,-1])[2] - dim(ia22_a12)[2]

upperC1 <- cbind2(ia22_a12, matrix(0, dim(ia22_a12)[1], diffCol))
lowerC1 <- cbind2(matrix(0, diffRow, dim(ia22_a12)[1]), diag(diffRow))
C1 <- rbind2(upperC1, lowerC1)

#naredi še identične matrike
I1 <- 
