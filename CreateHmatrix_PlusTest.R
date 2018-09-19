setwd("~/bin/AlphaMateLinux/OCSSloPop/CowSample/VseKrave/")
setwd("~/bin/AlphaRelate/testHMatrix//")
library(readr)

#1) a matrika za živali v optimizaciji + vse genotipizirane živali (v napovedi gPV)
a <- data.matrix(readr::read_table2("PedigreeNrm.txt", col_names=FALSE,  col_types = cols(X1 = col_character(), .default = col_double())))
Aanimals <- a[,1]
a <- a[,-1]
rownames(a) <- Aanimals
colnames(a) <- Aanimals
#2) g matrika za vse genotipizirane živali (v napovedi gPV)
g <- data.matrix(readr::read_table2("GenotypeNrm.txt", col_names=FALSE,  col_types = cols(X1 = col_character(), .default = col_double())))
Ganimals <- g[,1]
g <- g[,-1]
rownames(g) <- Ganimals
colnames(g) <- Ganimals

optInd <- read.table("IndOpt.txt")
genoInd <- read.table("IndGeno.txt")

#order of the animals
sorted <- sort(Aanimals)
NumOrder <- match(Aanimals, sorted)

###
#fiRST REARRANGE SO THAT THE ANIMALS FOLLOW IN A NUMERIC ORDER
a <- a[NumOrder, NumOrder]
allAnimals <- colnames(a)
a11Animals <- allAnimals[!(allAnimals %in% genoInd$V1)]
a22Animals <- allAnimals[allAnimals %in% genoInd$V1]
order <- match(c(a11Animals, a22Animals), allAnimals)


##REARRANGE A!!!!
#PutGenotypedAnimals in the right down corner!!!!
a <- a[order, order]


#3) pridobi a11, a12, a21 in a22
a11 <- as.matrix(a[!(rownames(a) %in% genoInd$V1), !(rownames(a) %in% genoInd$V1)])
#del a matrike za sorodstvo negenotipiziranih z genotipiziranimi živali
a12 <- as.matrix(a[!(rownames(a) %in% genoInd$V1), rownames(a) %in% genoInd$V1])
#del a matrike za sorodstvo genotipiziranih z negenotipiziranimi
a21 <- as.matrix(a[(rownames(a) %in% genoInd$V1), !(rownames(a) %in% genoInd$V1)])

#del a matrike za genotipizirane živali
a22 <- as.matrix(a[rownames(a) %in% genoInd$V1, rownames(a) %in% genoInd$V1])
a22animals <- rownames(a)[rownames(a) %in% genoInd$V1]
#inverza matrike
ia22 <- solve(a22)
#imena stolpcev = zivali
#colnames(ia22) <- a22animals



"""
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
""" 
  
  
##################################33
##################################33
##################################33
gb <- 0.95*g + 0.05*a22
d <- gb - a22
upLeft <- mmult(mmult(mmult(mmult(a12, ia22), d), ia22), a21)
upLeftO <- a12 %*% ia22 %*%  d %*% ia22 %*% a21
downLeft <- d %*% ia22 %*% a21
upRight <- a12 %*% ia22 %*% d
downRight <- d

HaddU <- cbind2(upLeft, upRight)
HaddD <- cbind2(downLeft, downRight)
Hadd <- rbind2(HaddU, HaddD)
dim(Hadd)
H <- a + Hadd
allAnimals <- as.numeric(rownames(a))
reorder <- match(sorted, allAnimals)
allAnimals[reorder]
H <- H[reorder, reorder]

indopt <- read.table("IndOpt.txt")
Hanimals <- colnames(H)
optOrder <- match(indopt$V1, Hanimals)
H <- H[optOrder, optOrder]

write.table(H, "Hmatrix.txt", quote=FALSE, col.names = FALSE, sep=" ")


library(Rcpp)
cppFunction('NumericMatrix mmult(const NumericMatrix& m1, const NumericMatrix& m2){
if (m1.ncol() != m2.nrow()) stop ("Incompatible matrix dimensions");
NumericMatrix out(m1.nrow(),m2.ncol());
NumericVector rm1, cm2;
for (size_t i = 0; i < m1.nrow(); ++i) {
    rm1 = m1(i,_);
    for (size_t j = 0; j < m2.ncol(); ++j) {
      cm2 = m2(_,j);
      out(i,j) = std::inner_product(rm1.begin(), rm1.end(), cm2.begin(), 0.);              
    }
  }
return out;
}')

###check whether blupf90 did the same matrix

Hf90 <- read.table("HinvOrig.txt")
hh <- Hf90[Hf90$V1 %in% allAnimals & Hf90$V2 %in% allAnimals,]
hh <- hh[order(hh$V1, hh$V2),]
hh <- hh[!(is.na(hh$V1)),]

hinv <- solve(H)
write.table(hinv1, "HINV.txt")

AINV <- read.table("AinvOrig.txt")
AINV <- AINV[AINV$V1 %in% allAnimals & AINV$V2 %in% allAnimals,]
AINV <- AINV[order(AINV$V1, AINV$V2),]

AINV
a <- read.table("PedigreeNrm.txt")
aA <- read.table("PedigreeNrm_ALL.txt")
colnames(aA) <- aA$V1
rownames(aA) <- aA$V1
aSmall <- aA[,-1][aA$V1 %in% allAnimals, aA$V1 %in% allAnimals]

ainv <- solve(a[,-1])
colnames(ainv) <- a$V1
rownames(ainv) <- a$V1
aaa <- ainv[rownames(ainv) %in% allAnimals, colnames(ainv) %in% allAnimals]


Gf <- read.table("G_Orig.txt")
g <- read.table("GenotypeNrm.txt")
G <- as.matrix(g)[,-1]
dimnames(G) <- list(as.vector(g[,1]), as.vector(g[,1]))
Gtable <- as.data.frame(as.table(G))

pedf90 <- read.table("Blupf90.ped")
alphaped <- read.table("PedigreeTest.txt", sep=",")

geno <- read.table("TestGeno.txt")
genoF <- as.matrix(geno[,-1])
(sum(genoF==0) / 340) * (sum(genoF==2) / 340) * 2
sum(sapply(geno[,-1], function(x) (sum(x==0) / 4) * (sum(x==2) / 4) * 2))


#centre
genoC <- geno
genoC$V2 <- genoC$V2 - ((colSums(geno[,-1]) / (nrow(geno)*2))*2)[1]
genoC$V3 <- genoC$V3 - ((colSums(geno[,-1]) / (nrow(geno)*2))*2)[2]
genoC <- genoC[,-1]

#scale
c <- sqrt(sum(sapply(geno[,-1], function(x) (  (sum(x) / (length(x)*2)) * (1 - ((sum(x) / (length(x)*2)))) * 2   ))))

genoZ <- genoC / c
genoZt <- t(genoZ)

Gm <- as.matrix(genoZ) %*% as.matrix(genoZt)

gb <- 0.95 * Gm + 0.05*a22
#scale- locus specific
c <- sqrt((sapply(geno[,-1], function(x) (  (sum(x) / (length(x)*2)) * (1 - ((sum(x) / (length(x)*2)))) * 2   ))))
cm <- rbind(c, c)
genoZ <- genoC / cm

genoZt <- t(genoZ)

Gm <- as.matrix(genoZ) %*% as.matrix(genoZt) / 2



Ginv <- read.table("G")

Gf[order(Gf$V2, Gf$V1),]
y <- matrix(0, nrow=4, ncol=4) 
y[1:4, 1:4] <- Gf$V3
Gf[Gf$V1 == Gf$V2,]
gfi <- solve(y)
ia22

solve(0.95 * gfi + 0.05*ia22)

solve(0.95 * solve(gb) + 0.05*ia22)
