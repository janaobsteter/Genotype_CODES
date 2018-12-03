library(readr)
library(methods)




#1) A matrix for the animals (in optimisation / estimation of BVs)
#the PedigreeNrm matrix was create with AlphaRelate 
a <- data.matrix(readr::read_table2("PedigreeNrm.txt", col_names=FALSE,  col_types = cols(X1 = col_character(), .default = col_double())))
#the animals' IDs
Aanimals <- a[,1]
#the animals' genotypes
a <- a[,-1]
#name rows and columns of A-matrix by animals IDs
rownames(a) <- Aanimals
colnames(a) <- Aanimals

#2) G matrix for all the genotyped animals (/ in estimation of BVs) 
#GenotypeNrm was created with AlphaRelate software
g <- data.matrix(readr::read_table2("GenotypeNrm.txt", col_names=FALSE,  col_types = cols(X1 = col_character(), .default = col_double())))
#get animals' IDs
Ganimals <- g[,1]
#get animals' genotyped
g <- g[,-1]
#name rows and columns of G-matrix by animals IDs
rownames(g) <- Ganimals
colnames(g) <- Ganimals


##########################################
#this part served only for rearranging the rows and columns of the matrices to order the animal IDs (numerical) and to put genotyped animals at the end (edges of the matrix)
############################################
#order of the animals
sorted <- sort(as.numeric(Aanimals))

#fiRST REARRANGE SO THAT THE ANIMALS FOLLOW IN A NUMERIC ORDER
a11Animals <- Aanimals[!(Aanimals %in% Ganimals)]
order <-match(c(sort(as.numeric(a11Animals)), Ganimals), Aanimals)


##REARRANGE A!!!!
#PutGenotypedAnimals in the right down corner!!!!
a <-a[order, order]

##################################################

#3) create a11, a12, a21 in a22 - BUT you only need a12 for the subsequent computation
#a11 <- as.matrix(a[!(rownames(a) %in% Ganimals), !(rownames(a) %in% Ganimals)])
#part of the Amatrix for relationship of un-genotyped with genotyped animals
a12 <- as.matrix(a[!(rownames(a) %in% Ganimals), rownames(a) %in% Ganimals])
#part of the Amatrix for relationship of genotyped with un-genotyped animals
#a21 <- as.matrix(a[(rownames(a) %in% Ganimals), !(rownames(a) %in% Ganimals)])

#part of the Amatrix for genotyped animals
a22 <- as.matrix(a[rownames(a) %in% Ganimals, rownames(a) %in% Ganimals])
#matrix inverse
ia22 <-solve(a22)

#column names = animals

##################################
##################################
##################################
#blending of the G-matrix
gb <- 0.95*g + 0.05*a22

#obtain alpha and beta - scale G to A
avgA22 <- mean(a22)
avgG <- mean(gb)
beta <- (mean(diag(a22)) - avgA22) / (mean(diag(gb)) - avgG)
alpha <- avgA22 - (beta * avgG)

Ga <- alpha + beta*gb
all.equal(mean(Ga), mean(a22))

#subtract due to "double counting"
d <- Ga - a22

#create parts of the matrix
downLeft <- d %*% ia22 %*% as.matrix(a[(rownames(a) %in% Ganimals), !(rownames(a) %in% Ganimals)])
upLeft <- a12 %*% ia22 %*% downLeft
upRight <- t(downLeft)
#downRight <- d

#combine parts of the matrix that is added to A matrix (Hadd; H = A + Hadd)
HaddU <-cbind2(upLeft, upRight)
HaddD <-cbind2(downLeft, d)
Hadd <-rbind2(HaddU, HaddD)

H <- a + Hadd
allAnimals <- as.numeric(rownames(a))
reorder <- match(sorted, as.numeric(rownames(H)))
H  <- H[reorder, reorder]



#this part is for when creating Hmatrix for OCS
'''
#extract only the part of Hmatrix corresponding to the animals in optimisation
indopt <- read.table("IndOpt.txt")
Hanimals <- colnames(H)
optOrder <- match(indopt$V1, Hanimals)
H <- H[optOrder, optOrder]
'''


write.table(H, "Hmatrix.txt", quote=FALSE, col.names = FALSE, sep=" ")
