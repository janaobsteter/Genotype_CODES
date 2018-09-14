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

genoInd <- read.table("IndGeno.txt")

#order of the animals
sorted <- sort(Aanimals)
NumOrder <- match(Aanimals, sorted)

###
#fiRST REARRANGE SO THAT THE ANIMALS FOLLOW IN A NUMERIC ORDER
a <- a[NumOrder, NumOrder]
allAnimals <- colnames(a)
a11Animals <- allAnimals[!(allAnimals %in% Ganimals)]
a22Animals <- allAnimals[allAnimals %in% Ganimals]
order <- match(c(a11Animals, a22Animals), allAnimals)


##REARRANGE A!!!!
#PutGenotypedAnimals in the right down corner!!!!
a <- a[order, order]


#3) pridobi a11, a12, a21 in a22
a11 <- as.matrix(a[!(rownames(a) %in% Ganimals), !(rownames(a) %in% Ganimals)])
#del a matrike za sorodstvo negenotipiziranih z genotipiziranimi živali
a12 <- as.matrix(a[!(rownames(a) %in% Ganimals), rownames(a) %in% Ganimals])
#del a matrike za sorodstvo genotipiziranih z negenotipiziranimi
a21 <- as.matrix(a[(rownames(a) %in% Ganimals), !(rownames(a) %in% Ganimals)])

#del a matrike za genotipizirane živali
a22 <- as.matrix(a[rownames(a) %in% Ganimals, rownames(a) %in% Ganimals])
a22animals <- rownames(a)[rownames(a) %in% Ganimals]
#inverza matrike
ia22 <- solve(a22)
#imena stolpcev = zivali

##################################33
##################################33
##################################33
gb <- 0.95*g + 0.05*a22
d <- gb - a22
upLeft <- a12 %*% ia22 %*%  d %*% ia22 %*% a21
downLeft <- d %*% ia22 %*% a21
upRight <- a12 %*% ia22 %*% d
downRight <- d

HaddU <- cbind2(upLeft, upRight)
HaddD <- cbind2(downLeft, downRight)
Hadd <- rbind2(HaddU, HaddD)
H <- a + Hadd
allAnimals <- as.numeric(rownames(a))
reorder <- match(sorted, allAnimals)
H <- H[reorder, reorder]


write.table(H, "Hmatrix.txt", quote=FALSE, col.names = FALSE, sep=" ")