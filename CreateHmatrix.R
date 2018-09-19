library(readr)
library(Rcpp)
library(FastGP)
library(methods)

setwd("~/bin/AlphaRelate/testHMatrix//")
#funkcija za multiplikacijo matrik
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

#genotipizirane živali - v bistvu Ganimals
#genoInd <- read.table("IndGeno.txt")


##########################################
#ta del je samo preurejanje stolpcev in vrstic, da dobim posameznike v numeričnem vrstnem redu + dam genotipizirane živali na konec
############################################
#order of the animals - 
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
##################################################

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
ia22 <- rcppeigen_invert_matrix(a22)
#imena stolpcev = zivali

##################################
##################################
##################################
#blend matrik
gb <- 0.95*g + 0.05*a22
#odštej zaradi "double counting"
d <- gb - a22
#ustvari dele matrike
upLeft <- mmult(mmult(mmult(mmult(a12, ia22), d), ia22), a21)
downLeft <- mmult(mmult(d, ia22), a21)
upRight <- mmult(mmult(a12, ia22), d)
downRight <- d

#poveži dele matrike, ki jo dodaš A matriki (Hadd \ H = A + Hadd)
HaddU <- cbind2(upLeft, upRight)
HaddD <- cbind2(downLeft, downRight)
Hadd <- rbind2(HaddU, HaddD)
H <- a + Hadd
allAnimals <- as.numeric(rownames(a))
reorder <- match(sorted, allAnimals)
H <- H[reorder, reorder]

#tukaj izvleči le del H matrike, ki pripada živalim, ki gredo v optimizacijo v AlphaMate
indopt <- read.table("IndOpt.txt")
Hanimals <- colnames(H)
optOrder <- match(indopt$V1, Hanimals)
H <- H[optOrder, optOrder]

write.table(H, "Hmatrix.txt", quote=FALSE, col.names = FALSE, sep=" ")