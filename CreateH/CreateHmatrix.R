library(readr)
library(methods)

f <- function () {
  
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
  rownames(g) <-Ganimals
  colnames(g) <-Ganimals
  
  #genoInd <- read.table("IndForGeno.txt")
  
  
  
  ###
  #fiRST REARRANGE SO THAT THE ANIMALS FOLLOW IN A NUMERIC ORDER
  #order of the animals
  sorted <- sort(as.numeric(Aanimals))
  
  ###
  #fiRST REARRANGE SO THAT THE ANIMALS FOLLOW IN A NUMERIC ORDER
  a11Animals <- Aanimals[!(Aanimals %in% Ganimals)]
  order <-match(c(sort(as.numeric(a11Animals)), Ganimals), Aanimals)
  
  
  ##REARRANGE A!!!!
  #PutGenotypedAnimals in the right down corner!!!!
  a <-a[order, order]
  
  
  
  #3) pridobi a11, a12, a21 in a22
  #a11 <- as.matrix(a[!(rownames(a) %in% Ganimals), !(rownames(a) %in% Ganimals)])
  #del a matrike za sorodstvo negenotipiziranih z genotipiziranimi živali
  a12 <- as.matrix(a[!(rownames(a) %in% Ganimals), rownames(a) %in% Ganimals])
  #del a matrike za sorodstvo genotipiziranih z negenotipiziranimi
  #a21 <- as.matrix(a[(rownames(a) %in% Ganimals), !(rownames(a) %in% Ganimals)])
  
  #del a matrike za genotipizirane živali
  a22 <- as.matrix(a[rownames(a) %in% Ganimals, rownames(a) %in% Ganimals])
  #a22animals <- rownames(a)[rownames(a) %in% Ganimals]
  #inverza matrike
  ia22 <-solve(a22)
  #imena stolpcev = zivali
  
  
  
  ##################################33
  gb <- 0.95*g + 0.05*a22
  
  #obtain alpha and beta - scale G to A
  avgA22 <- mean(a22)
  avgG <- mean(gb)
  beta <- (mean(diag(a22)) - avgA22) / (mean(diag(gb)) - avgG)
  alpha <- avgA22 - (beta * avgG)
  
  Ga <- alpha + beta*gb
  all.equal(mean(Ga), mean(a22))
  
  d <- Ga - a22
  downLeft <- d %*% ia22 %*% as.matrix(a[(rownames(a) %in% Ganimals), !(rownames(a) %in% Ganimals)])
  upLeft <- a12 %*% ia22 %*% downLeft
  upRight <- t(downLeft)
  #downRight <- d
  
  
  
  HaddU <-cbind2(upLeft, upRight)
  HaddD <-cbind2(downLeft, d)
  Hadd <-rbind2(HaddU, HaddD)
  
  
  
  
  H <- a + Hadd
  allAnimals <- as.numeric(rownames(a))
  reorder <- match(sorted, as.numeric(rownames(H)))
  H  <- H[reorder, reorder]
  
  indopt <-read.table("IndOpt.txt")
  optOrder <-match(indopt$V1, colnames(H))
  H <-H[optOrder, optOrder]
  
  write.table(H, "Hmatrix.txt", quote=FALSE, col.names = FALSE, sep=" ")
}


Rprof("Rprof.out", interval = 0.01, line.profiling = TRUE)
f()
Rprof(NULL)
print(summaryRprof("Rprof.out"), lines="show")

