library(pedigree)
library(pedigreemm)
 varE <- 40
 varA <- 20
 alpha <- varE / varA

blupAcc <- function(X, noAnObs, y, ped, alpha) {
  dat <- data.frame(id=1:nrow(ped), data=c(NA,NA,NA,y))
  dat$id = factor(dat$id, levels=ped$id)
  Z = model.matrix(dat$data~dat$id)
  XtZ = t(X)%*%Z
  ZtX = t(XtZ)
  Xty <- t(X) %*% y
  Zty <- t(Z) %*% y
  ZtZ <- t(Z) %*% Z
  XtX <- t(X) %*% X
  
  getAInv(ped)
  
  AinvAlpha <- Ainv*alpha
  gen <- ZtZ + AinvAlpha
  
  LSEU <- cbind(XtX, XtZ)
  LSEL <- cbind(ZtX, gen)
  LSE <- rbind(LSEU, LSEL)
  yT <- rbind(Xty, Zty)
  
  
  sol <- solve(LSE) %*% yT
  coeff <- solve(LSE)
  
  diagonal <- diag(coeff)[(ncol(X)+1):(ncol(X) + ncol(Z))]
  accuracies <- round(c(sqrt(1 - diagonal*alpha)),3)
  SEP <- round(c(sqrt(diagonal*varE)),3)
  return(cbind(accuracies, SEP))
}

ped <- data.frame(id=1:8, dam=c(NA,NA,NA,NA,2,2,5,6), sire=c(NA,NA,NA,1,3,1,4,3))
X <- matrix (nrow = 5 , ncol = 2, c(c(1,0,0,1,1), c(0,1,1,0,0)))
noAnObs <- 5
y <- as.vector(c(4.5, 2.9, 3.9, 3.5, 5.0))
blup <- function(X, Z, y, ped, alpha) {
  XtZ = t(X)%*%Z
  ZtX = t(XtZ)
  Xty <- t(X) %*% y
  Zty <- t(Z) %*% y
  ZtZ <- t(Z) %*% Z
  XtX <- t(X) %*% X

  makeAinv(ped)
  Ai <- read.table("Ainv.txt")
  nInd <- nrow(ped)
  Ainv <- matrix(0,nrow = nInd,ncol = nInd)
  Ainv[as.matrix(Ai[,1:2])] <- Ai[,3]
  dd <- diag(Ainv)
  Ainv <- Ainv + t(Ainv)
  diag(Ainv) <- dd
  
  AinvAlpha <- Ainv*2
  gen <- ZtZ + AinvAlpha
  
  LSEU <- cbind(XtX, XtZ)
  LSEL <- cbind(ZtX, gen)
  LSE <- rbind(LSEU, LSEL)
  yT <- rbind(Xty, Zty)

  
  sol <- solve(LSE) %*% yT
  coeff <- solve(LSE)
  
  diagonal <- diag(coeff)[3:length(diag(coeff))]
  accuracies <- round(c(sqrt(1 - diagonal*alpha)),3)
  return(accuracies)
}

alpha <- 2

AccDF <- data.frame(id = ped$id, basic = blupAcc(X, 5, y, ped, alpha))



#####################################################################################
#Mrode example - prepare blupf90 files
#ped code: 1 - both parents known,2 - one parent known, 3 - both parents unknown --> HOWEVER, unknown parents can be assigned parent group numbers
#ped <- data.frame(id=1:9, sire=c(0,0,0,1,3,1,4,3,3), dam=c(0,0,0,0,2,2,5,6,2), code=c(3,3,3,2,1,1,1,1,1))
#dat <- data.frame(animal = 1:8, p = c(0,0,0,4.5, 2.9, 3.9, 3.5, 5.0), sex = c(1,2,1,1,2,2,1,1))
#write.table(dat, '/home/jana/Documents/PhD/Mrode.dat', col.names=F, sep=" ", row.names=F, quote=F)
#write.table(ped, '/home/jana/Documents/PhD/Mrode.ped', col.names=F, sep=" ", row.names=F, quote=F)


#add full-sib to animal 5
#ped[9,] <- c(9, 3, 2, 1)
#female full-sib
#dat[9,] <- c(9, 3.5, 2)


#male full-sib
#dat[9,] <- c(9, 3.5, 2)
#what happens to accuracies


#add half-sib to animal 5 - FEMALE
X <- matrix (nrow = 6 , ncol = 2, c(c(1,0,0,1,1,0), c(0,1,1,0,0,1)))

y <- as.vector(c(4.5, 2.9, 3.9, 3.5, 5.0, 3.5))
ped <- data.frame(id=1:9, dam=c(NA,NA,NA,NA,2,2,5,6,5), sire=c(NA,NA,NA,1,3,1,4,3,1))



#add male full-sib

ped[9,] <- c(9, 3, 2)
X <- matrix (nrow = 6 , ncol = 2, c(c(1,0,0,1,1,1),c(0,1,1,0,0,0)))
Z <- matrix (nrow = 6 , ncol = 9, c(rep(0,18), 1, rep(c(rep(0, 6),1),5)))
y <- as.vector(c(4.5, 2.9, 3.9, 3.5, 5.0, 3.5))

blup(X, Z, y, ped, alpha)

AccDFTemp <- data.frame(id = ped$id, FeMaleHalfSib = blupAcc(X, 6, y, ped, alpha))
AccDF <- merge(AccDF, AccDFTemp, by='id', all=T)

#add half-sib to animal 5 - MALE
X <- matrix (nrow = 6 , ncol = 2, c(c(1,0,0,1,1,1), c(0,1,1,0,0,0)))
y <- as.vector(c(4.5, 2.9, 3.9, 3.5, 5.0, 3.5))
ped <- data.frame(id=1:9, dam=c(NA,NA,NA,NA,2,2,5,6,5), sire=c(NA,NA,NA,1,3,1,4,3,1))

AccDFTemp <- data.frame(id = ped$id, MaleHalfSib = blupAcc(X, 6, y, ped, alpha))
AccDF <- merge(AccDF, AccDFTemp, by='id', all=T)

#add full-sib to animal 5 - FEMALE
X <- matrix (nrow = 6 , ncol = 2, c(c(1,0,0,1,1,0), c(0,1,1,0,0,1)))
y <- as.vector(c(4.5, 2.9, 3.9, 3.5, 5.0, 3.5))
ped <- data.frame(id=1:9, dam=c(NA,NA,NA,NA,2,2,5,6,2), sire=c(NA,NA,NA1,3,1,4,3,3))

AccDFTemp <- data.frame(id = ped$id, FemaleFulSib = blupAcc(X, 6, y, ped, alpha))
AccDF <- merge(AccDF, AccDFTemp, by='id', all=T)

#add full-sib to animal 5 - MALE
X <- matrix (nrow = 6 , ncol = 2, c(c(1,0,0,1,1,1), c(0,1,1,0,0,0)))
y <- as.vector(c(4.5, 2.9, 3.9, 3.5, 5.0, 3.5))
ped <- data.frame(id=1:9, dam=c(NA,NA,NA,NA,2,2,5,6,2), sire=c(NA,NA,NA,1,3,1,4,3,3))

AccDFTemp <- data.frame(id = ped$id, MaleFulSib = blupAcc(X, 6, y, ped, alpha))
AccDF <- merge(AccDF, AccDFTemp, by='id', all=T)



#add another offspring - female
X <- matrix (nrow = 7 , ncol = 2, c(c(1,0,0,1,1,1, 0), c(0,1,1,0,0,0,1)))
y <- as.vector(c(4.5, 2.9, 3.9, 3.5, 5.0, 3.5, 4.7))
ped <- data.frame(id=1:10, dam=c(NA,NA,NA,NA,2,2,5,6,2,6), sire=c(NA,NA,NA,1,3,1,4,3,3,7))

AccDFTemp <- data.frame(id = ped$id, Offspring1 = blupAcc(X, 7, y, ped, alpha))
AccDF <- merge(AccDF, AccDFTemp, by='id', all=T)

#add another offspring - male
X <- matrix (nrow = 8 , ncol = 2, c(c(1,0,0,1,1,1, 0, 1), c(0,1,1,0,0,0,1,0)))
y <- as.vector(c(4.5, 2.9, 3.9, 3.5, 5.0, 3.5, 4.7, 3.1))
ped <- data.frame(id=1:11, dam=c(NA,NA,NA,NA,2,2,5,6,2,6,5), sire=c(NA,NA,NA,1,3,1,4,3,3,7,7))

AccDFTemp <- data.frame(id = ped$id, Offspring2 = blupAcc(X, 8, y, ped, alpha))
AccDF <- merge(AccDF, AccDFTemp, by='id', all=T)

#basic scenario 2
X <- matrix (nrow = 8 , ncol = 2, c(c(1,0,0,1,1,1, 0, 1), c(0,1,1,0,0,0,1,0)))
y <- as.vector(c(4.5, 2.9, 3.9, 3.5, 5.0, 3.5, 4.7, 3.1))
ped <- data.frame(id=1:11, dam=c(NA,NA,NA,NA,2,2,5,6,5,6,5), sire=c(NA,NA,NA,1,3,1,4,3,1,7,7))

AccDFTemp <- data.frame(id = ped$id, Basic2 = blupAcc(X, 8, y, ped, alpha))
AccDF <- merge(AccDF, AccDFTemp, by='id', all=T)

#add another offspring - male
X <- matrix (nrow = 9 , ncol = 2, c(c(1,0,0,1,1,1, 0, 1,1), c(0,1,1,0,0,0,1,0,0)))
y <- as.vector(c(4.5, 2.9, 3.9, 3.5, 5.0, 3.5, 4.7, 3.1,3.2))
ped <- data.frame(id=1:12, dam=c(NA,NA,NA,NA,2,2,5,6,5,6,5,5), sire=c(NA,NA,NA,1,3,1,4,3,1,7,7,3))

AccDFTemp <- data.frame(id = ped$id, Offspring3 = blupAcc(X, 9, y, ped, alpha))
AccDF <- merge(AccDF, AccDFTemp, by='id', all=T)

#add another offspring - female
X <- matrix (nrow = 10 , ncol = 2, c(c(1,0,0,1,1,1, 0, 1,1,0), c(0,1,1,0,0,0,1,0,0,1)))
y <- as.vector(c(4.5, 2.9, 3.9, 3.5, 5.0, 3.5, 4.7, 3.1,3.2, 4.4))
ped <- data.frame(id=1:13, dam=c(NA,NA,NA,NA,2,2,5,6,5,6,5,5,10), sire=c(NA,NA,NA,1,3,1,4,3,1,7,7,3,3))

AccDFTemp <- data.frame(id = ped$id, Offspring4 = blupAcc(X, 10, y, ped, alpha))
AccDF <- merge(AccDF, AccDFTemp, by='id', all=T)

#add another offspring - female
X <- matrix (nrow = 11 , ncol = 2, c(c(1,0,0,1,1,1, 0, 1,1,0,0), c(0,1,1,0,0,0,1,0,0,1,1)))
y <- as.vector(c(4.5, 2.9, 3.9, 3.5, 5.0, 3.5, 4.7, 3.1,3.2, 4.4, 3.0))
ped <- data.frame(id=1:14, dam=c(NA,NA,NA,NA,2,2,5,6,5,6,5,5,10, 10), sire=c(NA,NA,NA,1,3,1,4,3,1,7,7,3,3, 12))

AccDFTemp <- data.frame(id = ped$id, Offspring5 = blupAcc(X, 11, y, ped, alpha))
AccDF <- merge(AccDF, AccDFTemp, by='id', all=T)

#add unrelated male
X <- matrix (nrow = 12 , ncol = 2, c(c(1,0,0,1,1,1, 0, 1,1,0,0,1), c(0,1,1,0,0,0,1,0,0,1,1,0)))
y <- as.vector(c(4.5, 2.9, 3.9, 3.5, 5.0, 3.5, 4.7, 3.1,3.2, 4.4, 3.0, 4.3))
ped <- data.frame(id=1:15, dam=c(NA,NA,NA,NA,2,2,5,6,5,6,5,5,10, 10, NA), sire=c(NA,NA,NA,1,3,1,4,3,1,7,7,3,3, 12, NA))

AccDFTemp <- data.frame(id = ped$id, UnrelatedM = blupAcc(X, 12, y, ped, alpha))
AccDF <- merge(AccDF, AccDFTemp, by='id', all=T)

#add offspirng of unrelated male
X <- matrix (nrow = 13 , ncol = 2, c(c(1,0,0,1,1,1, 0, 1,1,0,0,1,1), c(0,1,1,0,0,0,1,0,0,1,1,0,0)))
y <- as.vector(c(4.5, 2.9, 3.9, 3.5, 5.0, 3.5, 4.7, 3.1,3.2, 4.4, 3.0, 4.3, 2.9))
ped <- data.frame(id=1:16, dam=c(NA,NA,NA,NA,2,2,5,6,5,6,5,5,10, 10, NA, 10), sire=c(NA, NA, NA,1,3,1,4,3,1,7,7,3,3, 12, NA,15))

AccDFTemp <- data.frame(id = ped$id, OffspringUnrel = blupAcc(X, 13, y, ped, alpha))
AccDF <- merge(AccDF, AccDFTemp, by='id', all=T)

#add offspirng of unrelated male - different EBV
X <- matrix (nrow = 13 , ncol = 2, c(c(1,0,0,1,1,1, 0, 1,1,0,0,1,1), c(0,1,1,0,0,0,1,0,0,1,1,0,0)))
y <- as.vector(c(4.5, 2.9, 3.9, 3.5, 5.0, 3.5, 4.7, 3.1,3.2, 4.4, 3.0, 4.3, 4.1))
ped <- data.frame(id=1:16, dam=c(NA,NA,NA,NA,2,2,5,6,5,6,5,5,10, 10, NA, 10), sire=c(NA,NA,NA,1,3,1,4,3,1,7,7,3,3, 12, NA,15))

AccDFTemp <- data.frame(id = ped$id, OffspringUnrelA = blupAcc(X, 13, y, ped, alpha))
AccDF <- merge(AccDF, AccDFTemp, by='id', all=T)

write.csv(AccDF, '/home/jana/Documents/PhD/Accuracy_EBV.csv', quote=F, row.names=F)

