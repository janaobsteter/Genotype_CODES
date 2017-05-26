library(pedigree)
X <- matrix (nrow = 5 , ncol = 2, c(1,0,0,1,1,0,1,1,0,0))
Z <- matrix (nrow = 5 , ncol = 8, c(rep(0,15), 1, rep(c(rep(0, 5),1),4)))
y <- as.vector(c(4.5, 2.9, 3.9, 3.5, 5.0))

XtZ = t(X)%*%Z
ZtX = t(XtZ)
Xty <- t(X) %*% y
Zty <- t(Z) %*% y
ZtZ <- t(Z) %*% Z
XtX <- t(X) %*% X

ped <- data.frame(id=1:8, dam=c(0,0,0,0,2,2,5,6), sire=c(0,0,0,1,3,1,4,3))
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

alpha <- 2
sol <- solve(LSE) %*% yT
coeff <- solve(LSE)

diagonal <- diag(coeff)[3:10]
accuracies <- round(c(sqrt(1 - diagonal*alpha)),3)

#####################################################################################
#Mrode example - prepare blupf90 files
#ped code: 1 - both parents known,2 - one parent known, 3 - both parents unknown --> HOWEVER, unknown parents can be assigned parent group numbers
ped <- data.frame(id=1:9, sire=c(0,0,0,1,3,1,4,3,3), dam=c(0,0,0,0,2,2,5,6,2), code=c(3,3,3,2,1,1,1,1,1))
dat <- data.frame(animal = 1:8, p = c(0,0,0,4.5, 2.9, 3.9, 3.5, 5.0), sex = c(1,2,1,1,2,2,1,1))
write.table(dat, '/home/jana/Documents/PhD/Mrode.dat', col.names=F, sep=" ", row.names=F, quote=F)
write.table(ped, '/home/jana/Documents/PhD/Mrode.ped', col.names=F, sep=" ", row.names=F, quote=F)

#add full-sib to animal 5
ped[9,] <- c(9, 3, 2, 1)
#female full-sib
dat[9,] <- c(9, 3.5, 2)
#male full-sib
dat[9,] <- c(9, 3.5, 2)
#what happens to accuracies


X <- matrix (nrow = 6 , ncol = 2, c(1,0,0,1,1,0,0,1,1,0,0,1))
Z <- matrix (nrow = 6 , ncol = 9, c(rep(0,18), 1, rep(c(rep(0, 6),1),5)))
y <- as.vector(c(4.5, 2.9, 3.9, 3.5, 5.0, 3.5))

XtZ = t(X)%*%Z
ZtX = t(XtZ)
Xty <- t(X) %*% y
Zty <- t(Z) %*% y
ZtZ <- t(Z) %*% Z
XtX <- t(X) %*% X

ped <- data.frame(id=1:9, dam=c(0,0,0,0,2,2,5,6,2), sire=c(0,0,0,1,3,1,4,3,3))
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

alpha <- 2
sol <- solve(LSE) %*% yT
coeff <- solve(LSE)

diagonal <- diag(coeff)[3:11]
accuracies <- round(c(sqrt(1 - diagonal*alpha)),3)
