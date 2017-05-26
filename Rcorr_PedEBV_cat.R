# rcorr <- function(x, r){ 
#   if (FALSE) {    
#     x <- rnorm(100)   
#     y <- rcorr(x=x, r=.5)  
#     cor(x,y) 
#     y <- rcorr(x=x, r=1)  
#     cor(x,y)  } 
#   r*x + rnorm(length(x), mean=0, sd=sqrt(1 - r**2)) 
# } 
# 
# a <- rcorr(1:5, 0.5)
# cor(a, 1:5)

rcorr = function (x,rho) {    
  theta = acos(rho)  
  x2 = rnorm(length(x),2,0.5)  
  X = cbind(x,x2)  
  Xctr  <- scale(X, center=TRUE, scale=FALSE)  
  Id   = diag(length(x))      
  Q    = qr.Q(qr(Xctr[ , 1, drop=FALSE]))        
  P    = tcrossprod(Q)          
  # = Q Q'        
  x2o  = (Id-P) %*% Xctr[ , 2]                  
  Xc2  = cbind(Xctr[ , 1], x2o)               
  Y    = Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  
  Y    = Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  
  x = Y[ , 2] + (1 / tan(theta)) * Y[ , 1]      
  return(x) 
}






#read in AlphaSim PedigreeAndGenetic values, ki si jim dodala kategorije
ped <- read.table("AlphaSimPed", header=T)
ped <- read.csv("/home/jana/bin/AlphaSim1.05Linux/SimulatedData/PedigreeAndGeneticValues_cat.txt", header=T)
#stolpci Gen, Ind, father, Mother, gvNormUnrest, cat
ped <- ped[,c(1,2,3,4,9,13)]

#tukaj nastavi tocnosti za posamezne kategorije
catAcc <- data.frame(cat = unique(ped$cat), acc = NA)
catAcc$acc[catAcc$cat=='nr']  <- 0.35
catAcc$acc[catAcc$cat=='telF']  <- 0.35
catAcc$acc[catAcc$cat=='pt']  <- 0.35
catAcc$acc[catAcc$cat=='k']  <- 0.35
catAcc$acc[catAcc$cat=='pBM']  <- 0.35
catAcc$acc[catAcc$cat=='bm']  <- 0.35
catAcc$acc[catAcc$cat=='potomciNP']  <- 0.35
catAcc$acc[catAcc$cat=='vhlevljeni'] <- 0.35
catAcc$acc[catAcc$cat=='mladi']  <- 0.35
catAcc$acc[catAcc$cat=='cak']  <- 0.35
catAcc$acc[catAcc$cat=='pb']  <- 0.35
catAcc$acc[catAcc$cat=='telM']  <- 0.35
catAcc$acc[catAcc$cat=='bik12']  <- 0.35
catAcc$acc[catAcc$cat=='pripust1']  <- 0.35
catAcc$acc[catAcc$cat=='pripust2']  <- 0.35
catAcc$acc[catAcc$cat=='izl']  <- 0

#nastavi accuracies po kategorijah
for (cat in catAcc$cat) {
  ped$EBV[ped$cat == cat] <- rcorr(ped$gvNormUnres1[ped$cat == cat], r=catAcc$acc[catAcc$cat==cat]) 
  print(paste(cat, cor(ped$gvNormUnres1[ped$cat == cat], ped$EBV[ped$cat == cat])), sep=" ")
}

#pridobi EBV iz TBV
cor(ped$gvNormUnres1, ped$EBV)
write.csv(ped[,c(1,2,3,4,5, 7)], "GenPed_EBV.txt", quote=F, row.names=F)

