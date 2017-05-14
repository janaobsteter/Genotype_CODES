rcorr <- function(x, r){ 
  if (FALSE) {    
    x <- rnorm(100)   
    y <- rcorr(x=x, r=.5)  
    cor(x,y) 
    y <- rcorr(x=x, r=1)  
    cor(x,y)  } 
  r*x + rnorm(length(x), mean=0, sd=sqrt(1 - r**2)) 
  } 



getBiCop <- function(n, rho, mar.fun=rnorm, x = NULL, ...) {
  if (!is.null(x)) {X1 <- x} else {X1 <- mar.fun(n, ...)}
  if (!is.null(x) & length(x) != n) warning("Variable x does not have the same length as n!")
  
  C <- matrix(rho, nrow = 2, ncol = 2)
  diag(C) <- 1
  
  C <- chol(C)
  
  X2 <- mar.fun(n)
  X <- cbind(X1,X2)
  
  # induce correlation (does not change X1)
  df <- X %*% C
  
  ## if desired: check results
  #all.equal(X1,X[,1])
  #cor(X)
  
  return(df)
}

correlatedValue = function(x, r){
  r2 = r**2
  ve = 1-r2
  SD = sqrt(ve)
  e  = rnorm(length(x), mean=0, sd=SD)
  y  = r*x + e
  return(y)
}

set.seed(5)

y = correlatedValue(x=ped$gvNormUnres1, r=.5)
cor(ped$gvNormUnres1, y)


ped <- read.table("/home/jana/bin/AlphaSim1.05Linux/SimulatedData/PedigreeAndGeneticValues.txt", header=T)
ped <- read.table("/home/jana/bin/AlphaSim1.05Linux/PedigreeAndGeneticValues_BURNIN.txt", header=T)
ped <- ped[,c(1,2,3,4,9)]
ped$EBV <- getBiCop(n = 67000, x=ped$gvNormUnres1, rho=0.9) #pridobi EBV iz TBV

ped$EBV <- rcorr(ped$gvNormUnres1, r=0.9) #pridobi EBV iz TBV
cor(ped$gvNormUnres1, ped$EBV)
write.csv(ped, "GenPed_EBV.txt", quote=F, row.names=F)

