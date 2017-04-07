rcorr <- function(x, r){ 
  if (FALSE) {    
    x <- rnorm(100)   
    y <- rcorr(x=x, r=.5)  
    cor(x,y) 
    y <- rcorr(x=x, r=1)  
    cor(x,y)  } 
  r*x + rnorm(length(x), mean=0, sd=sqrt(1 - r**2)) 
  } 
x <- rnorm(1000) 
y <- rcorr(x, r=0.6) 
cor(x, y)
plot(x, y)



exact_corr = function (x,rho) {   
  theta = acos(rho) 
  x2 = rnorm(length(x),2,0.5) 
  X = cbind(x,x2) 
  Xctr  <- scale(X, center=TRUE, scale=FALSE) 
  Id   <- diag(length(x))  
  Q    <- qr.Q(qr(Xctr[ , 1, drop=FALSE]))    
  P    <- tcrossprod(Q)        
  # = Q Q' 
  x2o  <- (Id-P) %*% Xctr[ , 2]       
  Xc2  <- cbind(Xctr[ , 1], x2o)      
  Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  
  Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  
  x <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]      
  return(x) } 

x=rnorm(500) 
y = exact_corr(x,0.5) 
cor(x,y)