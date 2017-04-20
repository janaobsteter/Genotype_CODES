rcorr <- function(x, r){ 
  if (FALSE) {    
    x <- rnorm(100)   
    y <- rcorr(x=x, r=.5)  
    cor(x,y) 
    y <- rcorr(x=x, r=1)  
    cor(x,y)  } 
  r*x + rnorm(length(x), mean=0, sd=sqrt(1 - r**2)) 
  } 




ped <- read.table("AlphaSimPed", header=T)
ped <- ped[,c(1,2,3,4,9)]
ped$EBV <- rcorr(ped$gvNormUnres1, r=0.9) #pridobi EBV iz TBV
cor(ped$gvNormUnres1, ped$EBV)
write.csv(ped, "GenPed_EBV.txt", quote=F, row.names=F)

