### descStat.R
###------------------------------------------------------------------------
### $Id$
###------------------------------------------------------------------------

## --- Unit tests ---

test(descStat) <- function() {
  
  ## Usage for one variable  
  v <- 1:10
  v[c(1, 2, 3)] <- NA
  z <- descStat(x=v)
  
  ## Check if we have the same unit tests
  (checkEquals(z["n"],             10, checkNames=FALSE))
  (checkEquals(z["obs"],            7, checkNames=FALSE))
  (checkEquals(z["mean"],           7, checkNames=FALSE))
  (checkEquals(z["median"],         7, checkNames=FALSE))
  (checkEquals(z["obs"],            7, checkNames=FALSE))
  (checkEquals(round(z["sd"], 1), 2.2, checkNames=FALSE))
  (checkEquals(round(z["cv"], 1), 0.3, checkNames=FALSE))
  (checkEquals(z["min"],            4, checkNames=FALSE))
  (checkEquals(z["max"],           10, checkNames=FALSE))
  
} ## --- Unit tests end ---

###------------------------------------------------------------------------
### descStat.R ends here
descStat <- structure(
  
  ## --- Function ---
  
  function # Descriptive statistics of numeric variable
  ### A simple wrapper function to get descriptive statistics of numeric variable.
  (
    x,         ##<< numeric, variable
    na.rm=TRUE ##<< logical, remove \code{NA} values
  ) {
    
    if(!is.numeric(x)) stop("'x' must be numeric")
    m <- mean(x, na.rm=na.rm)
    s <- sd(x, na.rm=na.rm)
    c(n=length(x),
      obs=sum(!is.na(x)),
      mean=m,
      median=median(x, na.rm=na.rm),
      sd=s,
      cv=s/m,
      min=min(x, na.rm=na.rm),
      max=max(x, na.rm=na.rm))
    
    ##value<< A vector with:
    ##  \item{n}{length of \code{x}} 
    ##  \item{obs}{number of values in \code{x}, i.e., length without \code{NA} values}
    ##  \item{mean}{mean}
    ##  \item{median}{median}
    ##  \item{sd}{standard deviation}
    ##  \item{cv}{coefficient of variation = sd / mean}
    ##  \item{min}{minimum}
    ##  \item{max}{maximum}
    
  },## --- Function end ---
  
  ## --- Examples ---
  
  ex=function() {
    
    ## Usage for one variable  
    test <- rnorm(n=100, mean=100, sd=10)
    test[c(1, 2, 3)] <- NA
    descStat(x=test)
    
    ## Can be used neatly with the summaryBy function from the doBy package
    if(require(package="doBy")) {
      data(dietox)
      summaryBy(Weight ~ Evit + Cu, data=dietox, FUN=descStat)  
    }
    
  })## --- Examples end ---


