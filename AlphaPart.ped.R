### pedS.R
###-----------------------------------------------------------------------------
### $Id$
###-----------------------------------------------------------------------------

### This is a copy of partAGV_stohastic demo code (only the data part)

## To keep environment clean
AlphaPart.ped.simulation <- function()
{

##install.packages(pkg=c("truncnorm"), dep=TRUE)
library(package="truncnorm")

### EXAMPLE PEDIGREE
###-----------------------------------------------------------------------------

## Generation 0
 id0 <- c("01", "02", "03", "04", "05", "06", "07", "08")
fid0 <- mid0 <- rep(NA, times=length(id0))
  h0 <- rep(c(1, 2), each=4)
  g0 <- rep(0, times=length(id0))

## Generation 1
 id1 <- c("11", "12", "13", "14", "15", "16", "17", "18")
fid1 <- c("02", "02", "02", "02", "06", "06", "06", "06")
mid1 <- c("01", "01", "03", "04", "05", "05", "07", "08")
  h1 <- h0
  g1 <- rep(1, times=length(id1))

## Generation 2
 id2 <- c("21", "22", "23", "24", "25", "26", "27", "28")
fid2 <- c("13", "13", "13", "13", "13", "13", "13", "13")
mid2 <- c("11", "12", "14", "14", "15", "16", "17", "18")
  h2 <- h0
  g2 <- rep(2, times=length(id2))

## Generation 3
 id3 <- c("31", "32", "33", "34", "35", "36", "37", "38")
fid3 <- c("24", "24", "24", "24", "24", "24", "24", "24")
mid3 <- c("21", "21", "22", "23", "25", "26", "27", "28")
  h3 <- h0
  g3 <- rep(3, times=length(id3))

## Generation 4
 id4 <- c("41", "42", "43", "44", "45", "46", "47", "48")
fid4 <- c("34", "34", "34", "34", "34", "34", "34", "34")
mid4 <- c("31", "32", "32", "33", "35", "36", "37", "38")
  h4 <- h0
  g4 <- rep(4, times=length(id4))

## Generation 5
 id5 <- c("51", "52", "53", "54", "55", "56", "57", "58")
fid5 <- c("44", "44", "44", "44", "44", "44", "44", "44")
mid5 <- c("41", "42", "43", "43", "45", "46", "47", "48")
  h5 <- h0
  g5 <- rep(5, times=length(id4))

ped <- data.frame( id=c( id0,  id1,  id2,  id3,  id4,  id5),
                  fid=c(fid0, fid1, fid2, fid3, fid4, fid5),
                  mid=c(mid0, mid1, mid2, mid3, mid4, mid5),
                  loc=c(  h0,   h1,   h2,   h3,   h4,  h5),
                  gen=c(  g0,   g1,   g2,   g3,   g4,  g5))
ped$sex <- 2
ped[ped$id %in% ped$fid, "sex"] <- 1
ped$loc.gen <- with(ped, paste(loc, gen, sep="-"))

### SIMULATE ADDITIVE GENETIC VALUES - STOHASTIC
###-----------------------------------------------------------------------------

## --- Parameters of simulation ---

## Additive genetic mean in founders by location
mu1 <- 2
mu2 <- 0

## Additive genetic variance in population
sigma2 <- 1
sigma  <- sqrt(sigma2)

## Threshold value for Mendelian sampling for selection - only values above this
##  will be accepted in simulation
t <- 0

## Set seed for simulation
set.seed(seed=19791123)

## --- Start of simulation ---

ped$agv1 <- NA ## Scenario (trait) 1: No selection in the second location
ped$agv2 <- NA ## Scenario (trait) 2:    Selection in the second location

## Generation 0  - founders (for simplicity set their values to the mean of location)
ped[ped$gen == 0 & ped$loc == 1, c("agv1", "agv2")] <- mu1
ped[ped$gen == 0 & ped$loc == 2, c("agv1", "agv2")] <- mu2

## Generation 1+ - non-founders
for(i in (length(g0)+1):nrow(ped)) { 
  ## Trait 1: measured in both locations, sample MST according to location - loc 1: only positive MST
  #agv1 - sample MST
  if(ped[i, "loc"] == 1) {
    w <- rtruncnorm(n=1, mean=0, sd=sqrt(sigma2/2), a=t)
  } else {
    w <- rnorm(n=1, mean=0, sd=sqrt(sigma2/2))
  }
  #compute AGV1 as PA + MST
  ped[i, "agv1"] <- round(0.5 * ped[ped$id %in% ped[i, "fid"], "agv1"] +
                          0.5 * ped[ped$id %in% ped[i, "mid"], "agv1"] +
                          w, digits=1)
  ## Trait 2: measured only in second location 
  #agv2
  if(ped[i, "loc"] == 2) {
    w <- rtruncnorm(n=1, mean=0, sd=sqrt(sigma2/2), a=t)
  } ## for location 1 take the same values as above
  ped[i, "agv2"] <- round(0.5 * ped[ped$id %in% ped[i, "fid"], "agv2"] +
                          0.5 * ped[ped$id %in% ped[i, "mid"], "agv2"] +
                          w, digits=1)
}

  ## --- Return ---
  ped

}

AlphaPart.ped <- AlphaPart.ped.simulation()

###-----------------------------------------------------------------------------
### pedS.R ends here
