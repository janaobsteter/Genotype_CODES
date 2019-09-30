library(AlphaSimR)

founderPop = runMacs(nInd=1000, nChr=10, segSites=1000)
SP = SimParam$new(founderPop)
SP$addTraitA(nQtlPerChr = 1000, mean = 0, var = 1)
SP$setGender("yes_sys")

pop = newPop(founderPop)
popMean = popVar = numeric(20)

for (cycle in 1:20) {
  pop = selectCross(pop=pop, nFemale = 500, nMale = 25, use = "gv", nCrosses = 1000)
  popMean[cycle] = meanG(pop)
  popVar[cycle] = varG(pop)
}


library(AlphaSimR)
library(ggplot2)

nChr = 10
nQtl = 200
nSnp = 500

nFounders = 50
nParents = 50
nCrosses = 100
nProgeny = 10
breedingCycles = 10

initMeanG = 10
initVarG = 1
h2 = 0.4

founderPop = runMacs(nInd = nFounders, nChr = nChr, segSites = nQtl + nSnp)


SP = SimParam$new(founderPop)
SP$addTraitA(nQtlPerChr = nQtl, mean = initMeanG, var = initVarG)
SP$setVarE(h2=h2)
SP$addSnpChip(nSnp)
output = data.frame(cycle=0:breedingCycles,
                    mean = numeric(breedingCycles+1),
                    var = numeric(breedingCycles +1 ))

initPop = newPop(founderPop)

initPop = randCross(initPop, nCrosses = nCrosses, nProgeny = nProgeny)

pop = initPop
