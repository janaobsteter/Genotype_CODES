

pedigree <- read.table("~/bin/AlphaRelate/example/Example1Pedigree/PedigreeInbreeding.txt")
pedigree <- read.table("~/Documents/PhD/Simulaton/InbreedingCLASS.txt")
pedigree <- read.table("~/Documents/PhD/Simulaton/INBREEDING.txt")
pedigree$V1 <- as.numeric(as.character(pedigree$V1))
pedigree$V2 <- as.numeric(pedigree$V2)
pedigree$V3 <- as.numeric(pedigree$V3)
pedigree$V3 <- as.character(pedigree$V3)


colnames(pedigree) <- c("Indiv", "inbreeding", "Rep")
colnames(pedigree) <- c("Indiv", "inbreeding", "Group")
#split he group into Rep and Scenatio
a <- matrix(unlist(strsplit(pedigree$Group, "[_]")), ncol=2, byrow=TRUE)
pedigree$Rep <- a[,1]
pedigree$Scenario <- a[,2]
pedigree$Indiv <- as.numeric(as.character(pedigree$Indiv))
PED <- pedigree

#pedigree <- pedigree[order(pedigree$Rep, pedigree$Indiv),]
pedigree <- pedigree[!is.na(pedigree$Indiv),]
#pedigree$Gen <- rep(1:60, 8640)
length(intersect(generation$Indiv, pedigree$Indiv))

generation <- read.table("~/bin/AlphaRelate/example/Example1Pedigree/Generation.txt")
generation <- read.table("~/Documents/PhD/Simulaton/GenerationCLASS.txt", header=TRUE)
colnames(generation) <- c("Indiv", "Gen", "Rep")
generation <- generation[generation$Indiv != "Indiv",]
generation$Indiv <- as.numeric(as.character(generation$Indiv))
generation$Rep <- as.numeric(generation$Rep)

#here split the pedigree since it is too large
for (scenario in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
  pedigreeSc <- subset(pedigree, pedigree$Scenario==scenario)
  pedGen <- merge(pedigreeSc, generation, by=c("Indiv"), all.x=TRUE)
}

generation <- generation[,c(1,2)] #rabiš samo indiv in Gen, je po vseh Rep enako
pedGen <- merge(pedigree, generation, by=c("Indiv"), all.x=TRUE)
genInb <- aggregate(pedGen$inbreeding ~ pedGen$Gen, FUN="mean")
colnames(genInb) <- c("Gen", "inbreeding")

pedGen$Rep <- as.numeric(pedGen$Rep)
pedGen$Indiv <- as.character(pedGen$Indiv)
pedGen$Gen <- as.numeric(pedGen$Gen)
ScenarioDF <- data.frame(Scenario="Class", Interval = NA, Rep=NA, Ne=NA)
for (rep in 0:10) {
  for (int in c(0,20,40)) {
    pedGen1 <- pedGen
    pedGen <- pedGen[(pedGen$Rep == rep) & (pedGen$Gen %in% int:(int+20)),]
    genInb <- aggregate(pedGen$inbreeding ~ pedGen$Gen, FUN="mean")
    colnames(genInb) <- c("Gen", "inbreeding")
    genInb$y = log(1 - genInb$inbreeding)
    genInb$t = genInb$Gen - min(genInb$Gen) + 1
    fit = MASS:::rlm(genInb$y ~ genInb$t, maxit=2000)
    # lahko bi uporabilli tudi lm(), ampak je rlm bolj robusten (pri rastlinskih simulacijah vidim precejsnje spremembe skozi cas, mislim, da pri tebi tukaj ne bi smelo biti vecjih razlik med rlm in lm)
    dF = 1 - exp(coef(fit)[2])
    Ne = 1 / (2 * dF)
    pedGen <- pedGen1
    ScenarioDF <- rbind(ScenarioDF, c("Class", int, rep, Ne))
  }
}
ScenarioDF <- ScenarioDF[-1,]
ScenarioDF$Interval <- as.factor(ScenarioDF$Interval)
ScenarioDF$Ne <- as.numeric(ScenarioDF$Ne)
aggregate(ScenarioDF$Ne ~ ScenarioDF$Interval, FUN="mean")

plotList <- list()
for (rep in 0:10) {
  df <- pedGen[pedGen$Rep == rep,]
  dfA <- aggregate(df$inbreeding ~ df $Gen, FUN="mean")
  colnames(dfA) <- c("Gen", "Inbreeding")
  ggplot(data=dfA, aes(x=Gen, y=Inbreeding)) + geom_point() 
}

#1-20
genInb1 <- genInb
genInb <- genInb[genInb$Gen %in% 1:20,]
genInb$y = log(1 - genInb$inbreeding)
genInb$t = genInb$Gen - min(genInb$Gen) + 1
fit = MASS:::rlm(genInb$y ~ genInb$t, maxit=2000)
# lahko bi uporabilli tudi lm(), ampak je rlm bolj robusten (pri rastlinskih simulacijah vidim precejsnje spremembe skozi cas, mislim, da pri tebi tukaj ne bi smelo biti vecjih razlik med rlm in lm)
dF = 1 - exp(coef(fit)[2])
Ne = 1 / (2 * dF)
print(Ne)
genInb <- genInb1

#20-40
genInb1 <- genInb
genInb <- genInb[genInb$Gen %in% 20:40,]
genInb$y = log(1 - genInb$inbreeding)
genInb$t = genInb$Gen - min(genInb$Gen) + 1
fit = MASS:::rlm(genInb$y ~ genInb$t, maxit=2000)
# lahko bi uporabilli tudi lm(), ampak je rlm bolj robusten (pri rastlinskih simulacijah vidim precejsnje spremembe skozi cas, mislim, da pri tebi tukaj ne bi smelo biti vecjih razlik med rlm in lm)
dF = 1 - exp(coef(fit)[2])
Ne = 1 / (2 * dF)
print(Ne)
genInb <- genInb1

#40-60
genInb1 <- genInb
genInb <- genInb[genInb$Gen %in% 40:60,]
genInb$y = log(1 - genInb$inbreeding)
genInb$t = genInb$Gen - min(genInb$Gen) + 1
fit = MASS:::rlm(genInb$y ~ genInb$t, maxit=2000)
# lahko bi uporabilli tudi lm(), ampak je rlm bolj robusten (pri rastlinskih simulacijah vidim precejsnje spremembe skozi cas, mislim, da pri tebi tukaj ne bi smelo biti vecjih razlik med rlm in lm)
dF = 1 - exp(coef(fit)[2])
Ne = 1 / (2 * dF)
print(Ne)
genInb <- genInb1
