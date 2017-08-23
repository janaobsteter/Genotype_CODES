ped = read.table('~/bin/AlphaSim1.05Linux/REAL20GenSel_GenBM_UpdatedRef//SimulatedData/PedigreeAndGeneticValues_cat.txt', header=T)
sex = read.table('~/bin/AlphaSim1.05Linux/REAL20GenSel_GenBM_UpdatedRef//SimulatedData/Gender.txt', header=T)
ped = read.table('~/bin/AlphaSim1.05Linux//REAL20GenSel_Gen/SimulatedData/PedigreeAndGeneticValues_cat.txt', header=T)
sex = read.table('~/bin/AlphaSim1.05Linux/REAL20GenSel_Gen/SimulatedData/Gender.txt', header=T)
pedsex <- merge(ped, sex, by='Indiv')
#prvih 20 generacij mora štimat, ker jih vleče iz burn-ina
pedsex20 <- pedsex[pedsex$Generation.x %in% 1:20,]
table(pedsex20$sex, pedsex20$Gender)
length(intersect(ped$Father, ped$Mother))

females = ped[ped$sex=='F',]
males = ped[ped$sex=='M',]
mothers = unique(ped$Mother)
fathers = unique(ped$Father)

length(intersect(males$Indiv, mothers))
length(intersect(females$Indiv, fathers))


##########################################################################################
#tukaj preveri starševstva
#########################################################################################
#splotaj število potomcev po očetu
father <- as.data.frame(table(ped$Father))
colnames(father) <- c("Indiv", "NoOff" )
father <- merge(father, ped, by="Indiv")
genFather <- aggregate(father$NoOff[father$cat=='gpb'] ~ father$Generation[father$cat=='gpb'], FUN="mean")
pbFather <- aggregate(father$NoOff[father$cat=='pripust1'] ~ father$Generation[father$cat=='pripust1'], FUN="mean")
colnames(genFather) <- c("Gen", "NoOff" )
colnames(pbFather) <- c("Gen", "NoOff" )
plot(genFather$NoOff ~genFather$Gen)
plot(pbFather$NoOff ~pbFather$Gen)
plot(father$NoOff[father$cat=='gpb'] ~father$Generation[father$cat=='gpb'])
plot(father$NoOff[father$cat=='pb'] ~father$Generation[father$cat=='pb'])

hist(as.data.frame(table(ped$Father))$Freq, breaks=100)
mean(table(ped$Father))
tail(table(ped$Father, ped$cat))


#splotaj število potomcev po mami
mother <- as.data.frame(table(ped$Mother))
hist(as.data.frame(table(ped$Mother))$Freq)
mean(table(ped$Mother))
tail(table(ped$Mother, ped$cat))

#preveri povprečno število potomcev po kategoriji očeta (in generaciji)
mean = 0
number = 0
for (father in intersect(ped$Father, ped$Indiv[ped$cat == 'gpb'])) {
  mean = mean + length(ped$Indiv[ped$Father==father])
  number = number + 1
}
mean / number


#na koliko potomcih so testirani pb
unique(ped$Generation[ped$cat=='gpb'])
mean = 0
for (pb in ped$Indiv[ped$cat=='gpb']) {
  Generation = ped$Generation[ped$Indiv == pb]
  mean = mean + length(ped$Indiv[ped$Father==pb & ped$Generation %in% Generation: (Generation+5)])
  print(c(pb, length(ped$Indiv[ped$Father==pb & ped$Generation %in% Generation: (Generation+5)]), Generation))
}
#preveri povprečno število potomcev po kategoriji matere (in generaciji)
mean = 0
number = 0
for (mother in intersect(ped$Mother, ped$Indiv[ped$cat == 'bm' ])) {#& ped$Generation==56
  #print(length(ped$Indiv[ped$Mother==mother]))
  #print(mother)
  mean = mean + length(ped$Indiv[ped$Mother==mother])
  number = number + 1
}
mean / number

#preveri kategorije mam potomcevNP
ped$cat[ped$Indiv %in% ped$Mother[ped$Indiv[ped$cat=='potomciNP']]] #več pBM kot pa potomcevNP - niso vse osemenitve uspešne

#preveri očete nr
table(ped$cat[ped$Indiv %in% ped$Father[ped$cat=='nr']]) #kategorije
table(ped$Father[ped$Father[ped$cat=='nr']]) #zastopanost
sum(table(ped$Father[ped$Father[ped$cat=='nr']]))
#preveri očete potomcev NP
table(ped$cat[ped$Indiv %in% ped$Father[ped$cat=='potomciNP']])
table(ped$Father[ped$Father[ped$cat=='potomciNP']])
sum(table(ped$Father[ped$Father[ped$cat=='potomciNP']]))



################################################
#preveri točnosti napovedi
sol <- read.table('~/bin/AlphaSim1.05Linux//renumbered_Solutions_42')
ped$EBV <- sol[,3]
cat <- 'gpb'
cor(ped$EBV[ped$cat==cat], ped$gvNormUnres1[ped$cat==cat])
genMean <- aggregate(ped$gvNormUnres1 ~ ped$Generation, FUN='mean')
genMean <- genMean[genMean$`ped$Generation`[41:60],]
plot(genMean$`ped$gvNormUnres1` ~ genMean$`ped$Generation`, type='line')


##########################3
#preveri intersect med gpb in pb -should be none!
ped = read.table('~/bin/AlphaSim1.05Linux/REAL20GenSel_Gen/SimulatedData/PedigreeAndGeneticValues_cat.txt', header=T)
gpb = ped$Indiv[ped$cat=='gpb']
pb = ped$Indiv[ped$cat=='pb']
intersect(gpb, pb)
length(gpb)
length(pb)



################################################
gpb <- ped[ped$cat=='gpb',]
gF <- subset(father)