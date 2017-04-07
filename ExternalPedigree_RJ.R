#najprej ena burn in generacija v AlphaSimu - obtain 20,000 TBV
##############################################
ped1 <- read.table("~/Documents/PhD/Simulaton/PedigreeAndGeneticValues_1x20000_init.txt", header=T)
#head(ped1)
#keep Ind, FID, MID and gvNormUnres1
ped1 <- ped1[,c(2,3,4,9)]
ped1$sex <- NA #spol
ped1$gen <- NA #starost
ped1$cat <- NA #kategorija

######################################################################################################################
#ŽENSKE
######################################################################################################################
#1:18000 so ženske
ped1$sex[1:18000] <- "F"
#imamo 15% populacije starost = 0, 14.5% starost 1, 8% starost 2, potem pa 62.5% krav - starost 3, 4, 5, 6, 7 vsaka 1/5
#izmed krav starost 3 leta vsako leto odberemo remont BM (najboljše)

###################################################
#dopiši starosti živalim
gen0 <- 0.15 #pod 1 letom
gen1 <- 0.145 # 1 - 2 leta
gen2 <- 0.08 # 2 - 3 leta
genK <- 0.625 * 0.2 # starejše (vsaka gen 1/5)

nF <- 18000 #celotno število ženskih živali
n <- 1
ped1$gen[n:(gen0*nF)] <- 0 
n <- n + gen0*nF
ped1$gen[n:(n+gen1*nF-1)] <- 1 
n <- n + gen1*nF
ped1$gen[n:(n+gen2*nF-1)] <- 2
n <- n + gen2*nF
ped1$gen[n:(n+genK*nF-1)] <- 3
n <- n + genK*nF
ped1$gen[n:(n+genK*nF-1)] <- 4
n <- n + genK*nF
ped1$gen[n:(n+genK*nF-1)] <- 5
n <- n + genK*nF
ped1$gen[n:(n+genK*nF-1)] <- 6
n <- n + genK*nF
ped1$gen[n:(n+genK*nF-1)] <- 7

table(ped1$gen)
sum(table(ped1$gen))
#ped1$gen <- as.factor(ped1$gen)


##########################################################################
#dopiši kategorije
t0 <- 0.15
t12 <- 0.11
pt12 <- 0.035
t24 <- 0.03
pt24 <- 0.05
k <- 0.625 #zaenkrat pusti vse krave --> kasneje najboljšim starosti 3, 4 in 5 let spremeni kategorijo v BM

b0 <- 0.05
b12 <- 0.03
b24 <- 0.005
pb <- 0.001

sum(t0, t12, pt12, t24, pt24, k)

nF <- 18000 #celotno število ženskih živali
n <- 1
ped1$cat[n:(t0*nF)] <- "T0"
n <- n + t0*nF
ped1$cat[n:(n+t12*nF-1)] <- "T12"
n <- n + t12*nF
ped1$cat[n:(n+pt12*nF-1)] <- "PT12"
n <- n + pt12*nF
ped1$cat[n:(n+t24*nF-1)] <- "T24"
n <- n + t24*nF
ped1$cat[n:(n+pt24*nF-1)] <- "PT24"
n <- n + pt24*nF
ped1$cat[n:(n+k*nF-1)] <- "K"

table(ped1$cat)
sum(table(ped1$gen))
table(ped1$gen, ped1$cat) #TUKAJ POGLEJ, če ti štimajo številke! Vsaka kategorijo samo eno leto, razen krave 5 let

######################################################################################################
#MOŠKI
######################################################################################################
#18000:20000 so moški
ped1$sex[18001:20000] <- "M"
#imamo 55% populacije starost = 0, 40% starost 1, 5% starost 2 = aktivna populacija (skupaj 99%)
#1% bikov pa je iz AI - niso v aktivni populaciji, se pa uporabljajo za osemenjevanje --> morajo pa biti v pedigreju
#plus še 1.2% čakajočih bikov iz 3,4,5 let nazaj

###################################################
"""#dopiši starosti živalim
gen0 <- 0.5 #pod 1 letom
gen1 <- 0.145 # 1 - 2 leta
gen2 <- 0.08 # 2 - 3 leta
genK <- 0.625 * 0.2 # starejše (vsaka gen 1/5)

nF <- 18000 #celotno število ženskih živali
n <- 1
ped1$gen[n:(gen0*nF)] <- 0 
n <- n + gen0*nF
ped1$gen[n:(n+gen1*nF-1)] <- 1 
n <- n + gen1*nF
ped1$gen[n:(n+gen2*nF-1)] <- 2
n <- n + gen2*nF
ped1$gen[n:(n+genK*nF-1)] <- 3
n <- n + genK*nF
ped1$gen[n:(n+genK*nF-1)] <- 4
n <- n + genK*nF
ped1$gen[n:(n+genK*nF-1)] <- 5
n <- n + genK*nF
ped1$gen[n:(n+genK*nF-1)] <- 6
n <- n + genK*nF
ped1$gen[n:(n+genK*nF-1)] <- 7

table(ped1$gen)
sum(table(ped1$gen))
#ped1$gen <- as.factor(ped1$gen)


##########################################################################
#dopiši kategorije
t0 <- 0.15
t12 <- 0.11
pt12 <- 0.035
t24 <- 0.03
pt24 <- 0.05
k <- 0.625 #zaenkrat pusti vse krave --> kasneje najboljšim starosti 3, 4 in 5 let spremeni kategorijo v BM

b0 <- 0.05
b12 <- 0.03
b24 <- 0.005
pb <- 0.001

sum(t0, t12, pt12, t24, pt24, k)

nF <- 18000 #celotno število ženskih živali
n <- 1
ped1$cat[n:(t0*nF)] <- "T0"
n <- n + t0*nF
ped1$cat[n:(n+t12*nF-1)] <- "T12"
n <- n + t12*nF
ped1$cat[n:(n+pt12*nF-1)] <- "PT12"
n <- n + pt12*nF
ped1$cat[n:(n+t24*nF-1)] <- "T24"
n <- n + t24*nF
ped1$cat[n:(n+pt24*nF-1)] <- "PT24"
n <- n + pt24*nF
ped1$cat[n:(n+k*nF-1)] <- "K"

table(ped1$cat)
sum(table(ped1$gen))
table(ped1$gen, ped1$cat) #TUKAJ POGLEJ, če ti štimajo številke! Vsaka kategorijo samo eno leto, razen krave 5 let
"""


###################################################################################################################3
###################################################################################################################3
###################################################################################################################3
###################################################################################################################3
#verjetno bolje, da štartaš iz nule - torej novorojene populacije
#nasimuliraj samo 2800 (F) + 1100 (M) novorojenih telet = 3900

#kaj narediš z novorojeno generacijo:
  #female: izloči 3.3%, osemeni 73.3%, prenesi 23.3.% telic
  #males: odberi najboljših 8 bikov za mlade, ostali iz testa za pripust, izloči 27% in prenesi 73% bikov

ped1 <- read.table("~/Documents/PhD/Simulaton/Pedigrees/Pedigree_10burnIn_Gen1.txt", header=T)
ped1 <- ped1[,c(1,2,3,4,9)]
nrow(ped1)
ped1$EBV <- exact_corr(ped1$gvNormUnres1, r=0.8) #pridobi EBV iz TBV
cor(ped1$gvNormUnres1, ped1$EBV)
write.csv(ped1, "~/Documents/PhD/Simulaton/Pedigrees/PedPython.txt", quote=F, row.names=F)

ngen <- 0
nF <- 2700
nM <- 1100
####################################################
#tu zdj pripraviš pedigree za naslednjo generacijo
ped1$sex <- NA
ped1$gen <- ngen
ped1$cat <- NA
ped1$age <- max(ped1$gen) - ped1$gen

ped1$sex[1:2700] <- "F"
ped1$sex[2701:3800] <- "M"

##################################################################################################
#ŽENSKE
##################################################################################################
#novorojena teleta - izloči 3.3% najslabpih po PV
#obrejo 73.3% telic, 23.3.% ostane telice
pedF <- subset(ped1, ped1$sex=="F")
pedF <- pedF[order(-pedF$EBV),]

pt12 <- nF*0.73333333333333333333# % novorojenih telic, ki jih do 1 leta osemeniš
tel12 <- nF*0.2333 # % novorojenih telic, ki jih do 1 leta osemeniš

nrow(pedF)
n <- 1
pedF$cat[1:pt12] <- "pt12"
n <- n + pt12 + 1 #NE VEM; ZAKAJ TUKAJ PLUS 1
pedF$cat[n:(n+tel12)] <- "t12"
n <- n + tel12 - 1 #NE VEM; ZAKAJ TUKAJ PLUS 1
n
pedF$cat[n:nF+1] <- "izl"
table(pedF$cat)
sum(table(pedF$cat))

pedF$active <- NA # tukaj še označi, ali so živali aktivne = 1 ali ne =2 (tudi kategorija "izl")
n <- 1
pedF$active[n:(pt12+tel12+1)] <- 1
n <- n + (pt12+tel12) 
n
pedF$active[(n+1):(nF+1)] <- 2
table(pedF$active)
table(pedF$cat,pedF$active) #vse izločene morajo biti imeti active=2 (check)

##########################################################################
#moški
###########################################################################
pedM <- subset(ped1, ped1$sex=="M")
nrow(pedM)
pedM <- pedM[order(pedM$EBV),] #sortiraj po EBV
pedM$cat <- "izl"
pedM$cat[1:8] <- "mladi" #najboljših 8 postane mladih bikov
pedM$cat[9:27] <- "pripust" #ostalih 19 iz testa (27 v testu) gre za pripust
pedM$cat[c(sample(28:1100, 783, replace=F))] <- "biki12" # tukaj random odberi bike, ki preživijo do naslednjega leta
table(pedM$cat)

#naštimaj še aktivnost živali
pedM$active <- 2
pedM$active[pedM$cat=="biki12"] <- 1
pedM$active[pedM$cat=="mladi"] <- 1
pedM$active[pedM$cat=="pripust"] <- 1
table(pedM$cat)
sum(table(pedM$cat))
table(pedM$active)
sum(table(pedM$active))
table(pedM$cat,pedM$active)

##########################################################################
#ustvari še novo generacijo
#na tem mestu še nimaš nobenih znanih staršev! --> naslednje leto bodo mame že znane
ngen <- ngen + 1
pedN <- data.frame(Indiv = 3801:(3800*2), Father = 0, Mother = 0, sex=rep(c("F", "M"), c(2700, 1100)), 
                   gen = ngen,  cat = "nr", active = 1)
table(pedN$sex)
table(pedN$cat)
table(pedN$active)

ped1 <- rbind(pedF, pedM)
ped1 <- ped1[,c(1,2,3,6,7,8,10)]
ped1 <- rbind(ped1, pedN)
ped1$age <- max(ped1$gen) - ped1$gen
table(ped1$age)
sum(ped1$active==1)

write.csv(ped1, "~/Documents/PhD/Simulaton/Pedigrees/External_PedigreeGen1.txt", row.names=F, quote=F)
write.csv(pedN[,1:3], "~/Documents/PhD/Simulaton/Pedigrees/External_PedigreeGen1.txt", row.names=F, quote=F)
write.table(pedN[,1:3], "~/Documents/PhD/Simulaton/Pedigrees/External_PedigreeGen1.txt", row.names=F, quote=F, col.names=F, sep=",")
write.table(pedN[,1:3], "~/bin/AlphaSim1.05Linux/External_PedigreeGen1.txt", row.names=F, quote=F, col.names=F, sep=" ")
#################################################################################################################
#################################################################################################################
#################################################################################################################