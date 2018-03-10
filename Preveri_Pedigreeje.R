ped <- read.csv("~/Class_20Rep.txt",sep=" ")
ped <- read.csv("~/Class_1Year.txt",sep=" ")



pedL <- ped[ped$Generation %in% 40:60,]

table(pedL$Father, pedL$Generation)

ped50 <- ped[ped$Generation==50,]
table(ped$Generation[ped$Indiv %in% ped50$Father])
table(ped50$Father)

FatherUse <- data.frame()
for (father in unique(pedL$Father)) {
  range <- max(ped$Generation[ped$Father == father]) - min(ped$Generation[ped$Father == father])
  FatherUse <- rbind(FatherUse, c(father, range))
}
colnames(FatherUse) <- c("Father", "NoOffspring")

FatherUse$NoOffspring <- as.numeric(FatherUse$NoOffspring)
hist(FatherUse$NoOffspring)
table(FatherUse$NoOffspring)


ped[ped$Indiv %in% FatherUse$Father[FatherUse$NoOffspring > 5],]
table(ped$cat[ped$Indiv %in% FatherUse$Father[FatherUse$NoOffspring == 5]])
table(ped$Generation[ped$Indiv %in% FatherUse$Father[FatherUse$NoOffspring == 5]])

ped[ped$Indiv %in% FatherUse$Father[FatherUse$NoOffspring == 0],]
table(ped$cat[ped$Indiv %in% FatherUse$Father[FatherUse$NoOffspring == 0]])
table(ped$Generation[ped$Indiv %in% FatherUse$Father[FatherUse$NoOffspring == 0]])


#
#tukaj pa rabiš preverit, koliko očetov je odbranih na eno generacijo
#še vedno
ped <- read.csv("~/Class_1Pb.txt",sep=" ")
ped <- read.csv("~/Gen_1Pb.txt",sep=" ")
ped <- read.csv("~/GenSLO_1Pb.txt",sep=" ")
ped <- read.csv("~/BmGen_1Pb.txt",sep=" ")

table(ped$Generation[ped$cat == 'pb'])
table(ped$Generation[ped$cat == 'mladi'])


