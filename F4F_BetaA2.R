betaPed <- read.table("~/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/IDBv03/All_BetaA2.ped")
betaPed <- betaPed[,-c(3,4,5,6)]

betaPed$SNP1 <- paste0(betaPed$V7, betaPed$V8)
betaPed$SNP2 <- paste0(betaPed$V9, betaPed$V10)
betaPed$SNP3 <- paste0(betaPed$V11, betaPed$V12)
betaPed$SNP4 <- paste0(betaPed$V13, betaPed$V14)
betaPed$SNP5 <- paste0(betaPed$V15, betaPed$V16)


betaPed[betaPed$SNP1 != betaPed$SNP2]
betaPed[betaPed$SNP2 != betaPed$SNP3]
betaPed[betaPed$SNP3 != betaPed$SNP4]
betaPed[betaPed$SNP4 != betaPed$SNP5]

########################################################
#kappa kazein
kapaAB_1 <- read.table("~/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/IDBv03/KapaAB_1.ped")
kapaAB_1 <- kapaAB_1[,-c(3,4,5,6)]
kapaAB_1$SNP1 <- paste0(kapaAB_1$V7, kapaAB_1$V8)
kapaAB_1$SNP2 <- paste0(kapaAB_1$V9, kapaAB_1$V10)
kapaAB_1$SNP3 <- paste0(kapaAB_1$V11, kapaAB_1$V12)
kapaAB_1$SNP4 <- paste0(kapaAB_1$V13, kapaAB_1$V14)

kapaAB_1[kapaAB_1$SNP1 != kapaAB_1$SNP2,]
kapaAB_1[kapaAB_1$SNP2 != kapaAB_1$SNP3,]
kapaAB_1[kapaAB_1$SNP3 != kapaAB_1$SNP4,]
kapaAB_1[kapaAB_1$SNP1 != kapaAB_1$SNP2,]


##
kapaAB_2 <- read.table("~/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/IDBv03/KapaAB_2.ped")
kapaAB_2 <- kapaAB_2[,-c(3,4,5,6)]
x <- 7
y <- 8
c <- ncol(kapaAB_2) +1
for (i in 1:8) {
  kapaAB_2[,c] <- paste0(kapaAB_2[,x], kapaAB_2[,y])
  x <- x+2
  y <- y+2
  c <- c+1
}

all.equal(kapaAB_2[,23] , kapaAB_2[,24] , kapaAB_2[,25])
