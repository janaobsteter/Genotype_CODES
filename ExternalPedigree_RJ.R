externalPed <- as.data.frame(matrix(ncol=3, nrow=25000))
colnames(externalPed) <- c('ID', 'MID', 'FID')
externalPed$Sex <- NA
externalPed$Sex[1:3000] <- "F"
externalPed$Sex[3001:6000] <- "M"

krave <- sample(1:3000, 500)
externalPed$kat <- NA
externalPed$kat[krave] <- "K"

while (n > ngen) {
  
}

ped <- as.data.frame(matrix(ncol=3, nrow=19000))
ped$Sex <- NA
ped$Sex[1:17500] <- "F"
ped$Sex[17501:19000] <- "M"

ped$kat <- NA
numbersF <- 1:17500

krave <- sample(numbersF, 10000)
numbersF <-  setdiff(numbersF, krave)
length(numbersF)

nr <- sample(numbersF, 3000)
numbersF <-  setdiff(numbersF, nr)
length(numbersF)

tel12 <- sample(numbersF, 2500)
numbersF <-  setdiff(numbersF, tel12)
length(numbersF)

tel24 <-sample(numbersF, 600)
numbersF <-  setdiff(numbersF, tel24)
length(numbersF)

pt12 <- sample(numbersF, 500)
numbersF <-  setdiff(numbersF, pt12)
length(numbersF)

pt24 <-sample(numbersF, 900)
numbersF <-  setdiff(numbersF, pt24)
length(numbersF)

numbersM <- 1:1500

bikinr <- sample(numbersM, 800)
numbersM <- setdiff(numbersM, bikinr)

biki12 <- sample(numbersM, 585)
numbersM <- setdiff(numbersM, biki12)

biki24 <- sample(numbersM, 100)
numbersM <- setdiff(numbersM, biki24)

pb <- sample(numbersM, 15)
numbersM <- setdiff(numbersM, pb)
