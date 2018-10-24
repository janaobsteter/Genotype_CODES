#join PedTotals into GenPed

ped1 <- read.csv("~/PedTotalhome.txt")
ped2 <- read.csv("~/PedTotalimport.txt")
#table(ped1$Generation)
#table(ped2$Generation)
#tail(ped1)
#tail(ped2)
#summary(ped1$Indiv[ped1$Generation==21])
#summary(ped2$Indiv[ped2$Generation==21])
ped1$gvNormUnres1[ped1$Generation==max(ped1$Generation)] <- runif(8640, 0, 5)
ped2$gvNormUnres1[ped2$Generation==max(ped1$Generation)] <- runif(8640, 0, 5)
ped1$EBV[ped1$Generation==21] <- runif(8640, 0, 5)
ped2$EBV[ped2$Generation==21] <- runif(8640, 0, 5)

#length(intersect(ped1$Indiv, ped2$Indiv))
#intersect(ped1$Indiv, ped2$Indiv)


genped <- read.csv("~/GenPed_EBV.txt")
#genpedH <- read.csv("~/GenPed_EBVhome.txt")
#genpedI <- read.csv("~/GenPed_EBVimport.txt")
#table(genpedH$Generation)
#table(genpedI$Generation)
##tail(genpedH)
#tail(genpedI)

ped1 <- ped1[,colnames(genped)]
ped2 <- ped2[,colnames(genped)]

write.csv(rbind(ped1, ped2), "/home/jana/GenPed_EBV.txt", quote=FALSE, row.names=FALSE)
