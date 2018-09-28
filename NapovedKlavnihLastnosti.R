#data za napoved plemenskih vrednosti za klavne lastnosti
ped <- as.data.frame(read.csv("~/Documents/NapovedPV_Klavne/ft_LS_10_17.pedig", header=FALSE, sep=" ", as.is=TRUE, colClasses=c(rep("character",5)))[,c(2)]) #zadnji stolpec do dodane sekvence živali (ponovljeno)
colnames(ped) <- "Info"
ped$Animal <- substr(ped$Info, 1, 6)
ped$Sire <- substr(ped$Info, 7, 12)
ped$Dam <- substr(ped$Info, 13, 18)

dat <- read.csv("~/Documents/NapovedPV_Klavne/ft_LS_10_17_data.csv", header=TRUE)
colnames(dat) <- c("TOPLA_MASA", "PGK_TM",  "PRIRAST", "INT_RASTI",  "KONF", "ZAMASC", "ok2", "ok3", "ok4", "STAROST", "STAROSTk", "BREED", "SPOL", "S_ROJ", "S_ZAKOLA", "L_ZAKOLA", "M_ZAKOLA", "K_L_ZAKOLA", "KLAVNICA", "HERD", "animal")


dat <- read.csv("~/Documents/NapovedPV_Klavne/dat_ft_LS_05_17.txt", sep=";")
dat <- dat[dat$PASMA==2,]
summary(dat$STAROST)
summary(dat$DAT_ROJSTVO)


datModel <- dat[,c("ZIV_ID_SEQ", "CREDA", "TOPLA_MASA", "K_L_ZAKOLA", "M_ZAKOLA", "STAROST", "INT_RASTI","SPOL")]
hist(datModel$STAROST)
table(datModel$SPOL)
summary(datModel$CREDA)
summary(datModel$ZIV_ID_SEQ)
hist(datModel$M_ZAKOLA)
hist(datModel$INT_RASTI)
table(datModel$K_L_ZAKOLA)

ped <- read.table("~/Documents/NapovedPV_Klavne/ped_ft_LS_05_17.txt", header=FALSE)
colnames(ped) <- c("Animal", "Sire", "Dam", "BirthDate", "Sex", "Breed")
ped$Sire[ped$Sire=="9999999"] <- 0
ped$Dam[ped$Dam=="9999999"] <- 0

write.table(ped [, 1:3], "~/Documents/NapovedPV_Klavne/Blupf90.ped", quote=FALSE, col.names = FALSE, row.names=FALSE, sep=" ")
write.table(datModel, "~/Documents/NapovedPV_Klavne/Blupf90.dat", quote=FALSE, col.names = FALSE, row.names=FALSE, sep=" ")
