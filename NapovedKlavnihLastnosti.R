#data za napoved plemenskih vrednosti za klavne lastnosti
ped <- as.data.frame(read.csv("~/Documents/NapovedPV_Klavne/ft_LS_10_17.pedig", header=FALSE, sep=" ", as.is=TRUE, colClasses=c(rep("character",5)))[,c(2)]) #zadnji stolpec do dodane sekvence živali (ponovljeno)
colnames(ped) <- "Info"
ped$Animal <- substr(ped$Info, 1, 6)
ped$Sire <- substr(ped$Info, 7, 12)
ped$Dam <- substr(ped$Info, 13, 18)

dat <- read.csv("~/Documents/NapovedPV_Klavne/ft_LS_10_17_data.csv", header=TRUE)
colnames(dat) <- c("TOPLA_MASA", "PGK_TM",  "PRIRAST", "INT_RASTI",  "KONF", "ZAMASC", "ok2", "ok3", "ok4", "STAROST", "STAROSTk", "BREED", "SPOL", "S_ROJ", "S_ZAKOLA", "L_ZAKOLA", "M_ZAKOLA", "K_L_ZAKOLA", "KLAVNICA", "HERD", "animal")
