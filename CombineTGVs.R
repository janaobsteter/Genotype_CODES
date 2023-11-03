tgvs1 <- read.csv("/home/v1jobste/JanaO/TGVsAll_10KRef_SU52_06032019.csv", sep=" ")
tgvs2 <- read.csv("/home/v1jobste/JanaO/TGVsAll_10KRef_SU53_06032019.csv", sep=" ")
tgvs3 <- read.csv("/home/v1jobste/JanaO/TGVsAll_10KRef_SU54_06032019.csv", sep=" ")

TGV <- rbind(tgvs1, tgvs2)
TGV <- rbind(TGV, tgvs3)

write.csv(TGV, "TGVSALL_SUx_06032019.csv", quote=FALSE, row.names=FALSE)
print(paste0("Created data frame TGVSALL_SUx_06032019.csv"))
