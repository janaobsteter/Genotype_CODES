m1 <- read.csv("MSTCat_SU55.csv")
m1$Strategy <- "SU55"
m2 <- read.csv("MSTCat_SU51.csv")[1541:3080,]
m2$Strategy <- "SU51"
m3 <- read.csv("MSTCat_SU15.csv")[3081:4620,]
m3$Strategy <- "SU15"

M <- rbind(m1, m2)
M <- rbind(M, m3)
print("Writing MSTCat.csv file")
write.csv(M, "MSTCat.csv", quote=FALSE, row.names=FALSE)



m1 <- read.csv("MST_SU55.csv")
m1$Strategy <- "SU55"
m2 <- read.csv("MST_SU51.csv")[2101:4200,]
m2$Strategy <- "SU51"
m3 <- read.csv("MST_SU15.csv")[4201:6300,]
m3$Strategy <- "SU15"

M <- rbind(m1, m2)
M <- rbind(M, m3)
print("Writing MST.csv file")
write.csv(M, "MST.csv",	quote=FALSE, row.names=FALSE)

