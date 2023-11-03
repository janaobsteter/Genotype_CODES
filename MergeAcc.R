acc15 <- read.csv("Accuracy_CatAgeSU15.csv")
acc51 <- read.csv("Accuracy_CatAgeSU51.csv")
acc55 <- read.csv("Accuracy_CatAgeSU55.csv")

ACCAge <- rbind(acc15, acc51)
ACCAge <- rbind(ACCAge, acc55)

write.csv(ACCAge, "ACCAge.csv")
acc15 <- read.csv("Accuracy_CatSU15.csv")
acc51 <- read.csv("Accuracy_CatSU51.csv")
acc55 <- read.csv("Accuracy_CatSU55.csv")

ACC <- rbind(acc15, acc51)
ACC <- rbind(ACC, acc55)
write.csv(ACC, "ACC.csv")

