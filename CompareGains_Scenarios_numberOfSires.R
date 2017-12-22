setwd("~/Essentials/")


ped4 <- read.table("PedigreeClass04.txt", header=TRUE)
ped8 <- read.table("PedigreeClass08.txt", header=TRUE)
ped12 <- read.table("PedigreeClass012.txt", header=TRUE)
ped41 <- read.table("Pedigree4_1.txt", header=TRUE)
ped81 <- read.table("Pedigree8_1.txt", header=TRUE)
ped121 <- read.table("Pedigree12_1.txt", header=TRUE)


ped4A <- aggregate(ped4$gvNormUnres1 ~ ped4$Generation, FUN="mean")
colnames(ped4A) <- c("Generation", "Ped4")
ped8A <- aggregate(ped8$gvNormUnres1 ~ ped8$Generation, FUN="mean")
colnames(ped8A) <- c("Generation", "Ped8")
ped12A <- aggregate(ped12$gvNormUnres1 ~ ped12$Generation, FUN="mean")
colnames(ped12A) <- c("Generation", "Ped12")
ped41A <- aggregate(ped41$gvNormUnres1 ~ ped41$Generation, FUN="mean")
colnames(ped41A) <- c("Generation", "Ped4a")
ped81A <- aggregate(ped81$gvNormUnres1 ~ ped81$Generation, FUN="mean")
colnames(ped81A) <- c("Generation", "Ped8a")
ped121A <- aggregate(ped121$gvNormUnres1 ~ ped121$Generation, FUN="mean")
colnames(ped121A) <- c("Generation", "Ped12a")

pedA <- merge(ped4A, ped8A, by="Generation")
pedA <- merge(pedA, ped12A, by="Generation")
pedA <- merge(pedA, ped41A, by="Generation")
pedA <- merge(pedA, ped81A, by="Generation")
pedA <- merge(pedA, ped121A, by="Generation")

library(ggplot2)
library(reshape)

pedAM <- melt(pedA, id.vars = "Generation")
Gainplot <- ggplot(data=pedAM, aes(x=Generation, y=value, group=variable, colour=variable)) + geom_path() + xlab("Mean Genetic Value")


var4 <- read.table("Variance4.txt", header=TRUE)
var8 <- read.table("Variance8.txt", header=TRUE)
var12 <- read.table("Variance12.txt", header=TRUE)

var4 <- var4[var4$QtnModel==1, c("Generation", "AdditGeneticVar1")]
var8 <- var8[var8$QtnModel==1, c("Generation", "AdditGeneticVar1")]
var12 <- var12[var12$QtnModel==1, c("Generation", "AdditGeneticVar1")]

colnames(var4) <- c("Generation", "Var4")
colnames(var8) <- c("Generation", "Var8")
colnames(var12) <- c("Generation", "Var12")

VAR <- merge(var4, var8, by="Generation")
VAR <- merge(VAR, var12, by="Generation")

VARm <- melt(VAR, id.vars = "Generation")
Varplot <- ggplot(data=VARm, aes(x=Generation, y=value, group=variable, colour=variable)) + geom_path() + ylim(0, 10) + geom_smooth() + ylab("Additive genetic variance")

library(Rmisc)
multiplot(Gainplot, Varplot, StGain)




#standardiziraj na generacijo 20
pedA20 <- pedA[pedA$Generation %in% 20:60,]
pedA20$st4 <- (pedA20$Ped4 - pedA20$Ped4[pedA20$Generation == 20]) / var4$Var4[var4$Generation==20]
pedA20$st8 <- (pedA20$Ped8 - pedA20$Ped8[pedA20$Generation == 20]) / var8$Var8[var8$Generation==20]
pedA20$st12 <- (pedA20$Ped12 - pedA20$Ped12[pedA20$Generation == 20]) / var12$Var12[var12$Generation==20]



pedA20st <- pedA20[,c("Generation", "st4", "st8", "st12")]
pedA20stM <- melt(pedA20st, id.vars = "Generation")
StGain <- ggplot(data=pedA20stM, aes(x=Generation, y=value, group=variable, colour=variable)) + geom_path() + ylab("Standardised mean genetic value")
