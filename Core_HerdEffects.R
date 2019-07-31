ped <- read.csv("~/PedPlusEffects.csv")
pedC <- ped[!is.na(ped$PhenoPE) & !is.na(ped$PhenoPEHY),]

ped 
library(dplyr)

#COR TGV - 
gvEbv <- pedC %>%
  group_by(Generation) %>%
  summarize(COR=cor(gvNormUnres1,EBV))
colnames(gvEbv) <- c("Generation", "gvEbv")

gvEbvPE <- pedC %>%
  group_by(Generation) %>%
  summarize(COR=cor(gvNormUnres1,EBVPE))
colnames(gvEbvPE) <- c("Generation", "gvEbvPE")

gvEbvPEHY <- pedC %>%
  group_by(Generation) %>%
  summarize(COR=cor(gvNormUnres1,EBVPEHY))
colnames(gvEbvPEHY) <- c("Generation", "gvEbvPEHY")

gvCor <- merge(gvEbv, gvEbvPE, by="Generation")
gvCor <- merge(gvCor, gvEbvPEHY, by="Generation")

gvCorm <- melt(gvCor, id.vars = "Generation")
qplot(data=gvCorm, x=Generation, y=value, colour=variable, group=variable, geom="point")

#Cor Pheno - 
phenoEbv <- pedC %>%
  group_by(Generation) %>%
  summarize(COR=cor(phenoNormUnres1,EBV))
colnames(phenoEbv)[2] <- "Pheno-EBV"

phenoEbvPE <- pedC %>%
  group_by(Generation) %>%
  summarize(COR=cor(phenoNormUnres1,EBVPE))
colnames(phenoEbvPE)[2] <- "Pheno-EbvPE"

phenoEbvPEHY <- pedC %>%
  group_by(Generation) %>%
  summarize(COR=cor(phenoNormUnres1,EBVPEHY))
colnames(phenoEbvPEHY)[2] <- "Pheno-EbvPEHY"

phenoCor <- merge(phenoEbv, phenoEbvPE, by="Generation")
phenoCor <- merge(phenoCor, phenoEbvPEHY, by="Generation")

head(phenoCor)
phenoCorm <- melt(phenoCor, id.vars = "Generation")
qplot(data=phenoCorm, x=Generation, y=value, colour=variable, group=variable, geom="point")


pedC$PhenominPE <- pedC$phenoNormUnres1 - pedC$permEnv_est
cor(pedC$phenoNormUnres1, pedC$EBV)
cor(pedC$PhenominPE, pedC$EBV)

cor(ped$HerdYEff_est, ped$HerdYEff, use="pairwise.complete.obs")
cor(ped$permEnv.x, ped$permEnv_est, use="pairwise.complete.obs")


b <- read.table("~/Blupf90.dat")
colnames(b) <- c("HerdY", "Fenotip", "Indiv", "Sex")
b <- b[,1:3]

PE <- unique(ped[,c("Indiv", "permEnv_est")])
b <- merge(b, PE, by="Indiv")

SOL <- unique(ped[,c("Indiv", "EBV")])
b <- merge(b, SOL, by="Indiv")


HY <- unique(ped[,c("HerdY", "HerdYEff_est")])
b <- merge(b, HY, by="HerdY")

head(b)
nrow(b)
length(unique(b$Indiv))
length(unique(b$HerdY))


b$EBVPE <- b$EBV + b$permEnv_est
b$EBVPEHY <- b$EBV + b$HerdYEff_est

cor(b$Fenotip, b$EBV)
cor(b$Fenotip, b$EBVPE)
cor(b$Fenotip, b$EBVPEHY)


b$HerdY <- as.character(b$HerdY)
b$Year <- sapply(strsplit(b$HerdY,"_"), `[`, 2)

#Cor Pheno - 
phenoEbv <- b %>%
  group_by(Year) %>%
  summarize(COR=cor(Fenotip,EBV))
colnames(phenoEbv)[2] <- "Pheno-EBV"

phenoEbvPE <- b %>%
  group_by(Year) %>%
  summarize(COR=cor(Fenotip,EBVPE))
colnames(phenoEbvPE)[2] <- "Pheno-EbvPE"

phenoEbvPEHY <- b %>%
  group_by(Year) %>%
  summarize(COR=cor(Fenotip,EBVPEHY))
colnames(phenoEbvPEHY)[2] <- "Pheno-EbvPEHY"

phenoCor <- merge(phenoEbv, phenoEbvPE, by="Year")
phenoCor <- merge(phenoCor, phenoEbvPEHY, by="Year")

head(phenoCor)
phenoCorm <- melt(phenoCor, id.vars = "Year")
qplot(data=phenoCorm, x=Year, y=value, colour=variable, group=variable, geom="point")
