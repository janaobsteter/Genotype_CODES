################
#Funkcija summarySE
#####################
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=TRUE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  colnames(datac)[colnames(datac) == "mean"] <- measurevar
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
######################################################################
args = commandArgs(trailingOnly=TRUE)
strategy = args[1]
date = args[2]
print(args)
print(strategy)
print(args)

setwd("~/Documents/PhD/Projects/inProgress/GenomicStrategies_Import/")



TGVs <- data.frame(Generation=0:60)
#ped <- read.table(paste0('~/PedOCS.txt'), header=TRUE)
pedT <- read.table('~/Import/PedigreeAndGeneticValues.txt', header=TRUE)
VarT <- read.csv('~/Import/GenicVariance_import.csv', header=T)[-1,]
colnames(VarT)[colnames(Var) == "Gen"] <- "Generation"
#Var <- Var[Var$QtnModel==1,c(1,3)]
#split populations
popSplit <- read.csv("~/Import/PopulationSplit.txt")
pedH <- ped[ped$Indiv %in% popSplit$ID[popSplit$Group=="home"],]
pedI <- ped[ped$Indiv %in% popSplit$ID[popSplit$Group=="import"],]
length(intersect(pedH$Indiv, pedI$Indiv))
#to standardise onto the generation 40 - which is the generation of comparison
TGVsAll <- data.frame()
for (group in c("home", "import")) {
  TGVs <- data.frame(Generation=40:60)
  ped <- pedT[pedT$Indiv %in% popSplit$ID[popSplit$Group==group],]
  ped <- ped[ped$Generation %in% 40:60,]
  #obtain mean and sd of genetic values
  TGV <- summarySE(data = ped, measurevar = "gvNormUnres1", groupvars = "Generation")[,c(1,3,4)]
  #variance of genetic values
  TGV$var <- (TGV$sd)^2
  colnames(TGV)[1] <- c("Generation")
  #standardise genetic standard devistion
  TGV$SDSt <- TGV$sd / TGV$sd[1]
  #standardise genetic values with genetic standard deviation
  TGV$zMean <- (TGV$gvNormUnres1 - TGV$gvNormUnres1[1]) / TGV$sd[1]
  TGVs <- merge(TGVs, TGV, by="Generation")
  #read in genic variance
  Var <- VarT[VarT$Group==group,]
  
  #Qtn model 1 is unrestricted 

  TGVs <- merge(TGVs, Var, by="Generation")
  #obtain genic standard deviation
  TGVs$SDGenic <- (sqrt(TGVs$GenVar))
  #standarise genic standard devistion
  TGVs$SDGenicSt <- TGVs$SDGenic / TGVs$SDGenic[1]
  #standardise genetic values with genic standard devistion
  TGVs$zMeanGenic <- (TGVs$gvNormUnres1 - TGVs$gvNormUnres1[1]) / TGVs$SDGenic[1]
  #reciprocated genic standard deviation
  TGVs$SDGenicStNeg <- 1 - (TGVs$SDGenic / TGVs$SDGenic[1])
  #genic variance standardised onto genetic variance
  koef <- TGVs$var[1] / TGVs$GenVar[1]
  TGVs$Genic_Genetic_VAR <- TGVs$GenVar * koef
  TGVs$Genic_Genetic_SD <- sqrt(TGVs$Genic_Genetic_VAR)
  #standardise genic_genetic standard deviation
  TGVs$Genic_Genetic_SDSt <- TGVs$Genic_Genetic_SD / TGVs$Genic_Genetic_SD[1]
  #standarise genetic values with genic_genetic standard deviation
  TGVs$zMeanGenic_Genetic <- (TGVs$gvNormUnres1 - TGVs$gvNormUnres1[1]) / TGVs$Genic_Genetic_SD[1]
  #TGVsAll$zSdGenic <- (sqrt(TGVsAll$AdditGenicVar1) - sqrt(TGVsAll$))
  TGVs$scenario <- scenario
  TGVs$Rep <- 0
  TGVs$Group = group
  #colnames(TGVs) < c("Generation", paste0("TGV_mean", scenario), paste0("TGV_sd", scenario), paste0("zMean_", scenario), paste0("GenicVar_", scenario), paste0("zMeanGenic_", scenario))
  TGVsAll <- rbind(TGVsAll, TGVs)
}

#preveri očete
pedC <- read.table("~/Import/PedigreeAndGeneticValues_cat.txt", header=TRUE)
table(popSplit$Group[popSplit$ID %in% pedC$Father[pedC$cat=="potomciNP"]])
table(popSplit$Group[popSplit$ID %in% pedC$Father[pedC$cat=="nr"]])
mean(pedC$gvNormUnres1[pedC$Indiv %in% pedC$Father[pedC$cat=="potomciNP"]])
mean(pedC$gvNormUnres1[pedC$Indiv %in% pedC$Father[pedC$cat=="nr"]])

colnames(popSplit)[2] <- "Indiv"
pedG <- merge(pedC, popSplit, by="Indiv")

mean(pedG$gvNormUnres1[pedG$cat=="pb" & pedG$Group=="import"])
mean(pedG$gvNormUnres1[pedG$cat=="pb" & pedG$Group=="home"])
ped20 <- pedG[pedG$Generation %in% 40:60,]
mean(ped20$gvNormUnres1[ped20$cat=="pb" & ped20$Group=="import"])
mean(ped20$gvNormUnres1[ped20$cat=="pb" & ped20$Group=="home"])
mean(pedG$gvNormUnres1[pedG$cat=="pripust1" & pedG$Group=="import"])
mean(pedG$gvNormUnres1[pedG$cat=="pripust1" & pedG$Group=="home"])
mean(pedG$gvNormUnres1[pedG$cat=="mladi" & pedG$Group=="import"])
mean(pedG$gvNormUnres1[pedG$cat=="mladi" & pedG$Group=="home"])
mean(pedG$gvNormUnres1[pedG$cat=="k" & pedG$Group=="import"])
mean(pedG$gvNormUnres1[pedG$cat=="k" & pedG$Group=="home"])

pedPB <- pedG[pedG$cat=="pb",]
PBa <- aggregate(pedPB$gvNormRestr1 ~ pedPB$Generation + pedPB$Group, FUN="mean")
colnames(PBa) <- c("Generation", "Group", "TGV")
ggplot(data=PBa, aes(x=Generation, y=TGV, group=Group, colour=Group)) + geom_line() 


colnames(popSplit) <- c("Group", "Indiv")
pedG <- merge(pedC, popSplit, by="Indiv")
mean(pedG$gvNormUnres1[pedG$cat=="pb" & pedG$Group=="import"])
mean(pedG$gvNormUnres1[pedG$cat=="pb" & pedG$Group=="home"])
table(pedG$Group[pedG$Indiv %in% pedG$Father[pedG$cat=="potomciNP"]])
table(pedG$Group[pedG$Indiv %in% pedG$Father[pedG$cat=="nr" & pedG$Group=="home"]])
table(pedG$cat[pedG$Indiv %in% pedG$Father[pedG$cat=="nr"]], pedG$Group[pedG$Indiv %in% pedG$Father[pedG$cat=="nr"]])

write.table(TGVsAll, paste0("TGVsAll_10KRef_", strategy, "_", date, ".csv"), quote=FALSE, row.names=FALSE)

library(ggplot2)
ggplot(data=TGVsAll, aes(x=Generation, y=zMean, colour=Group)) + geom_line() + 
  scale_y_continuous(breaks = round(seq(min(TGVsAll$zMean), max(TGVsAll$zMean), by = 0.5),1)) +
  theme(text = element_text(size=20))

#genicVariance
ggplot(data=TGVsAll, aes(x=Generation, y=GenVar, colour=Group)) + geom_line() + 
  scale_y_continuous() +
  theme(text = element_text(size=20))

# pedOCS <- read.table("~/PedOCS.txt", header=TRUE)
pedOCS40 <- ped[ped$Generation %in% 40:60,]
aggregate(pedOCS40$Father ~ pedOCS40$Generation, FUN = function (x) {length(unique(x))})

#all in one
TGVsAll <- data.frame()
for (group in c("home", "import")) {
  TGVs <- data.frame(Generation=0:60)
  ped <- read.table('~/Import/PedigreeAndGeneticValues.txt', header=TRUE)
  split <- read.csv("~/Import/PopulationSplit.txt")
  colnames(split) <- c("Group", "Indiv")
  ped <- merge(ped, split, by="Indiv")
  #ped <- ped[ped$Generation %in% 0:40,]
  #obtain mean and sd of genetic values
  TGV <- summarySE(ped, measurevar = "gvNormUnres1", groupvars = "Generation")[,c(1,3,4)]
  #variance of genetic values
  TGV$var <- (TGV$sd)^2
  colnames(TGV)[1] <- c("Generation")
  #standardise genetic standard devistion
  TGV$SDSt <- TGV$sd / TGV$sd[1]
  #standardise genetic values with genetic standard deviation
  TGV$zMean <- (TGV$gvNormUnres1 - TGV$gvNormUnres1[1]) / TGV$sd[1]
  TGVs <- merge(TGVs, TGV, by="Generation")
  #read in genic variance
  #Var <- read.table(paste0('~/VarOCS.txt'), header=T)
  
  #Qtn model 1 is unrestricted 
  
  TGVs <- merge(TGVs, Var, by="Generation")
  #obtain genic standard deviation
  TGVs$SDGenic <- (sqrt(TGVs$AdditGenicVar1))
  #standarise genic standard devistion
  TGVs$SDGenicSt <- TGVs$SDGenic / TGVs$SDGenic[1]
  #standardise genetic values with genic standard devistion
  TGVs$zMeanGenic <- (TGVs$gvNormUnres1 - TGVs$gvNormUnres1[1]) / TGVs$SDGenic[1]
  #reciprocated genic standard deviation
  TGVs$SDGenicStNeg <- 1 - (TGVs$SDGenic / TGVs$SDGenic[1])
  #genic variance standardised onto genetic variance
  koef <- TGVs$var[1] / TGVs$AdditGenicVar1[1]
  TGVs$Genic_Genetic_VAR <- TGVs$AdditGenicVar1 * koef
  TGVs$Genic_Genetic_SD <- sqrt(TGVs$Genic_Genetic_VAR)
  #standardise genic_genetic standard deviation
  TGVs$Genic_Genetic_SDSt <- TGVs$Genic_Genetic_SD / TGVs$Genic_Genetic_SD[1]
  #standarise genetic values with genic_genetic standard deviation
  TGVs$zMeanGenic_Genetic <- (TGVs$gvNormUnres1 - TGVs$gvNormUnres1[1]) / TGVs$Genic_Genetic_SD[1]
  #TGVsAll$zSdGenic <- (sqrt(TGVsAll$AdditGenicVar1) - sqrt(TGVsAll$))
  TGVs$scenario <- scenario
  TGVs$Rep <- 0
  TGVs$Group = group
  #colnames(TGVs) < c("Generation", paste0("TGV_mean", scenario), paste0("TGV_sd", scenario), paste0("zMean_", scenario), paste0("GenicVar_", scenario), paste0("zMeanGenic_", scenario))
  TGVsAll <- rbind(TGVsAll, TGVs)
}



PreseTheme <- theme_bw(base_size=18, base_family="sans")  + theme(legend.position="top")
PreseSize <- 16
p = ggplot(data=pedG, aes(y=Generation)) +
  geom_density_ridges(aes(x=gvNormUnres1, fill=paste(Generation, Group), colour=Group), rel_min_height = 0.01, alpha=.8) +
  scale_fill_cyclical(values = c(ColB, ColV, ColB2, ColV2) # c("#ff8080", "#8080ff", "#ff0000", "#0000ff"),
                      )  + scale_colour_manual(values=c("steelblue3", "deeppink3")) +
  labs(x="Genetic value",
       y="Generation") +
  theme_ridges(grid=FALSE)
p + PreseTheme

PreseTheme <- theme_bw(base_size=18, base_family="sans")  + theme(legend.position="top")
PreseSize <- 16
p = ggplot(data=pedG[pedG$cat=="pb",], aes(y=Generation)) +
  geom_density_ridges(aes(x=gvNormUnres1, fill=paste(Generation, Group), colour=Group), rel_min_height = 0.01, alpha=.8) +
  scale_fill_cyclical(values = c(ColB, ColV, ColB2, ColV2) # c("#ff8080", "#8080ff", "#ff0000", "#0000ff"),
                      ) + scale_colour_manual(values=c("steelblue3", "deeppink3")) +
  labs(x="Genetic value",
       y="Generation") +
  theme_ridges(grid=FALSE)
p + PreseTheme
