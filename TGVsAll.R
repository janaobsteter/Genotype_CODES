
args = commandArgs(trailingOnly=TRUE)
strategy = args[1]
data = args[2]

################
#Funkcija summarySE
#####################
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
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
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
######################################################################

TGVsAll <- data.frame()
for (rep in 0:20) {
  for (scenario in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
    WorkingDir = paste0("/home/v1jobste/JanaO/10K_Ref_1Year/", scenario, rep, "/SimulatedData/")
    TGVs <- data.frame(Generation=40:60)
    ped <- read.table(paste0(WorkingDir,'/PedigreeAndGeneticValues.txt'), header=T)
    #to standardise onto the generation 40 - which is the generation of comparison
    ped <- ped[ped$Generation %in% 40:60,]
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
    Var <- read.table(paste0(WorkingDir,'TotalGenicAndGeneticVariancesPerGeneration.txt'), header=T)
    #Qtn model 1 is unrestricted 
    Var <- Var[Var$QtnModel==1,c(1,3)]
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
    TGVs$Rep <- rep
    #colnames(TGVs) < c("Generation", paste0("TGV_mean", scenario), paste0("TGV_sd", scenario), paste0("zMean_", scenario), paste0("GenicVar_", scenario), paste0("zMeanGenic_", scenario))
    TGVsAll <- rbind(TGVsAll, TGVs)
  }
}

write.table(TGVsAll,paste0("TGVsAll_", strategy, "_", date, ".csv", quote=FALSE, row.names=FALSE)
