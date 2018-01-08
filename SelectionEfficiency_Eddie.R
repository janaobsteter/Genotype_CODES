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
for (rep in repetition) {
   for (scenario in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
    WorkingDir = paste0("/home/v1jobste/JanaO/", scenario, rep, "/SimulatedData/")
    TGVs <- data.frame(Generation=40:60)
    ped <- read.table(paste0(WorkingDir,'/PedigreeAndGeneticValues.txt'), header=T)
    ped <- ped[ped$Generation %in% 40:60,]
    TGV <- summarySE(ped, measurevar = "gvNormUnres1", groupvars = "Generation")[,c(1,3,4)]
    colnames(TGV)[1] <- c("Generation")
    TGV$zMean <- (TGV$gvNormUnres1 - TGV$gvNormUnres1[1]) / TGV$sd[1]
    TGVs <- merge(TGVs, TGV, by="Generation")
    Var <- read.table(paste0(WorkingDir,'TotalGenicAndGeneticVariancesPerGeneration.txt'), header=T)
    Var <- Var[Var$QtnModel==1,c(1,3)]
    TGVs <- merge(TGVs, Var, by="Generation")
    TGVs$zSdGenic <- (sqrt(TGVs$AdditGenicVar1)) 
    TGVs$SDGenicSt <- TGVs$AdditGenicVar1 / TGVs$AdditGenicVar1[1] 
    TGVs$SDSt <- TGVs$sd / TGVs$sd[1] 
    TGVs$zMeanGenic <- (TGVs$gvNormUnres1 - TGVs$gvNormUnres1[1]) / TGVs$AdditGenicVar1[1]
    #TGVsAll$zSdGenic <- (sqrt(TGVsAll$AdditGenicVar1) - sqrt(TGVsAll$))
    TGVs$scenario <- scenario
    TGVs$Rep <- rep
    #colnames(TGVs) < c("Generation", paste0("TGV_mean", scenario), paste0("TGV_sd", scenario), paste0("zMean_", scenario), paste0("GenicVar_", scenario), paste0("zMeanGenic_", scenario))
    TGVsAll <- rbind(TGVsAll, TGVs)
  }
}

write.table(TGVsAll, paste0(WorkingDir, "/TGVsAll.csv"), quote=FALSE, row.names=FALSE)