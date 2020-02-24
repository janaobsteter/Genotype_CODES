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
args = commandArgs(trailingOnly=TRUE)
strategy = args[1]
date = args[2]
print(args)
print(strategy)
print(args)

TGVsAll <- data.frame()
for (rep in 0:1) {
  for (import in c(75, 100)) {
    for (scenario in c("Class")) {
      WorkingDir = paste0("/home/v1jobste/JanaO/10K/",strategy, "_import/", scenario, rep, "_", import, "/SimulatedData/")
      UpDir = paste0("/home/v1jobste/JanaO/10K/",strategy, "_import/", scenario, rep, "_", import, "/")
      ped <- read.table(paste0(WorkingDir,'/PedigreeAndGeneticValues.txt'), header=TRUE)
      popsplit <- read.csv(paste0(UpDir, "PopulationSplit.txt"))
      colnames(popsplit) <- c("Group", "Indiv")
      #to standardise onto the generation 40 - which is the generation of comparison
      ped <- ped[ped$Generation %in% 40:60,]
      ped <- unique(merge(ped, popsplit, by="Indiv", all.x=TRUE))
      ped$import <- 100
      #obtain mean and sd of genetic values
      TGV <- summarySE(ped, measurevar = "gvNormUnres1", groupvars = c("Generation", "Group"))[,c(1,2,4,5)]
      #variance of genetic values
      TGV$var <- (TGV$sd)^2
      colnames(TGV)[1] <- c("Generation")

      TGVs <- data.frame()
      #read in genic variance
      Var <- read.csv(paste0(UpDir,'GenicVariance_import.csv'), header=T)[-1,]
      for (group in c("home", "import")) {
        TGVGroup <- TGV[TGV$Group == group,]
        TGVGroup <- TGVGroup[order(TGVGroup$Generation),]
        TGVGroup$SDSt <- TGVGroup$sd / TGVGroup$sd[1]
        #standardise genetic values with genetic standard deviation
        TGVGroup$zMean <- (TGVGroup$gvNormUnres1 - TGVGroup$gvNormUnres1[1]) / TGVGroup$sd[1]
        #TGVs <- merge(TGVs, TGV, by="Generation")
        
        #variation
        VarGroup <- Var[Var$Group == group,]
        colnames(VarGroup)[5] <- "Generation"
        TGVGroup <- merge(TGVGroup, VarGroup, by=c("Generation", "Group"))
        
        #obtain genic standard deviation
        TGVGroup$SDGenic <- (sqrt(TGVGroup$GenVar))
        
        #standarise genic standard devistion
        TGVGroup$SDGenicSt <- TGVGroup$SDGenic / TGVGroup$SDGenic[1]
        #standardise genetic values with genic standard devistion
        TGVGroup$zMeanGenic <- (TGVGroup$gvNormUnres1 - TGVGroup$gvNormUnres1[1]) / TGVGroup$SDGenic[1]
        #reciprocated genic standard deviation
        TGVGroup$SDGenicStNeg <- 1 - (TGVGroup$SDGenic / TGVGroup$SDGenic[1])
        #genic variance standardised onto genetic variance
        koef <- TGVGroup$var[1] / TGVGroup$GenVar[1]
        TGVGroup$Genic_Genetic_VAR <- TGVGroup$GenVar * koef
        TGVGroup$Genic_Genetic_SD <- sqrt(TGVGroup$Genic_Genetic_VAR)
        #standardise genic_genetic standard deviation
        TGVGroup$Genic_Genetic_SDSt <- TGVGroup$Genic_Genetic_SD / TGVGroup$Genic_Genetic_SD[1]
        #standarise genetic values with genic_genetic standard deviation
        TGVGroup$zMeanGenic_Genetic <- (TGVGroup$gvNormUnres1 - TGVGroup$gvNormUnres1[1]) / TGVGroup$Genic_Genetic_SD[1]
        #TGVsAll$zSdGenic <- (sqrt(TGVsAll$AdditGenicVar1) - sqrt(TGVsAll$))
	TGVs <- rbind(TGVs, TGVGroup )        
      }
      TGVs$scenario <- scenario
      TGVs$Rep <- rep
      TGVs$Strategy <- strategy
      TGVs$Import <- import


      #colnames(TGVs) < c("Generation", paste0("TGV_mean", scenario), paste0("TGV_sd", scenario), paste0("zMean_", scenario), paste0("GenicVar_", scenario), paste0("zMeanGenic_", scenario))
      TGVsAll <- rbind(TGVsAll, TGVs)
    }
  }
}


write.table(TGVsAll, paste0("TGVsAll_import_", strategy, "_", date, ".csv"), quote=FALSE, row.names=FALSE)
