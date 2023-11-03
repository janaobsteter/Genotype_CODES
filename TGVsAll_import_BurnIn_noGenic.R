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
library(reshape)
args = commandArgs(trailingOnly=TRUE)
date = args[1]
trait1 = args[2]
trait2 = args[3]
import = args[4]
print(args)
print(args)

homedir <- getwd()

TGVsAll <- data.frame()
for (rep in 0:4) {
  for (scenario in c("Class")) {
    WorkingDir = paste0(homedir, "/BurnIn_TwoPop_", rep, "_", trait1, trait2, "/SimulatedData/")
    UpDir = paste0(homedir, "/BurnIn_TwoPop_", rep, "_", trait1, trait2, "/")
    ped <- read.table(paste0(WorkingDir,'/PedigreeAndGeneticValues.txt'), header=TRUE)
    popsplit <- read.csv(paste0(UpDir, "PopulationSplit.txt"))
    colnames(popsplit) <- c("Group", "Indiv")
    #to standardise onto the generation 40 - which is the generation of comparison
    ped <- ped[ped$Generation %in% 20:40,]
    print("Merging ped and pop")
    ped <- unique(merge(ped, popsplit, by="Indiv", all.x=TRUE))
    #obtain mean and sd of genetic values
    TGV <- ped[,c("Generation", "Group", paste0("gvNormUnres", trait1), paste0("gvNormUnres", trait2))]
    TGVm <- melt(TGV, id.vars = c("Generation", "Group"))
    TGVsum <- summarySE(TGVm, measurevar = "value", groupvars = c("Generation", "variable", "Group"))[,c(1,2,3,5,6)]

    #variance of genetic values
    TGVsum$var <- (TGVsum$sd)^2
    colnames(TGVsum)[1] <- c("Generation")

    TGVs <- data.frame()
    #read in genic variance
    Var <- read.csv(paste0(UpDir,'GenicVariance_import.csv'), header=T)[-1,]

    for (group in c("home", "import")) {
      for (trait in c(trait1, trait2)) {
        TGVGroup <- TGVsum[TGVsum$Group == group & TGVsum$variable == paste0('gvNormUnres', trait),]
        TGVGroup <- TGVGroup[order(TGVGroup$Generation),]
        TGVGroup$SDSt <- TGVGroup$sd / TGVGroup$sd[1]
        #standardise genetic values with genetic standard deviation
        TGVGroup$zMean <- (TGVGroup$value - TGVGroup$value[1]) / TGVGroup$sd[1]
        #TGVs <- merge(TGVs, TGV, by="Generation")
        TGVs <- rbind(TGVs, TGVGroup )     
      }
    }
    
    TGVs$scenario <- scenario
    TGVs$Rep <- rep
    TGVs$Import <- import


    #colnames(TGVs) < c("Generation", paste0("TGV_mean", scenario), paste0("TGV_sd", scenario), paste0("zMean_", scenario), paste0("GenicVar_", scenario), paste0("zMeanGenic_", scenario))
    TGVsAll <- rbind(TGVsAll, TGVs)
  }
}



write.table(TGVsAll, paste0("TGVsAll_import_BurnIn_",  date, ".csv"), quote=FALSE, row.names=FALSE)
