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

TGVsAll1 <- data.frame()
for (strategy in c("SU55", "SU51", "SU15")) {
  for (scenario in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
    for (rep in 0:19) {
      base <- TGVsAll$zMean[TGVsAll$scenario=="Class" & TGVsAll$Strategy=="SU55" & TGVsAll$Rep==rep]
      TGVs <- TGVsAll[TGVsAll$Strategy==strategy & TGVsAll$scenario==scenario & TGVsAll$Rep==rep,]
      TGVs$GeneticVarSt <- TGVs$var / TGVs$var[1]
      TGVs$GenicVarSt <- TGVs$AdditGenicVar1 / TGVs$AdditGenicVar1[1]

      
      #add percentage
 #     TGVs$per_zMean <- TGVs$zMean / base

      TGVsAll1 <- rbind(TGVsAll1, TGVs)  
    }
  }
}



write.table(TGVsAll, paste0("TGVsAll_10KRef_", strategy, "_", date, ".csv"), quote=FALSE, row.names=FALSE)

library(ggplot2)
ggplot(data=TGVsAll, aes(x=Generation, y=zMean, colour=scenario)) + geom_line()

pedOCS <- read.table("~/PedOCS.txt", header=TRUE)
pedOCS40 <- ped[ped$Generation %in% 40:60,]
aggregate(pedOCS40$Father ~ pedOCS40$Generation, FUN = function (x) {length(unique(x))})
