library(tidyr)

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

homedir <- getwd()
print(paste("Home dir is", homedir, sep=" "))

print("Creating empty data frame")
ACC <- data.frame(Degree=NA, Rep=NA, Gen=NA, Cor=NA)
for (scenario in c(50, 55)) {
  for (rep in c(0:9)) {
    setwd(paste0(homedir, "/Gen", rep , "_", scenario, "OCS/"))
    print(paste0(homedir, "/Gen", rep , "_", scenario, "OCS/"))
    ped <- read.csv("GenPed_EBV.txt")[,1:5]
    for (gen in 41:max(ped$Generation)) {
        print(gen)
        
        #preberi kategorije
        cat <- read.csv(paste0("Categories_gen", gen, "DF.csv"))
        C <- gather(cat)
        colnames(C) <- c("Cat", "Indiv")
        C <- C[!(is.na(C$Indiv)),]
        
        #preberi solutions
        sol <- read.table(paste0("renumbered_Solutions_", gen))[,c(2,3)]
        colnames(sol) <- c("Indiv", "EBV")
        
        if (nrow(C) == nrow(sol)) {
          #združi
          P <- merge(ped, C, by="Indiv", all.y=TRUE)
          P <- merge(P, sol, by="Indiv", all.x=TRUE)
          
          #samo zdnjih 6 let - toliko let uporabljaš / consider bike
          P <- P[P$Generation %in% ((max(P$Generation) - 6):max(P$Generation)),]
          
          #tukaj pridobi kandidate za očete
          kand <- P[P$Cat %in% c("vhlevljeni", "mladi", "cak", "genTest", "gpb"),]
          
          ACC <- rbind(ACC, c(scenario, rep, gen, cor(kand$gvNormUnres1, kand$EBV)))
        }
        
      }
    }
  setwd(homedir)
  print(homedir)
  print("Writing csv")
  write.csv(ACC, paste0("./Accuracies/ACC_OCS_", scenario, ".csv"), quote=FALSE, row.names=FALSE)
}
setwd(homedir)
write.csv(ACC, "ACC_OCS_5055.csv", quote=FALSE, row.names=FALSE)
