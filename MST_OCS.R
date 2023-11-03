
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

library(tidyr)
MST <- data.frame()
MSTCat <- data.frame()

for (scenario in c(15, 30, 45, 60, 75)) {
  for (rep in 0:19) {
    setwd(paste0(homedir, "/Gen", rep , "_", scenario, "OCS/"))
    print(paste0(homedir, "/Gen", rep , "_", scenario, "OCS/"))
    if (file.info("GenPed_EBV.txt")$size > 0) {
      ped <- read.csv("GenPed_EBV.txt")
      print("Reading GenPed")
      #C <- gather(read.csv("Categories_gen60DF.csv"))
      #C <- C[!(is.na(C$value)),]
      #colnames(C) <- c("cat", "Indiv")
      #ped <- merge(ped, C, by="Indiv", all.x=TRUE)
    }  else {
      ped <- read.table("./SimulatedData/PedigreeAndGeneticValues_cat.txt", header=TRUE)
      print("Reading Pedigree")
    }
    ped$MST <- NA
    
    fathers <-  unique(ped$Father[ped$Generation %in% 40:60])
    mothers <-  unique(ped$Mother[ped$Generation %in% 40:60])
    pedFather <- ped[ped$Indiv %in% fathers,c("Indiv", "gvNormUnres1")]
    pedMother <- ped[ped$Indiv %in% mothers,c("Indiv", "gvNormUnres1")]
    colnames(pedFather) <- c("Father", "gvF")
    colnames(pedMother) <- c("Mother", "gvM")
    
    ped1 <- merge(ped, pedFather, by="Father", all.x=TRUE)
    nrow(ped1)
    ped1 <- merge(ped1, pedMother, by="Mother", all.x=TRUE)
    nrow(ped1)
    ped1 <- ped1[ped1$Generation %in% 40:60,]
    nrow(ped1)
    
    ped1$MST <- ped1$gvNormUnres1 -  ((ped1$gvF +  ped1$gvM) / 2)
    
    mstA <- summarySE(data=ped1, measurevar = "MST", groupvars = "Generation")[,c(1,3,4)]
    mstA$Scenario <- scenario
    mstA$Rep <- rep
        
    #only fathers
    ped2 <- ped[ped$Generation %in% 40:60 & ped$Indiv %in% fathers,]
    nrow(ped2)
    ped2 <- merge(ped2, pedFather, by="Father", all.x=TRUE)
    nrow(ped2)
    ped2 <- merge(ped2, pedMother, by="Mother", all.x=TRUE)
    nrow(ped2)
    
    ped2$MST <- ped2$gvNormUnres1 -  ((ped2$gvF +  ped2$gvM) / 2)
    
    mstCatm <- aggregate(ped2$MST ~ ped2$Generation, FUN="mean")
    colnames(mstCatm) <- c("Generation", "mean")
    mstCatsd <- aggregate(ped2$MST ~ ped2$Generation, FUN="sd") 
    colnames(mstCatsd) <- c("Generation", "sd")
    mstCat <- merge(mstCatm, mstCatsd, by="Generation") 

    mstCat$Scenario <- scenario
    mstCat$Rep <- rep
    
    
    MST <- rbind(MST, mstA)
    MSTCat <- rbind(MSTCat, mstCat)
    
    setwd(homedir)
  }
}
write.csv(MST, "MST_OCS_40.csv", quote=FALSE, row.names=FALSE)
write.csv(MSTCat, "MSTCat_OCS_40.csv", quote=FALSE, row.names=FALSE)
