
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

library(tidyr)
homedir <- getwd()

MST <- data.frame()
MSTCat <- data.frame()

for (strategy in c("SU55", "SU51", "SU15")) {
  MSTStrategy <- data.frame()  	
  for (scenario in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
    for (rep in 0:19) {
      setwd(paste0(homedir, "/", strategy, "/", scenario, rep, "/"))
      print(paste0(homedir, "/", strategy, "/", scenario, rep, "/"))
      if (file.info("GenPed_EBV.txt")$size > 0) {
        ped <- read.csv("GenPed_EBV.txt")
        print("Reading GenPed")
        C <- gather(read.csv("Categories_gen60DF.csv"))
        C <- C[!(is.na(C$value)),]
        colnames(C) <- c("cat", "Indiv")
        ped <- merge(ped, C, by="Indiv", all.x=TRUE)
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
      ped1$PA <- (ped1$gvF +  ped1$gvM) / 2
      
      mstAa <- summarySE(data=ped1, measurevar = "MST", groupvars = "Generation")[,c(1,3,4)]
      colnames(mstAa) <- c("Generation", "MSTmean", "MSTsd")
      mstAb <- summarySE(data=ped1, measurevar = "PA", groupvars = "Generation")[,c(1,3,4)]
      colnames(mstAb) <-	c("Generation",	"PAmean", "PAsd")
      mstA <- merge(mstAa, mstAb, by="Generation")


      pedPB <- ped1[ped1$cat %in% c("pb", "gpb"),]
      PBmst <- summarySE(data=pedPB, measurevar = "MST", groupvars = "Generation")[,c(1,3,4)]
      colnames(PBmst) <- c("Generation", "MSTmean_PB", "MSTsd_PB")
      PBpa <- summarySE(data=pedPB, measurevar = "PA", groupvars = "Generation")[,c(1,3,4)]
      colnames(PBpa) <-        c("Generation", "PAmean_PB", "PAsd_PB")
      mstA <- merge(mstA, PBmst, by="Generation")
      mstA <- merge(mstA, PBpa, by="Generation")

      mstA$Scenario <- scenario
      mstA$Rep <- rep
      mstA$Strategy <- strategy


      
      mstCata <- summarySE(data=ped1, measurevar = "MST", groupvars = "cat")[,c(1,3,4)]
      colnames(mstCata) <-  c("Cat",	"MSTmean", "MSTsd")
      mstCatb <- summarySE(data=ped1, measurevar = "PA", groupvars = "cat")[,c(1,3,4)]
      colnames(mstCatb) <- c("Cat", "PAmean", "PAsd")
      mstCat <- merge(mstCata, mstCatb, by="Cat")

      mstCat$Scenario <- scenario
      mstCat$Rep <- rep
      
      MST <- rbind(MST, mstA)
      MSTCat <- rbind(MSTCat, mstCat)
      MSTStrategy <- rbind(MSTStrategy, mstA)
      
      setwd(homedir)
    }
  }
    write.csv(MSTStrategy,paste0( "MST_", strategy, ".csv"), quote=FALSE, row.names=FALSE)
}
write.csv(MST, "MST.csv", quote=FALSE, row.names=FALSE)
write.csv(MSTCat, "MSTCat.csv", quote=FALSE, row.names=FALSE)

