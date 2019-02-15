
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


ACC <- data.frame(Degree=NA, Rep=NA, Gen=NA, Cor=NA)
for (scenario in c(15, 30, 45, 60, 75)) {
  for (rep in 0:19) {
    setwd(paste0(homedir, "/Gen", rep , "_", scenario, "OCS/"))
    ped <- read.csv("GenPed_EBV.txt")[,1:5]
    
    for (gen in 41:60) {
      sol <- read.table(paste0("renumbered_Solutions_", gen - 1))[,c(2,3)]
      colnames(sol) <- c("Indiv", "EBV")
      fathers <- ped$Father[ped$Generation==gen]
      fatherGen <- ped[ped$Indiv %in% fathers,]
      fatherSol <- merge(fatherGen, sol, by="Indiv", all.x=TRUE)    
      ACC <- rbind(ACC, c(scenario, rep, gen, cor(fatherSol$gvNormUnres1, fatherSol$EBV)))
      }
    }
}
write.csv(ACC, "ACC_OCS.csv", quote=FALSE, row.names=FALSE)

##########################
#MST
###############################



library(tidyr)
homedir <- getwd()


for (strategy in c("SU55", "SU51", "SU15")) {
  MST <- data.frame()
  for (scenario in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
    for (rep in 0:19) {
      setwd(paste0(homedir, "/", strategy, "/", scenario, rep, "/"))
      print(paste0(homedir, "/", strategy, "/", scenario, rep, "/"))
      if (file.info("GenPed_EBV.txt")$size > 0) {
        ped <- read.csv("GenPed_EBV.txt")
        print("Reading GenPed")
        C <- gather(read.csv("Categories_gen60DF.csv"))
        C <- C[!(is.na(C$value)),]
        colnames(C) <- c("Cat", "Indiv")
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
      
      ped1 <- ped[ped$Generation %in% 40:60,]
      nrow(ped1)
      ped1 <- merge(ped1, pedFather, by="Father", all.x=TRUE)
      nrow(ped1)
      ped1 <- merge(ped1, pedMother, by="Mother", all.x=TRUE)
      nrow(ped1)
      
      ped1$MST <- ped1$gvNormUnres1 -  ((ped1$gvF +  ped1$gvM) / 2)
      
      mstA <- summarySE(data=ped1, measurevar = "MST", groupvars = "Generation")[,c(1,3,4)]
      mstA$Scenario <- scenario
      mstA$Rep <- rep
      
      mstCat <- summarySE(data=ped1, measurevar = "MST", groupvars = c("Cat"))#[,c(1,3,4)]
      mstCat$Scenario <- scenario
      mstCat$Rep <- rep
      
      MST <- rbind(MST, mstA)
      MSTCat <- rbind(MSTCat, mstCat)
      
      setwd(homedir)
    }
    write.csv(MST,paste0( "MST_", scenario, ".csv"), quote=FALSE, row.names=FALSE)
    write.csv(MSTCat ,paste0( "MSTCat_", scenario, ".csv"), quote=FALSE, row.names=FALSE)
  }
write.csv(MST, "MST_OCS.csv", quote=FALSE, row.names=FALSE)
write.csv(MSTCat, "MSTCat.csv", quote=FALSE, row.names=FALSE)

########################################################
#OCS
########################################################
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
    
    ped1 <- ped[ped$Generation %in% 40:60,]
    nrow(ped1)
    ped1 <- merge(ped1, pedFather, by="Father", all.x=TRUE)
    nrow(ped1)
    ped1 <- merge(ped1, pedMother, by="Mother", all.x=TRUE)
    nrow(ped1)
    
    ped1$MST <- ped1$gvNormUnres1 -  ((ped1$gvF +  ped1$gvM) / 2)
    
    mstA <- summarySE(data=ped1, measurevar = "MST", groupvars = "Generation")[,c(1,3,4)]
    mstA$Scenario <- scenario
    mstA$Rep <- rep
        
    #only fathers
    ped2 <- ped[ped$Indiv %in% fathers,]
    nrow(ped2)
    ped2 <- merge(ped2, pedFather, by="Father", all.x=TRUE)
    nrow(ped2)
    ped2 <- merge(ped2, pedMother, by="Mother", all.x=TRUE)
    nrow(ped2)
    
    ped2$MST <- ped2$gvNormUnres1 -  ((ped2$gvF +  ped2$gvM) / 2)
    
    mstCat <- summarySE(data=ped2, measurevar = "MST", groupvars = c("Generation"))#[,c(1,3,4)]
    mstCat$Scenario <- scenario
    mstCat$Rep <- rep
    
    
    MST <- rbind(MST, mstA)
    MSTCat <- rbind(MSTCat, mstCat)
    
    setwd(homedir)
    }
  }
  write.csv(MST,paste0( "MST_", scenario, "OCS.csv"), quote=FALSE, row.names=FALSE)
  write.csv(MSTCat ,paste0( "MSTCat_", scenario, "OCS.csv"), quote=FALSE, row.names=FALSE)
}
write.csv(MST, "MST_OCS.csv", quote=FALSE, row.names=FALSE)
write.csv(MSTCat, "MSTCat_OCS.csv", quote=FALSE, row.names=FALSE)



