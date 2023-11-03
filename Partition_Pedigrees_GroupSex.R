#########################################################################
### A script to partition the breeding values of the import scenarios
#########################################################################

args = commandArgs(trailingOnly=TRUE)
date = args[1]

## --- Libraries --- ##
library(reshape)
library(AlphaPart)

## --- Functions --- ##
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
#######################################################################

# Initialise a dataframe to hold the partitions information
partitions <- data.frame()
homedir = getwd()
trait1 = 1

# Read in the pedigree for each of the replicates (rep), for each of the correlated trait (trait 2 has r = 0.9 and trait 3 r = 0.8)
# ClassGen scenarios implements genomic selection with a 10-year delay in the home population, GenGen implements genomic selection in both populations simulteneously
# Import means import percentage importforcows_importforbulldams
# The directories are set for the Eddie server
for (rep in c(11, 12)) {
  for (trait2 in c(2,3)) {
    for (scenarioHome in c("Gen")) {
      for (import in c("0_100", "100_100", "50_50", "25_25", "10_10", "0_0")) {
        scenario = paste0(scenarioHome, "Gen")
        # Set direcotires
        WorkingDir = paste0(homedir, "/10K/SU55_import/", scenario, rep, "_", import, trait1, trait2, "/SimulatedData/")
        UpDir = paste0(homedir, "/10K/SU55_import/", scenario, rep, "_", import, trait1, trait2, "/")
        # If the pedigree file exists
        if (file.exists(paste0(WorkingDir,'/PedigreeAndGeneticValues.txt'))) {
        # Read in pedigree
	  if (file.exists(paste0(WorkingDir,'/PedigreeAndGeneticValues.txt'))) {
          ped <- read.table(paste0(WorkingDir,'/PedigreeAndGeneticValues.txt'), header=TRUE)
          # Read in the file holding the information about the population of the individuals
          popsplit <- read.csv(paste0(UpDir, "PopulationSplit.txt"))
          colnames(popsplit) <- c("Group", "Indiv")
          sex <- read.table(paste0(WorkingDir, "Gender.txt"), header=T)
	  ped <- merge(ped, sex[,c("Indiv", "Gender")], by="Indiv")
          # Keep only the pedigree from generation 20 onward (1 - 20 is fill in: random breeding)
          ped <- pedSetBase(ped, keep=ped$Generation > 19, unknown=0, colId = "Indiv", colFid = "Father", colMid = "Mother")
          # Merge ped in population information
          ped <- unique(merge(ped, popsplit, by="Indiv", all.x=TRUE))
          # Standardize the trait 1
          ped$z1 <- ped$gvNormUnres1 - mean(ped$gvNormUnres1[ped$Generation == 20 & ped$Group == "home"])
	  ped$GroupSex <- paste0(ped$Group, ped$Gender)
	  ped$GenerationSex <- paste0(ped$Generation, ped$Gender)
          # Partition the breeding values
          part = AlphaPart(x = ped, colPath= "Group", colBV = "z1", colId = "Indiv", colFid = "Father", colMid = "Mother")
          # Summarise the breeding values
          sumPartHome <- summary(part, by="GenerationSex", subset = part$z1$Group == "home")
          sumPartHome <- sumPartHome$z1
          sumPartHome$Population <- "home"
          sumPartHome$Rep <- rep
          sumPartHome$Import <- import
          sumPartHome$Scenario <- scenario
          sumPartHome$Cor <- ifelse(trait2 == 2, 0.9, 0.8)
          sumPartImport <- summary(part, by="GenerationSex", subset = part$z1$Group == "import")
          sumPartImport <- sumPartImport$z1
          sumPartImport$Population <- "import"
          sumPartImport$Rep <- rep
          sumPartImport$Import <- import
          sumPartImport$Scenario <- scenario
          sumPartImport$Cor <- ifelse(trait2 == 2, 0.9, 0.8)
          # Save the partitions
          partitions <- rbind(partitions, sumPartHome)
          partitions <- rbind(partitions, sumPartImport)
	}
        }
      }
    }
  }
}

write.table(partitions, paste0("Partition_Import_", date, ".csv"), quote=F, row.names=F)
        
        
        
        
        
        
