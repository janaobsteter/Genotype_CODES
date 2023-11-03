
#!/usr/bin/env Rscript
library(tidyr)
library(dplyr)
#library(ggplot2)
library(pedigree)


args = commandArgs(trailingOnly=TRUE)
rep = args[1]
scenario = args[2]
import = args[3]
trait = args[4]

dir <- paste0(scenario, rep, "_", import, "1", trait)

ped <- read.table(paste0(dir, "/SimulatedData/PedigreeAndGeneticValues.txt"), header=T)
split <- read.csv(paste0(dir, "/PopulationSplit.txt"), header=T)
giH <- read.table(paste0(dir,"/GenIntshome.txt"), header=T)
giI <- read.table(paste0(dir,"/GenIntsimport.txt"), header=T)
het <- read.csv(paste0(dir, "/MeanHet_Marker_import.csv"))[-1,]

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}


# Generation interval
GI <- bind_rows(giH, giI) %>%  mutate(Path = paste0(line, "_", sex)) %>% 
  filter(Gen > 30) %>% 
  group_by(Group) %>% 
  summarize(GI = mean(genInt))

# Heterozygosity
Fhet <- het %>% group_by(Group) %>% 
  filter(Gen > 30) %>% 
  mutate(y1 = log(1 - (1 - HetM))) %>% 
  summarise(dF = 1 - exp(coef(glm(y1 ~ Gen))[2]),
            Ne = 1 / (2 * dF)) %>% 
  full_join(GI) %>% 
  mutate(dF_Gen = dF * GI,
         Ne_Gen = 1 / (2 * dF_Gen),
         scenario = scenario,
         rep = rep,
         import = import,
         trait = trait,
         method = "Heterozygosity")

# Pedigree
colnames(split) <- c("Group", "Indiv")
ped <- full_join(ped, split, by = "Indiv")
ped$Group <- as.factor(ped$Group)

ped$Fped <- calcInbreeding(ped[, c("Indiv", "Father", "Mother")])
Fped <- ped %>%  group_by(Generation, Group) %>% summarise(Fmean = mean(Fped)) %>% 
  filter(Generation > 30) %>% 
  mutate(y1 = log(1 - Fmean)) %>% 
  group_by(Group) %>% 
  summarise(dF = 1 - exp(coef(glm(y1 ~ Generation))[2]),
            Ne = 1 / (2 * dF)) %>% 
  full_join(GI) %>% 
  mutate(dF_Gen = dF * GI,
         Ne_Gen = 1 / (2 * dF_Gen),
         scenario = scenario,
         rep = rep,
         import = import,
         trait = trait,
         method = "Pedigree") 

write.csv(bind_rows(Fhet, Fped), paste0("Inbreeding_", scenario, rep, "_", import, "1", trait, ".csv"), quote=F, row.names=F)
