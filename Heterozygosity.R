library(readr)
library(reshape2)

scenarioHet <- data.frame(strategy=NA, scenario=NA, rep=NA, gen=NA, het=NA)
MeanHet <- data.frame()

for (strategy in c("10K_Ref_20Rep", "10K_Ref_1Pb", "10K_Ref_1Year")) {
  for (scenario in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
    for (rep in 0:19) {
      scenarioHet <- data.frame(strategy=NA, scenario=NA, rep=NA, gen=NA, het=NA)
      M <- readr::read_table(paste0("/home/v1jobste/JanaO/", strategy, "/", scenario, rep, "/SimulatedData/AllIndividualsSnpChips/Last20Gen.txt"),
                             col_names=FALSE,  col_types = cols(X1 = col_character(), .default = col_integer()))
        for (gen in 1:20) {
          genM <- M[((gen-1)*8640 +1):(gen*8640),]
        
          meanhet <- mean(apply(X = genM, 2,  FUN = function(z) sum(z == 1))) / ncol(genM)
          scenarioHet <- rbind(scenarioHet, c(strategy, scenario, rep, gen, apply(X = genM, 2,  FUN = function(z) sum(z == 1)) ))
          #MeanHet <- rbind(MeanHet, c(strategy, scenario, rep, gen, meanhet))
          #write.csv(scenarioHet, paste0("/home/v1jobste/JanaO/Heterozygosity/ScenarioHET", strategy, "_", scenario, rep, gen, ".csv"), quote=FALSE)
      }
      write.csv(scenarioHet, paste0("/home/v1jobste/JanaO/Heterozygosity/ScenarioHET", strategy, "_", scenario,"_Rep", rep, ".csv"), quote=FALSE)
    }
    #write.csv(scenarioHet, paste0("/home/v1jobste/JanaO/Heterozygosity/ScenarioHET", strategy, "_", scenario, ".csv"), quote=FALSE)
  }
  #write.csv(scenarioHet, "/home/v1jobste/JanaO/Heterozygosity/ScenarioHET.csv", quote=FALSE)
}

#write.csv(MeanHet, "/home/v1jobste/JanaO/MeanHET.csv", quote=FALSE)



