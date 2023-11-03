
library(readr)
library(reshape2)

MeanHet <- matrix(ncol=5)

for (rep in c(11, 13, 14:19)) {
#  scenarioHet <- matrix(ncol=20004)
  M <- readr::read_table(paste0("/home/v1jobste/JanaO/SU55/_scenario_", rep, "/SimulatedData/QTNGenotype_Last20Gen.txt"),
                         col_names=FALSE,  col_types = cols(X1 = col_character(), .default = col_integer()))
    for (gen in 1:20) {
      genM <- M[((gen-1)*8640 +1):(gen*8640),]
    
      MeanHet <- rbind(MeanHet, c("SU55", "_scenario_", rep, gen, mean(apply(X = genM, 2,  FUN = function(z) sum(z == 1)) / ncol(genM))))
#      scenarioHet <- rbind(scenarioHet, c("SU55", "_scenario_", rep, gen, apply(X = genM, 2,  FUN = function(z) sum(z == 1)) / nrow(genM) ))
  }
#  write.csv(scenarioHet, paste0("/home/v1jobste/JanaO/Heterozygosity/ScenarioHET_SU55__scenario__Rep", rep, ".csv"), quote=FALSE)
}



write.csv(MeanHet, paste0("/home/v1jobste/JanaO/Heterozygosity/MeanHET1_QTNs_SU55__scenario_.csv"), quote=FALSE)
