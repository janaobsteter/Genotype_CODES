
library(readr)
library(reshape2)

scenarioHet <- matrix(ncol=20004)


for (rep in 0:19) {
  genM <- readr::read_table(paste0("/home/v1jobste/JanaO/FillInBurnIn", rep, "/SimulatedData/AllIndividualsSnpChips/FirstGen.txt"),
                         col_names=FALSE,  col_types = cols(X1 = col_character(), .default = col_integer()))
    
  scenarioHet <- rbind(scenarioHet, c("_strategy_", "_scenario_", rep, apply(X = genM, 2,  FUN = function(z) sum(z == 1)) / nrow(genM) ))
}



write.csv(scenarioHet, "/home/v1jobste/JanaO/Heterozygosity/HETFirstGen_20BurnIns.csv", quote=FALSE)

