library(readr)
library(reshape2)

scenarioHet <- data.frame()
MeanHet <- data.frame()

for (strategy in c("10K_Ref_20Rep")) {
  for (scenario in c("Class")) {
    for (rep in 0) {
        for (gen in 1) {
          M <- readr::read_table(paste0("/home/v1jobste/JanaO/", strategy, "/", scenario, rep, "/SimulatedData/AllIndividualsSnpChips/Last20Gen.txt"),
                                 col_names=FALSE,  col_types = cols(X1 = col_character(), .default = col_integer()))[((gen-1)*8640 +1):(gen*8640),]
    
          
          meanhet <- mean(apply(X = M, 2,  FUN = function(z) sum(z == 1))) / ncol(M)
          scenarioHet <- rbind(scenarioHet, c(strategy, scenario, rep, gen, apply(X = M, 2,  FUN = function(z) sum(z == 1)) ))
          MeanHet <- rbind(MeanHet, c(strategy, scenario, rep, gen, meanhet))
        }
    }
write.csv(scenarioHet, paste0("/home/v1jobste/JanaO/ScenarioHET", strategy, "_", scenario, ".csv"), quote=FALSE)
  }
write.csv(scenarioHet, "/home/v1jobste/JanaO/ScenarioHET.csv", quote=FALSE)
}

write.csv(MeanHet, "/home/v1jobste/JanaO/MeanHET.csv", quote=FALSE)

