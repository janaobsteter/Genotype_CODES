
library(readr)
library(reshape2)


for (rep in 10:19) { ##10:19
    MeanHet <- matrix(ncol=5)
    M <- readr::read_table(paste0("/home/v1jobste/JanaO/_strategy_/_scenario_", rep, "/SimulatedData/AllIndividualsSnpChips/Chip1Genotype_Last20Gen.txt"),
                         col_names=FALSE,  col_types = cols(X1 = col_character(), .default = col_integer()))
    for (gen in 1:20) {
        genM <- M[((gen-1)*8640 +1):(gen*8640),]
    
        MeanHet <- rbind(MeanHet, c("_strategy_", "_scenario_", rep, gen, mean(apply(X = genM, 2,  FUN = function(z) sum(z == 1)) / ncol(genM))))
    }
    write.csv(MeanHet, paste0("/home/v1jobste/JanaO/Heterozygosity/MeanHET_Marker__strategy___scenario__Rep", rep, ".csv"), quote=FALSE)
}



#write.csv(MeanHet, paste0("/home/v1jobste/JanaO/Heterozygosity/MeanHET__strategy___scenario_.csv"), quote=FALSE)
