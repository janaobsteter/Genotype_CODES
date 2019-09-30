stFactor = data.frame(Re[ = NA, h2 = NA, Program = NA, Trait = NA, BV = NA, Variable = NA, Value = NA])


for (rep in c(0)) {
  for (h2 in c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.99)) {
    #program 1
    PedEval1 <- read.table(paste0("Rep", rep, "/H_", h2, "/PN1/PedEval1_", h2, ".csv"), header=TRUE)
    #program 2
    PedEval2 <- read.table(paste0("Rep", rep, "/H_", h2, "/PN2/PedEval2_", h2, ".csv"), header=TRUE)
    
    # ---- Partitioning the trend PN1 ----
    stFactor <- rbind(stFactor, c(rep, h2, "Program1", "Trait1", "Ebv", "Sd", sd(PedEval1$EbvT1[PedEval1$Generation == 0])
    stFactor <- rbind(stFactor, c(rep, h2, "Program2", "Trait1", "Ebv", "Sd",  sd(PedEval2$EbvT1[PedEval2$Generation == 0])
    stFactor <- rbind(stFactor, c(rep, h2, "Program1", "Trait2", "Ebv",  "Sd", sd(PedEval1$EbvT2[PedEval1$Generation == 0])
    stFactor <- rbind(stFactor, c(rep, h2, "Program2", "Trait2", "Ebv",  "Sd", sd(PedEval2$EbvT2[PedEval2$Generation == 0])
                                  
                                  
    stFactor <- rbind(stFactor, c(rep, h2, "Program1", "Trait1", "Ebv",  "Mean", mean(PedEval1$EbvT1[PedEval1$Generation == 0]))
    stFactor <- rbind(stFactor, c(rep, h2, "Program2", "Trait1", "Ebv", "Mean",mean(PedEval2$EbvT1[PedEval2$Generation == 0]))
    stFactor <- rbind(stFactor, c(rep, h2, "Program1", "Trait2", "Ebv", "Mean",mean(PedEval1$EbvT2[PedEval1$Generation == 0]))
    stFactor <- rbind(stFactor, c(rep, h2, "Program2", "Trait2", "Ebv", "Mean",mean(PedEval2$EbvT1[PedEval2$Generation == 0]))
                    
                    
                    #Tbv
    stFactor <- rbind(stFactor, c(rep, h2, "Program1", "Trait1", "Tbv", "Sd",sd(PedEval1$TbvT1[PedEval1$Generation == 0])
    stFactor <- rbind(stFactor, c(rep, h2, "Program2", "Trait1", "Tbv", "Sd",sd(PedEval2$TbvT1[PedEval2$Generation == 0])
    stFactor <- rbind(stFactor, c(rep, h2, "Program1", "Trait2", "Tbv", "Sd",sd(PedEval1$TbvT2[PedEval1$Generation == 0])
    stFactor <- rbind(stFactor, c(rep, h2, "Program2", "Trait2", "Tbv", "Sd",sd(PedEval2$TbvT2[PedEval2$Generation == 0])
                                                                                                                                      
                                                                                                                                      
    stFactor <- rbind(stFactor, c(rep, h2, "Program1", "Trait1", "Tbv", "Mean", mean(PedEval1$TbvT1[PedEval1$Generation == 0]))
    stFactor <- rbind(stFactor, c(rep, h2, "Program2", "Trait1", "Tbv", "Mean",  mean(PedEval2$TbvT1[PedEval2$Generation == 0]))
    stFactor <- rbind(stFactor, c(rep, h2, "Program1", "Trait2", "Tbv",  "Mean", mean(PedEval1$TbvT2[PedEval1$Generation == 0]))
    stFactor <- rbind(stFactor, c(rep, h2, "Program2", "Trait2", "Tbv",  "Mean", mean(PedEval2$TbvT1[PedEval2$Generation == 0]))
                    
  }
  
}

write.table(stFactor, "StandFactors.csv", quote=FALSE, row.names=FALSE)
                      

                      
