## this is a script to combine the heritability partitions
######################################################
#EBV
######################################################
library(AlphaPart)
library(mltools)


partCor <- data.frame(H2=NA, Program = NA, Trait = NA, Path = NA, Cor = NA)

for (h2 in c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5)) {
  #program 1
  PedEval1 <- read.table(paste0("H_", h2, "/PN1/PedEval1_", h2, ".csv"), header=TRUE)
  #program 2
  PedEval2 <- read.table(paste0("H_", h2, "/PN2/PedEval2_", h2, ".csv"), header=TRUE)
  
  # ---- Partitioning the trend PN1 ----
  PedEval1$EbvI = 0.5 * (PedEval1$EbvT1 + PedEval1$EbvT2)
  PedEval1$EbvT1_s <- (PedEval1$EbvT1 - mean(PedEval1$EbvT1[PedEval1$Generation == 0])) / sd(PedEval1$EbvT1[PedEval1$Generation == 0])
  PedEval1$EbvT2_s <- (PedEval1$EbvT2 - mean(PedEval1$EbvT2[PedEval1$Generation == 0])) / sd(PedEval1$EbvT2[PedEval1$Generation == 0])
  PedEval1$EbvI_s <- 0.5 * (PedEval1$EbvT1_s + PedEval1$EbvT2_s)
  
  # ---- Partitioning the trend PN2 ----
  PedEval2$EbvI = 0.5 * (PedEval2$EbvT1 + PedEval2$EbvT2)
  PedEval2$EbvT1_s <- (PedEval2$EbvT1 - mean(PedEval2$EbvT1[PedEval2$Generation == 0])) / sd(PedEval2$EbvT1[PedEval2$Generation == 0])
  PedEval2$EbvT2_s <- (PedEval2$EbvT2 - mean(PedEval2$EbvT2[PedEval2$Generation == 0])) / sd(PedEval2$EbvT2[PedEval2$Generation == 0])
  PedEval2$EbvI_s <- 0.5 * (PedEval2$EbvT1_s + PedEval2$EbvT2_s)
  
  PedEval1$ProgramGender = paste(PedEval1$Program, PedEval1$Gender, sep = "-")
  PedEval2$ProgramGender = paste(PedEval2$Program, PedEval2$Gender, sep = "-")
  PedEval1$GenerationProgram <- paste(PedEval1$Generation, PedEval1$Program, sep="-")
  PedEval2$GenerationProgram <- paste(PedEval2$Generation, PedEval2$Program, sep="-")
  Part1 = AlphaPart(x = as.data.frame(PedEval1), sort = FALSE,
                    colId = "IId", colFid = "FId", colMid = "MId",
                    colPath = "ProgramGender", colAGV = c("EbvT1_s", "EbvT2_s", "EbvI_s"))
  
  Part2 = AlphaPart(x = as.data.frame(PedEval2), sort = FALSE,
                    colId = "IId", colFid = "FId", colMid = "MId",
                    colPath = "ProgramGender", colAGV = c("EbvT1_s", "EbvT2_s", "EbvI_s"))
  
  
  ########################################################
  #TGV
  #######################################################
  # ---- Partitioning the trend PN1 ----
  PedEval1$TbvI = 0.5 * (PedEval1$TbvT1 + PedEval1$TbvT2)
  PedEval1$TbvT1_s <- (PedEval1$TbvT1 - mean(PedEval1$TbvT1[PedEval1$Generation == 0])) / sd(PedEval1$TbvT1[PedEval1$Generation == 0])
  PedEval1$TbvT2_s <- (PedEval1$TbvT2 - mean(PedEval1$TbvT2[PedEval1$Generation == 0])) / sd(PedEval1$TbvT2[PedEval1$Generation == 0])
  PedEval1$TbvI_s <- 0.5 * (PedEval1$TbvT1_s + PedEval1$TbvT2_s)
  
  # ---- Partitioning the trend PN2 ----
  PedEval2$TbvI = 0.5 * (PedEval2$TbvT1 + PedEval2$TbvT2)
  PedEval2$TbvT1_s <- (PedEval2$TbvT1 - mean(PedEval2$TbvT1[PedEval2$Generation == 0])) / sd(PedEval2$TbvT1[PedEval2$Generation == 0])
  PedEval2$TbvT2_s <- (PedEval2$TbvT2 - mean(PedEval2$TbvT2[PedEval2$Generation == 0])) / sd(PedEval2$TbvT2[PedEval2$Generation == 0])
  PedEval2$TbvI_s <- 0.5 * (PedEval2$TbvT1_s + PedEval2$TbvT2_s)
  
  
  PedEval1$ProgramGender = paste(PedEval1$Program, PedEval1$Gender, sep = "-")
  PedEval2$ProgramGender = paste(PedEval2$Program, PedEval2$Gender, sep = "-")
  Part1g = AlphaPart(x = as.data.frame(PedEval1), sort = FALSE,
                    colId = "IId", colFid = "FId", colMid = "MId",
                    colPath = "ProgramGender", colAGV = c("TbvT1_s", "TbvT2_s", "TbvI_s"))
  Part2g = AlphaPart(x = as.data.frame(PedEval2), sort = FALSE,
                    colId = "IId", colFid = "FId", colMid = "MId",
                    colPath = "ProgramGender", colAGV = c("TbvT1_s", "TbvT2_s", "TbvI_s"))
  
  
  
  
  Part11cor <- merge(Part1$EbvT1_s, Part1g$TbvT1_s, by=c("Generation", "IId", "FId", "MId", "Program", "ProgramGender"))
  colnames(Part11cor) <- gsub("-", ".", colnames(Part11cor))
  
  Part12cor <- merge(Part1$EbvT2_s, Part1g$TbvT2_s, by=c("Generation", "IId", "FId", "MId", "Program", "ProgramGender"))
  colnames(Part12cor) <- gsub("-", ".", colnames(Part12cor))
  
  #program 1
  paths1 <- unique(gsub("-", ".", Part11cor$ProgramGender))
  for (path in paths1) {
    #trait 1
    tmpcor <- cor(Part11cor[[paste0("EbvT1_s_", path)]], Part11cor[[paste0("TbvT1_s_", path)]])
    partCor <- rbind(partCor, c(h2, "PN1", "T1", path, tmpcor))
    #trait 2
    tmpcor <- cor(Part12cor[[paste0("EbvT2_s_", path)]], Part12cor[[paste0("TbvT2_s_", path)]])
    partCor <- rbind(partCor, c(h2, "PN1", "T2", path, tmpcor))
    
  }  
  
  
  Part21cor <- merge(Part2$EbvT1_s, Part2g$TbvT1_s, by=c("Generation", "IId", "FId", "MId", "Program", "ProgramGender"))
  colnames(Part21cor) <- gsub("-", ".", colnames(Part21cor))
  
  Part22cor <- merge(Part2$EbvT2_s, Part2g$TbvT2_s, by=c("Generation", "IId", "FId", "MId", "Program", "ProgramGender"))
  colnames(Part22cor) <- gsub("-", ".", colnames(Part22cor))
  
  #program 2
  paths2 <- unique(gsub("-", ".", Part21cor$ProgramGender))
  for (path in paths2) {
    #trait 1
    tmpcor <- cor(Part21cor[[paste0("EbvT1_s_", path)]], Part21cor[[paste0("TbvT1_s_", path)]])
    partCor <- rbind(partCor, c(h2, "PN2", "T1", path, tmpcor))
    #trait 2
    tmpcor <- cor(Part22cor[[paste0("EbvT2_s_", path)]], Part22cor[[paste0("TbvT2_s_", path)]])
    partCor <- rbind(partCor, c(h2, "PN2", "T2", path, tmpcor))
    
  }
}
  

write.table(partCor, paste0("AccuracyPart.csv"), quote=FALSE, row.names=FALSE)



