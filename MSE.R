## this is a script to combine the heritability partitions
######################################################
#EBV
######################################################
library(AlphaPart)
library(mltools)

MSE <- data.frame(H2=NA, Program = NA, Trait = NA, Path = NA, Population = NA,  MSE = NA)


for (h2 in c(0.05, 0.1, 0.15, 0.2, 0.3, 0.35, 0.4, 0.5)) {
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
  
  
  #poreži živali iz nukleusa
  #Program 1
  #EBV
  Part1PN <- Part1
  Part1PN$EbvT1_s <- Part1PN$EbvT1_s[Part1PN$EbvT1_s$Program == "PN1",]
  Part1PN$EbvT2_s <- Part1PN$EbvT2_s[Part1PN$EbvT2_s$Program == "PN1",]
  Part1PN$EbvI_s <- Part1PN$EbvI_s[Part1PN$EbvI_s$Program == "PN1",]
  
  Part1GN <- Part1
  Part1GN$EbvT1_s <- Part1GN$EbvT1_s[Part1GN$EbvT1_s$Program == "GN",]
  Part1GN$EbvT2_s <- Part1GN$EbvT2_s[Part1GN$EbvT2_s$Program == "GN",]
  Part1GN$EbvI_s <- Part1GN$EbvI_s[Part1GN$EbvI_s$Program == "GN",]
  
  #TGV
  
  Part1gPN <- Part1g
  Part1gPN$TbvT1_s <- Part1gPN$TbvT1_s[Part1gPN$TbvT1_s$Program == "PN1",]
  Part1gPN$TbvT2_s <- Part1gPN$TbvT2_s[Part1gPN$TbvT2_s$Program == "PN1",]
  Part1gPN$TbvI_s <- Part1gPN$TbvI_s[Part1gPN$TbvI_s$Program == "PN1",]
  
  Part1gGN <- Part1g
  Part1gGN$TbvT1_s <- Part1gGN$TbvT1_s[Part1gGN$TbvT1_s$Program == "GN",]
  Part1gGN$TbvT2_s <- Part1gGN$TbvT2_s[Part1gGN$TbvT2_s$Program == "GN",]
  Part1gGN$TbvI_s <- Part1gGN$TbvI_s[Part1gGN$TbvI_s$Program == "GN",]

  
  #program 2
  #EBV
  Part2PN <- Part2
  Part2PN$EbvT2_s <- Part2PN$EbvT2_s[Part2PN$EbvT2_s$Program == "PN2",]
  Part2PN$EbvT2_s <- Part2PN$EbvT2_s[Part2PN$EbvT2_s$Program == "PN2",]
  Part2PN$EbvI_s <- Part2PN$EbvI_s[Part2PN$EbvI_s$Program == "PN2",]
  
  Part2GN <- Part2
  Part2GN$EbvT2_s <- Part2GN$EbvT2_s[Part2GN$EbvT2_s$Program == "GN",]
  Part2GN$EbvT2_s <- Part2GN$EbvT2_s[Part2GN$EbvT2_s$Program == "GN",]
  Part2GN$EbvI_s <- Part2GN$EbvI_s[Part2GN$EbvI_s$Program == "GN",]
  
  #TGV
  Part2gPN <- Part2g
  Part2gPN$TbvT2_s <- Part2gPN$TbvT2_s[Part2gPN$TbvT2_s$Program == "PN2",]
  Part2gPN$TbvT2_s <- Part2gPN$TbvT2_s[Part2gPN$TbvT2_s$Program == "PN2",]
  Part2gPN$TbvI_s <- Part2gPN$TbvI_s[Part2gPN$TbvI_s$Program == "PN2",]
  
  Part2gGN <- Part2g
  Part2gGN$TbvT2_s <- Part2gGN$TbvT2_s[Part2gGN$TbvT2_s$Program == "GN",]
  Part2gGN$TbvT2_s <- Part2gGN$TbvT2_s[Part2gGN$TbvT2_s$Program == "GN",]
  Part2gGN$TbvI_s <- Part2gGN$TbvI_s[Part2gGN$TbvI_s$Program == "GN",]
  
  
  #for each trait
  #program 1
  #for GN anumals
  T11gn = cbind(Part1GN$EbvT1_s, Part1gGN$TbvT1_s)
  colnames(T11gn) <- gsub("-", ".", colnames(T11gn))
  T12gn = cbind(Part1GN$EbvT2_s, Part1gGN$TbvT2_s)
  colnames(T12gn) <- gsub("-", ".", colnames(T12gn))
  #for PN anumals
  T11pn = cbind(Part1PN$EbvT1_s, Part1gPN$TbvT1_s)
  colnames(T11pn) <- gsub("-", ".", colnames(T11pn))
  T12pn = cbind(Part1PN$EbvT2_s, Part1gPN$TbvT2_s)
  colnames(T12pn) <- gsub("-", ".", colnames(T12pn))
  
  #########
  #program 2
  #for GN animals
  T21gn = cbind(Part2GN$EbvT1_s, Part2gGN$TbvT1_s)
  colnames(T21gn) <- gsub("-", ".", colnames(T21gn))
  T22gn = cbind(Part2GN$EbvT2_s, Part2gGN$TbvT2_s)
  colnames(T22gn) <- gsub("-", ".", colnames(T22gn))
  #for PN animals
  T21pn = cbind(Part2PN$EbvT1_s, Part2gPN$TbvT1_s)
  colnames(T21pn) <- gsub("-", ".", colnames(T21pn))
  T22pn = cbind(Part2PN$EbvT2_s, Part2gPN$TbvT2_s)
  colnames(T22pn) <- gsub("-", ".", colnames(T22pn))
  
  #mse
  # program 1 
  paths1 <- gsub("-", ".", unique(Part1$EbvT1_s$ProgramGender))
  for (path in paths1) {
    #GN animals
    T11gnmse <- mse(T11gn[[paste0("EbvT1_s_", path)]], T11gn[[paste0("TbvT1_s_", path)]])
    T12gnmse <- mse(T12gn[[paste0("EbvT2_s_", path)]], T12gn[[paste0("TbvT2_s_", path)]])
    MSE <- rbind(MSE, c(h2, "PN1", "1", path, "GN1", T11gnmse))
    MSE <- rbind(MSE, c(h2, "PN1", "2", path, "GN1", T12gnmse))
    
    #PN animals
    T11pnmse <- mse(T11pn[[paste0("EbvT1_s_", path)]], T11pn[[paste0("TbvT1_s_", path)]])
    T12pnmse <- mse(T12pn[[paste0("EbvT2_s_", path)]], T12pn[[paste0("TbvT2_s_", path)]])
    MSE <- rbind(MSE, c(h2, "PN1", "1", path, "PN1", T11pnmse))
    MSE <- rbind(MSE, c(h2, "PN1", "2", path, "PN1", T12pnmse))
  }
  
  # program 2
  paths2 <- gsub("-", ".", unique(Part2$EbvT1_s$ProgramGender))
  for (path in paths2) {
    #for GN animal
    T21gnmse <- mse(T21gn[[paste0("EbvT1_s_", path)]], T21gn[[paste0("TbvT1_s_", path)]])
    T22gnmse <- mse(T22gn[[paste0("EbvT2_s_", path)]], T22gn[[paste0("TbvT2_s_", path)]])
    MSE <- rbind(MSE, c(h2, "PN2", "1", path, "GN2", T21gnmse))
    MSE <- rbind(MSE, c(h2, "PN2", "2", path, "GN2", T22gnmse))    
    
    T21pnmse <- mse(T21pn[[paste0("EbvT1_s_", path)]], T21pn[[paste0("TbvT1_s_", path)]])
    T22pnmse <- mse(T22pn[[paste0("EbvT2_s_", path)]], T22pn[[paste0("TbvT2_s_", path)]])
    MSE <- rbind(MSE, c(h2, "PN2", "1", path, "PN2", T21pnmse))
    MSE <- rbind(MSE, c(h2, "PN2", "2", path, "PN2", T22pnmse))
  }


}
write.table(MSE, "MSE.csv", quote=FALSE, row.names=FALSE)
