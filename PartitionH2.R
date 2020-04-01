## this is a script to combine the heritability partitions
######################################################
#EBV
######################################################
library(AlphaPart)
library(mltools)


for (rep in 0:9) {
for (h2 in c(0.25)) { ##c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 0.99)) {
  #program 1
  PedEval1 <- read.table(paste0("Rep", rep, "/H_", h2, "/PN1/PedEval1_", h2, ".csv"), header=TRUE)
  PedEval1 <- PedEval1[PedEval1$Generation %in% 20:41,]
  PedEval1$Program[PedEval1$Program == "BurnIn"] <- "GN"
  #program 2
  PedEval2 <- read.table(paste0("Rep", rep, "/H_", h2, "/PN2/PedEval2_", h2, ".csv"), header=TRUE)
  PedEval2 <- PedEval2[PedEval2$Generation %in% 20:41,]
  PedEval2$Program[PedEval2$Program == "BurnIn"] <- "GN"


  
  # ---- Partitioning the trend PN1 ----
  PedEval1$EbvI = 0.5 * (PedEval1$EbvT1 + PedEval1$EbvT2)

  PedEval1$EbvT1_s <- (PedEval1$EbvT1 - mean(PedEval1$EbvT1[PedEval1$Generation == 20])) / sd(PedEval1$TbvT1[PedEval1$Generation == 20])
  PedEval1$EbvT2_s <- (PedEval1$EbvT2 - mean(PedEval1$EbvT2[PedEval1$Generation == 20])) / sd(PedEval1$TbvT2[PedEval1$Generation == 20])
  PedEval1$EbvI_s <- 0.5 * (PedEval1$EbvT1_s + PedEval1$EbvT2_s)
  
  # ---- Partitioning the trend PN2 ----
  PedEval2$EbvI = 0.5 * (PedEval2$EbvT1 + PedEval2$EbvT2)
  PedEval2$EbvT1_s <- (PedEval2$EbvT1 - mean(PedEval2$EbvT1[PedEval2$Generation == 20])) / sd(PedEval2$TbvT1[PedEval2$Generation == 20])
  PedEval2$EbvT2_s <- (PedEval2$EbvT2 - mean(PedEval2$EbvT2[PedEval2$Generation == 20])) / sd(PedEval2$TbvT2[PedEval2$Generation == 20])
  PedEval2$EbvI_s <- 0.5 * (PedEval2$EbvT1_s + PedEval2$EbvT2_s)
  
  PedEval1$ProgramGender = paste(PedEval1$Program, PedEval1$Gender, sep = "-")
  PedEval2$ProgramGender = paste(PedEval2$Program, PedEval2$Gender, sep = "-")
  PedEval1$GenerationProgram <- paste(PedEval1$Generation, PedEval1$Program, sep="-")
  PedEval2$GenerationProgram <- paste(PedEval2$Generation, PedEval2$Program, sep="-")
  print("Partition EBVs")
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
#  PedEval1$TbvI = 0.5 * (PedEval1$TbvT1 + PedEval1$TbvT2)
#  for (pop in c("GN", "PN1")) {
#	  PedEval1$TbvT1_s[PedEval1$Program == pop] <- 
#		(PedEval1$TbvT1[PedEval1$Program == pop] - mean(PedEval1$TbvT1[PedEval1$Generation == 21 & PedEval1$Program == pop])) / sd(PedEval1$TbvT1[PedEval1$Generation == 21 & PedEval1$Program == pop])
#	  PedEval1$TbvT2_s[PedEval1$Program == pop] <- 
#		(PedEval1$TbvT2[PedEval1$Program == pop] - mean(PedEval1$TbvT2[PedEval1$Generation == 21 & PedEval1$Program == pop])) / sd(PedEval1$TbvT2[PedEval1$Generation == 21 & PedEval1$Program == pop])
 # }
  PedEval1$TbvT1_s <- (PedEval1$TbvT1 -  mean(PedEval1$TbvT1[PedEval1$Generation == 20])) / sd(PedEval1$TbvT1[PedEval1$Generation == 20])
  PedEval1$TbvT2_s <- (PedEval1$TbvT2 -  mean(PedEval1$TbvT2[PedEval1$Generation == 20])) / sd(PedEval1$TbvT2[PedEval1$Generation == 20])
  PedEval1$TbvI_s <- 0.5 * (PedEval1$TbvT1_s + PedEval1$TbvT2_s)
  
  # ---- Partitioning the trend PN2 ----
 # PedEval2$TbvI = 0.5 * (PedEval2$TbvT1 + PedEval2$TbvT2)
#  for (pop in c("GN", "PN2")) {
#          PedEval2$TbvT1_s[PedEval2$Program == pop] <- 
#		(PedEval2$TbvT1[PedEval2$Program == pop] - mean(PedEval2$TbvT1[PedEval2$Generation == 21 & PedEval2$Program == pop])) / sd(PedEval2$TbvT1[PedEval2$Generation == 21 & PedEval2$Program == pop])
#          PedEval2$TbvT2_s[PedEval2$Program == pop] <- 
#		(PedEval2$TbvT2[PedEval2$Program == pop] - mean(PedEval2$TbvT2[PedEval2$Generation == 21 & PedEval2$Program == pop])) / sd(PedEval2$TbvT2[PedEval2$Generation == 21 & PedEval2$Program == pop])
#  }
  PedEval2$TbvT1_s <- (PedEval2$TbvT1 -  mean(PedEval2$TbvT1[PedEval1$Generation == 20])) / sd(PedEval2$TbvT1[PedEval1$Generation == 20])
  PedEval2$TbvT2_s <- (PedEval2$TbvT2 -  mean(PedEval2$TbvT2[PedEval1$Generation == 20])) / sd(PedEval2$TbvT2[PedEval1$Generation == 20])
  PedEval2$TbvI_s <- 0.5 * (PedEval2$TbvT1_s + PedEval2$TbvT2_s)
  
  
  PedEval1$ProgramGender = paste(PedEval1$Program, PedEval1$Gender, sep = "-")
  PedEval2$ProgramGender = paste(PedEval2$Program, PedEval2$Gender, sep = "-")
  print("Partition genetic values")
  
  Part1g = AlphaPart(x = as.data.frame(PedEval1), sort = FALSE,
                    colId = "IId", colFid = "FId", colMid = "MId",
                    colPath = "ProgramGender", colAGV = c("TbvT1_s", "TbvT2_s", "TbvI_s"))
  Part2g = AlphaPart(x = as.data.frame(PedEval2), sort = FALSE,
                    colId = "IId", colFid = "FId", colMid = "MId",
                    colPath = "ProgramGender", colAGV = c("TbvT1_s", "TbvT2_s", "TbvI_s"))
  
  
  print("poreži živali iz nukleusa")
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
  print("porezi zivali iz nukleusa, PN2")

  #EBV
  print("TO je PN2, EBV")
  Part2PN <- Part2
  Part2PN$EbvT1_s <- Part2PN$EbvT1_s[Part2PN$EbvT1_s$Program == "PN2",]
  Part2PN$EbvT2_s <- Part2PN$EbvT2_s[Part2PN$EbvT2_s$Program == "PN2",]
  Part2PN$EbvI_s <- Part2PN$EbvI_s[Part2PN$EbvI_s$Program == "PN2",]
  
  print("TO je GN, EBV")
  Part2GN <- Part2
  print("THis is T1")
  Part2GN$EbvT1_s <- Part2GN$EbvT1_s[Part2GN$EbvT1_s$Program == "GN",]
  print("This is T2")
  Part2GN$EbvT2_s <- Part2GN$EbvT2_s[Part2GN$EbvT2_s$Program == "GN",]
  print("THis is index")
  Part2GN$EbvI_s <- Part2GN$EbvI_s[Part2GN$EbvI_s$Program == "GN",]
  
  #TGV
  print("To je PN1, TGV")
  Part2gPN <- Part2g
  Part2gPN$TbvT1_s <- Part2gPN$TbvT1_s[Part2gPN$TbvT1_s$Program == "PN2",]
  Part2gPN$TbvT2_s <- Part2gPN$TbvT2_s[Part2gPN$TbvT2_s$Program == "PN2",]
  Part2gPN$TbvI_s <- Part2gPN$TbvI_s[Part2gPN$TbvI_s$Program == "PN2",]
  
  print("To je GN, EBV")
  Part2gGN <- Part2g
  print("This is T1")
  Part2gGN$TbvT1_s <- Part2gGN$TbvT1_s[Part2gGN$TbvT1_s$Program == "GN",]
  print("THis is T2")
  Part2gGN$TbvT2_s <- Part2gGN$TbvT2_s[Part2gGN$TbvT2_s$Program == "GN",]
  print("THis is indwex")
  Part2gGN$TbvI_s <- Part2gGN$TbvI_s[Part2gGN$TbvI_s$Program == "GN",]
  
  
  #fparticije
  #program1
  #EBV
  print("THis are partitions")
  print("Gn1 sum")
  Part1GNSummary = summary(object = Part1GN, by = "Generation")
  print("pn1 sum")
  Part1PNSummary = summary(object = Part1PN, by = "Generation")
  #TGV
  print("GN1 sum TGV")
  Part1gGNSummary = summary(object = Part1gGN, by = "Generation")
  print("PN1 sum TGV")
  Part1gPNSummary = summary(object = Part1gPN, by = "Generation")

  #program2  
  #EBV
  print("GN2")
  Part2GNSummary = summary(object = Part2GN, by = "Generation")
  print("PN2")
  Part2PNSummary = summary(object = Part2PN, by = "Generation")
  #TGV
  print("GN2 TGV")
  Part2gGNSummary = summary(object = Part2gGN, by = "Generation")
  print("PN2 TGV")
  Part2gPNSummary = summary(object = Part2gPN, by = "Generation")
  

  #bind partitions
  #EbvT1, EbvT2 and EbvI
  #TbvT1, TbvT2 and TbvI for TGV
  partitionDF1 = data.frame()
  partitionDF2 = data.frame()
  print("Entering the loop")
  print("Entering the LOOP")
  for (trait in c("T1", "T2", "I")) {
    print(trait)
    #GN
    print("THis is GN, ebv")
    #Program 1
    t1 <- Part1GNSummary[[paste0("Ebv", trait, "_s")]]$abs
    t1$Program <- "PN1"
    t1$Trait <- trait
    t1$value <- "Ebv"
    t1$Population <- "GN1"
    partitionDF1 <- rbind(partitionDF1, t1)
    
    print("THis is GN, TGV")
    t1g <- Part1gGNSummary[[paste0("Tbv", trait, "_s")]]$abs
    t1g$Program <- "PN1"
    t1g$Trait <- trait
    t1g$value <- "Tbv"
    t1g$Population <- "GN1"
    partitionDF1 <- rbind(partitionDF1, t1g)
    
    #Program 2
    print("THis is GN2")
    t2 <- Part2GNSummary[[paste0("Ebv", trait, "_s")]]$abs
    t2$Program <- "PN2"
    t2$Trait <- trait
    t2$value <- "Ebv"
    t2$Population <- "GN2"
    partitionDF2 <- rbind(partitionDF2, t2)
    
    print("This is GN2, tgv")
    t2g <- Part2gGNSummary[[paste0("Tbv", trait, "_s")]]$abs
    t2g$Program <- "PN2"
    t2g$Trait <- trait
    t2g$value <- "Tbv"
    t2g$Population <- "GN2"
    partitionDF2 <- rbind(partitionDF2, t2g)
    
    rm(t1, t1g, t2, t2g)
    
    #PN
    #Program 1
    print("This is PN1, ebv")
    t1 <- Part1PNSummary[[paste0("Ebv", trait, "_s")]]$abs
    t1$Program <- "PN1"
    t1$Trait <- trait
    t1$value <- "Ebv"
    t1$Population <- "PN1"
    partitionDF1 <- rbind(partitionDF1, t1)
    
    print("This is PN1, tgv")
    t1g <- Part1gPNSummary[[paste0("Tbv", trait, "_s")]]$abs
    t1g$Program <- "PN1"
    t1g$Trait <- trait
    t1g$value <- "Tbv"
    t1g$Population <- "PN1"
    partitionDF1 <- rbind(partitionDF1, t1g)
    
    #Program 2
    print("This is PN2, ebv")
    t2 <- Part2PNSummary[[paste0("Ebv", trait, "_s")]]$abs
    t2$Program <- "PN2"
    t2$Trait <- trait
    t2$value <- "Ebv"
    t2$Population <- "PN2"
    partitionDF2 <- rbind(partitionDF2, t2)
    
    print("This is Pn2, tgv")
    t2g <- Part2gPNSummary[[paste0("Tbv", trait, "_s")]]$abs
    t2g$Program <- "PN2"
    t2g$Trait <- trait
    t2g$value <- "Tbv"
    t2g$Population <- "PN2"
    partitionDF2 <- rbind(partitionDF2, t2g)
    
  }
  write.table(partitionDF1, paste0("PartitionPN1_", h2, "_", rep, ".csv"), quote=FALSE, row.names=FALSE)
  write.table(partitionDF2, paste0("PartitionPN2_", h2, "_", rep, ".csv"), quote=FALSE, row.names=FALSE)
}
}


