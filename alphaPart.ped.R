### alphaPart.ped.R


AlphaPart.ped <- data.frame(IId = c("A", "B", "C", "T", "D", "E", "U", "V"),
                            FId = c("", "", "B", "B", "", "D", "D", "E"),
                            MId = c("", "", "A", "", "", "C", "", ""),
                            Generation = c(1, 1, 2, 2, 2, 3, 3, 4),
                            Country = c("X", "Y", "X", "Y", "Y", "X", "Y", "X"),
                            Gender = c("F", "M", "F", "F", "M", "M", "F", "F"),
                            EBV = c(100, 105, 104, 102, 108, 107, 107, 109))

save(AlphaPart.ped, file = "/home/jana/AlphaPartition/alphaPart/data/AlphaPart.ped.rda")
###-----------------------------------------------------------------------------
### alphaPart.ped.R ends here
