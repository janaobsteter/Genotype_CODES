### AlphaPart.ped.R
###-----------------------------------------------------------------------------

AlphaPart.ped <- data.frame(IId = c("A", "B", "C", "T", "D", "E", "U", "V"),
                            FId = c("", "", "B", "B", "", "D", "D", "E"),
                            MId = c("", "", "A", "", "", "C", "", ""),
                            Gender = c("F", "M", "F", "F", "M", "M", "F", "F"),
                            Country = c("X", "Y", "C", "Y", "Y", "X", "Y", "X"),
                            Generation = c(1, 1, 2, 2, 2, 3, 3, 4))

###-----------------------------------------------------------------------------
### AlphaPart.ped.R ends here

save(AlphaPart.ped, file = "/home/jana/alphapartition/data/AlphaPart.ped.rda")

ped0 <- data.frame(id=c(1, 2, 3,  4, 5, 6, 7,  8, 9, 10),
                   fid=c(0, 0, 0,  1, 1, 1, 3,  4, 0, 5),
                   mid=c(0, 0, 0,  2, 0, 2, 2,  2, 5, 0),
                   by =c(NA, 0, 1, NA, 3, 3, 3, 3, 4, 4))

ped1 <- pedFixBirthYear(x=ped0, interval=1)
ped2 <- pedFixBirthYear(x=ped1, interval=1, down=TRUE)

pedSetBase (AlphaPart.ped, keep = AlphaPart.ped$Generation > 1)
