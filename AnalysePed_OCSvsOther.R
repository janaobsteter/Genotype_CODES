pedOCS <- read.table("~/Ped_OCS.txt", header=TRUE)
ped <- read.csv("~/Ped_Gen0.txt", header=TRUE, sep=" ")
ped1Pb <- read.csv("~/Ped_Gen0_1Pb.txt", header=TRUE, sep=" ")
ped1Year <- read.csv("~/Ped_Gen0_1Year.txt", header=TRUE, sep=" ")


pedOCS_a <- aggregate(pedOCS$gvNormRestr1 ~ pedOCS$Generation, FUN="mean")
ped_a <- aggregate(ped$gvNormRestr1 ~ ped$Generation, FUN="mean")
ped1Pb_a <- aggregate(ped1Pb$gvNormRestr1 ~ ped1Pb$Generation, FUN="mean")
ped1Year_a <- aggregate(ped1Year$gvNormRestr1 ~ ped1Year$Generation, FUN="mean")

pedOCS_a$strategy <- "OCS"
ped_a$strategy <- "SU 55"
ped1Pb_a$strategy <- "SU 15"
ped1Year_a$strategy <- "SU 51"

library(data.table)
Avg <- as.data.frame(rbindlist(list(pedOCS_a, ped_a, ped1Pb_a, ped1Year_a)))
colnames(Avg) <- c("Generation", "TGV", "Strategy")

library(ggplot2)
Avg$Generation <- as.numeric(Avg$Generation)
Avg$TGV <- as.numeric(Avg$TGV)
Avg$Strategy <- as.factor(Avg$Strategy)
ggplot(data=Avg, aes(x = Avg$Generation, y=Avg$TGV, group=Avg$Strategy, colour=Avg$Strategy, fill=Avg$Strategy)) + geom_path() + scale_colour_discrete() + 
  xlab("Generation") + ylab("Average True Genetic Value") + theme(legend.title = element_text("Strategy"))


ggplot(data=Avg[Avg$Strategy=="OCS",], aes(x = Avg$Generation[Avg$Strategy=="OCS"], y=Avg$TGV[Avg$Strategy=="OCS"]), group=Avg$Strategy[Avg$Strategy=="OCS"], colour=Avg$Strategy[Avg$Strategy=="OCS"], fill=Avg$Strategy[Avg$Strategy=="OCS"]) + geom_path() + scale_color_brewer()
ggplot(data=Avg[Avg$Strategy=="55",], aes(x = Avg$Generation[Avg$Strategy=="55"], y=Avg$TGV[Avg$Strategy=="55"]), group=Avg$Strategy[Avg$Strategy=="55"], colour=Avg$Strategy[Avg$Strategy=="55"], fill=Avg$Strategy[Avg$Strategy=="55"]) + geom_path() + scale_color_brewer()
ggplot(data=Avg[Avg$Strategy=="51",], aes(x = Avg$Generation[Avg$Strategy=="51"], y=Avg$TGV[Avg$Strategy=="51"]), group=Avg$Strategy[Avg$Strategy=="51"], colour=Avg$Strategy[Avg$Strategy=="51"], fill=Avg$Strategy[Avg$Strategy=="51"]) + geom_path() + scale_color_brewer()
ggplot(data=Avg[Avg$Strategy=="15",], aes(x = Avg$Generation[Avg$Strategy=="15"], y=Avg$TGV[Avg$Strategy=="15"]), group=Avg$Strategy[Avg$Strategy=="15"], colour=Avg$Strategy[Avg$Strategy=="15"], fill=Avg$Strategy[Avg$Strategy=="15"]) + geom_path() + scale_color_brewer()

varOCS <- read.table("~/Var_OCS.txt", header=TRUE)
var <- read.table("~/Var.txt", header=TRUE)
var1Pb <- read.table("~/Var_1Pb.txt", header=TRUE)
var1Year <- read.table("~/Var_1Year.txt", header=TRUE)

varOCS <- varOCS[varOCS$QtnModel==1, c(1,3)]
var <- var[var$QtnModel==1, c(1,3)]
var1Pb <- var1Pb[var1Pb$QtnModel==1, c(1,3)]
var1Year <- var1Year[var1Year$QtnModel==1, c(1,3)]

varOCS$strategy <- "OCS"
var$strategy <- "SU 55"
var1Pb$strategy <- "SU 15"
var1Year$strategy <- "SU 51"
 
AvgVar <- as.data.frame(rbindlist(list(varOCS, var, var1Pb, var1Year)))

ggplot(data=AvgVar, aes(x = AvgVar$Generation, y=AvgVar$AdditGenicVar1, group=AvgVar$strategy, colour=AvgVar$strategy, fill=AvgVar$strategy)) + geom_path() + scale_colour_discrete() + 
  xlab("Generation") + ylab("Additive genetic variance") + theme(legend.title = element_text("Strategy"))
