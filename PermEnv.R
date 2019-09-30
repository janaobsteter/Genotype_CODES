################
#Funkcija summarySE
#####################
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  colnames(datac)[colnames(datac)=="measurevar"] <- "mean"
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
#######################################################
#######################################################
ped <- read.csv("/home/jana/TGVsAll_permEnv_SU55_02082019.csv", sep=" ")

head(ped)

library(ggplot2)
ped$PlotGroup <- paste0(ped$scenario, ped$Repeats, "_", ped$Rep)
ped$PlotGroup1 <- paste0(ped$scenario, ped$Repeats)
head(ped)
ped$Repeats <- as.factor(ped$Repeats)
ggplot(data = ped, aes(x=Generation, y=zMean, colour=Repeats, linetype = scenario, group = PlotGroup)) + geom_line()

pedA <- summarySE(data = ped, measurevar = 'zMean', groupvars = c('Generation', 'scenario', 'Repeats'))
head(pedA)
pedA$PlotGroup <- paste0(pedA$scenario, pedA$Repeats)
pedA$max = pedA$mean + pedA$ci
pedA$min = pedA$mean - pedA$ci
ggplot(data=pedA, aes(x=Generation, y = mean, group=PlotGroup, colour=Repeats, linetype=scenario)) + geom_line() +
       geom_ribbon(aes(ymin=min, ymax=max,x=Generation, colour=Repeats),  alpha=0.3)

head(pedA)

# add economic values
pedAA = summarySE(data=pedA, measurevar = 'mean', groupvars = c('scenario', 'Repeats'))[,c(1,2,4,5)]
pedAA$NoLact = 10825
pedAA$NoGeno[pedAA$scenario == "Class"] = 0
pedAA$NoGeno[pedAA$scenario == "Gen"] = 2050

## price per unit of genetic gain


##
nz = data.frame(Control = c(1, 2, 3, 4, 5, 8, 9, 11), PerATcontrol = c(5.49, 5.49, 4.91, 4.37, 4.04, 3.37, 3.30, 3.13)*0.589,
                PerAcontrol = c(5.69, 5.69, 5.11, 4.57, 4.24, 3.57, 3.5, 3.33)*0.589)

library(reshape)
nzM = melt(nz, id.vars = "Control")
nzM$valueE = nzM$value * 0.589
ggplot(data=nzM, aes(x=Control, y=valueE, group=variable, colour=variable)) + geom_line()

nzM$valueCumE <- nzM$valueE * nzM$Control
ggplot(data=nzM, aes(x=Control, y=valueCumE, group=variable, colour=variable)) + geom_line()

#the price of milk ontrol
library(ggplot2)

ir <- data.frame(Control = c(11, 8, 6, 4), Price = c(17, 14, 12, 10))
qplot(data=ir, x=Control, y=Price, geom="point")
ggplot(data=ir, aes(x=Control, y=Price)) + geom_line() + geom_point()

model = lm(Price ~ Control, data=ir)
pir = ir
pir$Price = ''
pir$Control = c(1,2,3,4)

predict(model, pir)
#number of animals with genotype per year = 10,825

#number of genotpes per year = ~2,050
ir$perAcontrol = (ir$Price) / ir$Control
ir$Cum <- ir$perAcontrol * ir$Control
qplot(data=ir, x=Control, y=perAcontrol, geom="line")
qplot(data=ir, x=Control, y=Cum, geom="line")

ir2 <- data.frame(Control = c(4,5,7,6,8,10,11), 
                  PriceDIYT = c(11.5, 13.5 ,NA, 15., 19, 22, NA),
                  PriceTechnician = c(11.5, 12.5, 13.5, NA, NA, 16, 19),
                  HerdFeeManual = 60, DiscketteFee = 33)
ir2$PriceDIY = ir2$PriceDIYT / ir2$Control
ir2$TPriceTech = (ir2$PriceTechnician*100 + ir2$HerdFeeManual + ir2$DiscketteFee *ir2$Control)/(100 * ir2$Control)

nz$Provider = "NZCrv"
colnames(nz)[2:3] = c("ATPrice", "APrice")
nzM = melt(nz, id.vars = c("Control", "Provider"))
head(ir)
irnew = ir[,c(1,3)]
colnames(irnew) = c("Control", "APrice")
irnew$Provider = "ICBF"
irnewM = melt(irnew, id.vars = c("Control", "Provider"))

ir2new = ir2[,c(1,6, 7)]
ir2new$Provider = "MU Ireland"
colnames(ir2new)[3] = "PriceTechnician"
ir2newM = melt(ir2new, id.vars = c("Control", "Provider"))

all = rbind(nzM, irnewM, ir2newM)
all$Type = paste(all$Provider, all$variable, sep="_")

all <- all[!is.na(all$value),]
all$Control <- as.numeric(all$Control)
all$value <- as.numeric(all$value)
ggplot(data=all, aes(x=Control, y=value, group=Type, colour=Type)) + geom_path()


x <-  aggregate(all$value ~all$Control, FUN="mean")
colnames(x) <- c("Control", "Value")
x$Rel <- NA
x$Pot = NA
for (row in 2: nrow(x)) {
  print(row)
  x$Rel[row] <- x$Value[row] / x$Value[row-1]
  x$Pot[row] = x$Value[row-1] * 0.96

}
x
sum(x$Value - x$Pot, na.rm=TRUE)

x$Pot = x$Value * 0.95
xM = melt(x[, c(1,2,4)], id.vars = "Control")
xM$Control[xM$variable == "Pot"] = xM$Control[xM$variable == "Pot"] +1

ggplot(data=xM, aes(x=Control, y=value, colour=variable, group=variable)) + geom_line()
