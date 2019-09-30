library(Rmisc)
sc = read.csv("~/Documents/Projects/inProgress/Phenotyping/Scenarios.csv")
update = read.csv("~/Documents/Projects/inProgress/Phenotyping/YearsUpdate.csv")
head(sc)
sc <- sc[,-(ncol(sc))]

scM <- melt(sc, id.vars = c("P_G", "Reference", "Population", "ReduceControls"))
updateM <- melt(update, id.vars = c("P_G", "Reference", "Population", "ReduceControls", "YearsToUpdate"))
scM$variable <- gsub(pattern = "X", replacement = "", x = scM$variable)
head(scM)
table(scM$variable)
scM$variable <- as.numeric(scM$variable)
colnames(scM)[5] <- "Year"
head(scM)

scM$Group <- paste(scM$P_G, scM$Population, scM$ReduceControls, sep="_")
updateM$Group <- paste(update$P_G, update$Population, update$ReduceControls, sep="_")
updateM$Year = 5
updateM$value = 50000
updateM$YearsToUpdate <- as.numeric(updateM$YearsToUpdate)
updateMP = updateM
updateMP$PercentageUpdate = 1 / updateMP$YearsToUpdate
updateMP$PercentageUpdate = round(updateMP$PercentageUpdate, 2)
updateMP$value <- 75000
#add reference
scM$value[scM$Reference == "yes" & scM$Population == "Total"] <- scM$value[scM$Reference == "yes" & scM$Population == "Total"] + 10000



yes11 <- ggplot(data=scM[scM$Reference=="yes" & scM$P_G=='1:1',], aes(x=Year, y=value, group=Group, colour=Population)) + 
  geom_path() + ylab("Number of animals") + 
  geom_hline(yintercept = 10000, colour="grey30") + 
  geom_hline(yintercept = 5000, colour="grey70") + 
  facet_grid(~ ReduceControls) + 
  theme_minimal() + 
  geom_text(data = updateM[updateM$P_G == "1:1",], label = updateM$YearsToUpdate[updateM$P_G == "1:1"], colour="black") +
  geom_text(data = updateMP[updateMP$P_G == "1:1",], label = updateMP$PercentageUpdate[updateMP$P_G == "1:1"], colour="forestgreen")


yes21 <- ggplot(data=scM[scM$Reference=="yes" & scM$P_G=='2:1',], aes(x=Year, y=value, group=Group, colour=Population)) + 
  geom_path() + ylab("Number of animals") + 
  geom_hline(yintercept = 10000, colour="grey30") + 
  geom_hline(yintercept = 5000, colour="grey70") + 
  facet_grid(~ ReduceControls) +   
  theme_minimal() + 
  geom_text(data = updateM[updateM$P_G == "2:1",], label = updateM$YearsToUpdate[updateM$P_G == "2:1"], colour="black") +
  geom_text(data = updateMP[updateMP$P_G == "2:1",], label = updateMP$PercentageUpdate[updateMP$P_G == "2:1"], colour="forestgreen")


yes12 <- ggplot(data=scM[scM$Reference=="yes" & scM$P_G=='1:2',], aes(x=Year, y=value, group=Group, colour=Population)) + 
  geom_path() + ylab("Number of animals") + 
  geom_hline(yintercept = 10000, colour="grey30") + 
  geom_hline(yintercept = 5000, colour="grey70") + 
  theme_minimal() + 
  facet_grid(~ ReduceControls) + 
  geom_text(data = updateM[updateM$P_G == "1:2",], label = updateM$YearsToUpdate[updateM$P_G == "1:2"], colour="black") +
  geom_text(data = updateMP[updateMP$P_G == "1:2",], label = updateMP$PercentageUpdate[updateMP$P_G == "1:2"], colour="forestgreen")


multiplot(yes11, yes21, yes12)


no11 <- ggplot(data=scM[scM$Reference=="no" & scM$P_G=='1:1',], aes(x=Year, y=value, group=Group, colour=Population)) + 
  geom_path() + ylab("Number of animals") + 
  geom_hline(yintercept = 10000, colour="grey30") + 
  geom_hline(yintercept = 5000, colour="grey70") + 
  facet_grid(~ ReduceControls)
no21 <- ggplot(data=scM[scM$Reference=="no" & scM$P_G=='2:1',], aes(x=Year, y=value, group=Group, colour=Population)) + 
  geom_path() + ylab("Number of animals") + 
  geom_hline(yintercept = 10000, colour="grey30") + 
  geom_hline(yintercept = 5000, colour="grey70") + 
  facet_grid(~ ReduceControls)
no12 <- ggplot(data=scM[scM$Reference=="no" & scM$P_G=='1:2',], aes(x=Year, y=value, group=Group, colour=Population)) + 
  geom_path() + ylab("Number of animals") + 
  geom_hline(yintercept = 10000, colour="grey30") + 
  geom_hline(yintercept = 5000, colour="grey70") + 
  facet_grid(~ ReduceControls)


multiplot(no11, no21, no12)



library(ggplot2)
p <- ggplot(mtcars, aes(mpg, wt)) + geom_point()
p <- p + facet_grid(. ~ cyl)

ann_text <- data.frame(mpg = 15,wt = 5,lab = "Text",
                       cyl = factor(8,levels = c("4","6","8")))
p + geom_text(data = ann_text,label = ann_text$lab)


######
#Price
price = read.csv("~/Documents/Projects/inProgress/Phenotyping/PriceControl.csv")
priceM <- melt(price, id.vars = "No.Controls")

ggplot(priceM, aes(x=No.Controls, y=value, group=variable, colour=variable)) + 
  geom_line() +
  scale_y_continuous(breaks = 1:20) + 
  scale_x_continuous(breaks = 1:11) + 
  ylab("Price [euro]") + 
  xlab("Number of recordings")

saved = read.csv("~/Documents/Projects/inProgress/Phenotyping/PG.csv")
saved$Group = paste0(saved$Remove..controls, saved$Phenotype.Genotype)
saved$Genotypes = as.numeric(as.character(saved$Genotypes))

ggplot(data=saved, aes(x=Remove..controls, y=Genotypes, group=Phenotype.Genotype, colour=Phenotype.Genotype)) + geom_line() + 
  xlab("Recordings removed")
