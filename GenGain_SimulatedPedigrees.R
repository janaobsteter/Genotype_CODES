<<<<<<< HEAD
pedClas <- read.csv("~/bin/AlphaSim1.05Linux/REAL20GenSel_Class/SimulatedData/PedigreeAndGeneticValues_cat.txt", sep=" ")
=======
pedClas <- read.csv("~/bin/AlphaSim1.05Linux/REAL20GenSel_Class1/SimulatedData/PedigreeAndGeneticValues_cat.txt", sep=" ")
>>>>>>> 1b8d6b61038ae92f43525513dcd93532f3369678
gengainM <- aggregate(pedClas$gvNormUnres1 ~pedClas$Generation, pedClas, function(x) (mean = mean(x)))
colnames(gengainM) <- c("Gen", "Mean")
gengainSD <- aggregate(pedClas$gvNormUnres1 ~pedClas$Generation, pedClas, function(x) (sd = sd(x)))
colnames(gengainSD) <- c("Gen", "SD")
gengain <- merge(gengainM, gengainSD, by="Gen")
gengain <- gengain[40:60,]
gengain$stMean <- (gengain$Mean - gengain$Mean[1]) / gengain$SD[1]
<<<<<<< HEAD
gengainClass <- gengain
gengainClass$scenario <- "Conventional"
lm(gengainClass$Mean ~ gengainClass$Gen)

edClas <- read.csv("~/bin/AlphaSim1.05Linux/REAL20GenSel_GenSLO/SimulatedData/PedigreeAndGeneticValues_cat.txt", sep=" ")
=======
gengain$stSD <- (gengain$SD / gengain$SD[1])
gengainClass <- gengain

gengainClass$scenario <- "Conventional"
lm(gengainClass$Mean ~ gengainClass$Gen)

pedClas <- read.csv("~/bin/AlphaSim1.05Linux/REAL20GenSel_GenSLO/SimulatedData/PedigreeAndGeneticValues_cat.txt", sep=" ")
>>>>>>> 1b8d6b61038ae92f43525513dcd93532f3369678
gengainM <- aggregate(pedClas$gvNormUnres1 ~pedClas$Generation, pedClas, function(x) (mean = mean(x)))
colnames(gengainM) <- c("Gen", "Mean")
gengainSD <- aggregate(pedClas$gvNormUnres1 ~pedClas$Generation, pedClas, function(x) (sd = sd(x)))
colnames(gengainSD) <- c("Gen", "SD")
gengain <- merge(gengainM, gengainSD, by="Gen")
gengain <- gengain[40:60,]
gengain$stMean <- (gengain$Mean - gengain$Mean[1]) / gengain$SD[1]
<<<<<<< HEAD
gengainGenSLO <- gengain
=======
gengain$stSD <- (gengain$SD / gengain$SD[1])
gengainGenSLO <- gengain

>>>>>>> 1b8d6b61038ae92f43525513dcd93532f3369678
gengainGenSLO$scenario <- "Genomic A"
lm(gengainGenSLO$Mean ~ gengainGenSLO$Gen)

pedClas <- read.csv("~/bin/AlphaSim1.05Linux/REAL20GenSel_GenSLO_BmGen/SimulatedData/PedigreeAndGeneticValues_cat.txt", sep=" ")
gengainM <- aggregate(pedClas$gvNormUnres1 ~pedClas$Generation, pedClas, function(x) (mean = mean(x)))
colnames(gengainM) <- c("Gen", "Mean")
gengainSD <- aggregate(pedClas$gvNormUnres1 ~pedClas$Generation, pedClas, function(x) (sd = sd(x)))
colnames(gengainSD) <- c("Gen", "SD")
gengain <- merge(gengainM, gengainSD, by="Gen")
gengain <- gengain[40:60,]
gengain$stMean <- (gengain$Mean - gengain$Mean[1]) / gengain$SD[1]
<<<<<<< HEAD
gengainGenSLO_BmGen <- gengain
=======
gengain$stSD <- (gengain$SD / gengain$SD[1])
gengainGenSLO_BmGen <- gengain

>>>>>>> 1b8d6b61038ae92f43525513dcd93532f3369678
gengainGenSLO_BmGen$scenario <- "Genomic C"
lm(gengainGenSLO_BmGen$Mean ~ gengainGenSLO_BmGen$Gen)

pedClas <- read.csv("~/bin/AlphaSim1.05Linux/REAL20GenSel_GenSplosnaPop/SimulatedData/PedigreeAndGeneticValues_cat.txt", sep=" ")
gengainM <- aggregate(pedClas$gvNormUnres1 ~pedClas$Generation, pedClas, function(x) (mean = mean(x)))
colnames(gengainM) <- c("Gen", "Mean")
gengainSD <- aggregate(pedClas$gvNormUnres1 ~pedClas$Generation, pedClas, function(x) (sd = sd(x)))
colnames(gengainSD) <- c("Gen", "SD")
gengain <- merge(gengainM, gengainSD, by="Gen")
gengain <- gengain[40:60,]
gengain$stMean <- (gengain$Mean - gengain$Mean[1]) / gengain$SD[1]
<<<<<<< HEAD
gengainGenSplosnaPop <- gengain
gengainGenSplosnaPop$scenario <- "Genomic B"
lm(gengainGenSplosnaPop$Mean ~ gengainGenSplosnaPop$Gen)

pedClas <- read.csv("~/bin/AlphaSim1.05Linux/REAL20GenSel_Gen/SimulatedData/PedigreeAndGeneticValues_cat.txt", sep=" ")
=======
gengain$stSD <- (gengain$SD / gengain$SD[1])
gengainGenSplosnaPop <- gengain

gengainGenSplosnaPop$scenario <- "Genomic B"
lm(gengainGenSplosnaPop$Mean ~ gengainGenSplosnaPop$Gen)

pedClas <- read.csv("~/bin/AlphaSim1.05Linux/REAL20GenSel_Gen2/SimulatedData/PedigreeAndGeneticValues_cat.txt", sep=" ")
>>>>>>> 1b8d6b61038ae92f43525513dcd93532f3369678
gengainM <- aggregate(pedClas$gvNormUnres1 ~pedClas$Generation, pedClas, function(x) (mean = mean(x)))
colnames(gengainM) <- c("Gen", "Mean")
gengainSD <- aggregate(pedClas$gvNormUnres1 ~pedClas$Generation, pedClas, function(x) (sd = sd(x)))
colnames(gengainSD) <- c("Gen", "SD")
gengain <- merge(gengainM, gengainSD, by="Gen")
gengain <- gengain[40:60,]
gengain$stMean <- (gengain$Mean - gengain$Mean[1]) / gengain$SD[1]
<<<<<<< HEAD
=======
gengain$stSD <- (gengain$SD / gengain$SD[1])
>>>>>>> 1b8d6b61038ae92f43525513dcd93532f3369678
gengainGen <- gengain
gengainGen$scenario <- "Genomic D"
lm(gengainGen$Mean ~ gengainGen$Gen)


<<<<<<< HEAD
pedClas <- read.csv("~/bin/AlphaSim1.05Linux//SimulatedData/PedigreeAndGeneticValues_cat.txt", sep=" ")
gengainM <- aggregate(pedClas$gvNormUnres1 ~pedClas$Generation, pedClas, function(x) (mean = mean(x)))
colnames(gengainM) <- c("Gen", "Mean")
gengainSD <- aggregate(pedClas$gvNormUnres1 ~pedClas$Generation, pedClas, function(x) (sd = sd(x)))
colnames(gengainSD) <- c("Gen", "SD")
gengain <- merge(gengainM, gengainSD, by="Gen")
gengain <- gengain[40:60,]
gengain$stMean <- (gengain$Mean - gengain$Mean[1]) / gengain$SD[1]
gengainGen2 <- gengain
gengainGen2$scenario <- "Genomic D2"
lm(gengainGen$Mean ~ gengainGen$Gen)


genGain <- rbind(gengainClass, gengainGenSLO, gengainGenSplosnaPop, gengainGenSLO_BmGen, gengainGen, gengainGen2)
ggplot(genGain, aes(x=genGain$Gen, y=genGain$stMean, colour=genGain$scenario)) + geom_path() + scale_color_grey("") + 
  xlab("Generation") + ylab("Average true genetic value") 
=======


genGain$scenario <- as.character(genGain$scenario)
genGain <- rbind(gengainClass, gengainGenSLO, gengainGenSplosnaPop, gengainGenSLO_BmGen, gengainGen)
genGainPlot <- ggplot(genGain, aes(x=genGain$Gen, y=genGain$stMean, group=genGain$scenario, colour=scenario, linetype=scenario)) + geom_line(aes(linetype=genGain$scenario), size=1) + 
  xlab("Generation") + ylab("True genetic value")  + scale_linetype_manual("Scenario", values=c("solid", "dashed", "dotted", "dotdash", "twodash"), labels=c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D")) + 
  scale_colour_manual("Scenario", values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1"), labels=c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D")) + 
  ggtitle("a")


lm <- ggplot(data = genGain, aes(x=genGain$stSD, y=genGain$stMean, colour=genGain$scenario)) + geom_path() + scale_x_reverse() + 
  geom_smooth(method='lm', se=FALSE) + 
  scale_color_hue("Shema", labels=c("Conventional", "GenomicSLO", "GenBulls on Other Cows", "GenBulls on Bull Dams", "GenBulls on All Cows")) + 
  xlab("Genic sd") + ylab("Mean genetic gain") + ggtitle("Genic variance standardisation")


ggplot(data = genGain, aes(x=genGain$Gen, y=genGain$stSD, colour=genGain$scenario, linetype=scenario)) +  geom_line(aes(linetype=genGain$scenario), size=1) + 
  xlab("Generation") + ylab("True genetic value")  + scale_linetype_manual("Scenario", values=c("solid", "dashed", "dotted", "dotdash", "twodash"), labels=c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D")) + 
  scale_colour_manual("Scenario", values=c("forestgreen", "dodgerblue2", "purple", "red3", "orange1"), labels=c("Conventional", "Genomic A", "Genomic B", "Genomic C", "Genomic D")) + 
  xlab("Generation") + ylab("Genetic SD") 
>>>>>>> 1b8d6b61038ae92f43525513dcd93532f3369678

########################################################
########################################################
GI <- data.frame()
<<<<<<< HEAD
genInt <- read.csv("~/bin/AlphaSim1.05Linux/REAL20GenSel_Class/GenInts.txt", sep=" ")
=======
genInt <- read.csv("~/bin/AlphaSim1.05Linux/REAL20GenSel_Class1/GenInts.txt", sep=" ")
>>>>>>> 1b8d6b61038ae92f43525513dcd93532f3369678
genInt$label <- paste0(genInt$line, genInt$sex)
theme_set(theme_gray(base_size = 14))
plotClas <- ggplot(genInt, aes(x=Gen, y=genInt, group=label, colour=label)) + geom_path() + xlab("Generation") + ylab("Generation interval [years]") + 
  scale_color_grey("Selection path", labels=c("dam>dam", "dam>sire", "sire>dam", "sire>sire")) + ggtitle("Conventional")
genInt$Scenario <- "Class"
GI <- rbind(GI, genInt)

genInt <- read.csv("~/bin/AlphaSim1.05Linux/REAL20GenSel_GenSLO/GenInts.txt", sep=" ")
genInt$label <- paste0(genInt$line, genInt$sex)
theme_set(theme_gray(base_size = 14))
plotGenSLO <- ggplot(genInt, aes(x=Gen, y=genInt, group=label, colour=label)) + geom_path() + xlab("Generation") + ylab("Generation interval [years]") + 
  scale_color_grey("Selection path", labels=c("dam>dam", "dam>sire", "sire>dam", "sire>sire")) + ggtitle("Genomic A")
genInt$Scenario <- "Genomic A"
GI <- rbind(GI, genInt)

genInt <- read.csv("~/bin/AlphaSim1.05Linux/REAL20GenSel_GenSplosnaPop/GenInts.txt", sep=" ")
genInt$label <- paste0(genInt$line, genInt$sex)
theme_set(theme_gray(base_size = 14))
plotSplosnaPop <- ggplot(genInt, aes(x=Gen, y=genInt, group=label, colour=label)) + geom_path() + xlab("Generation") + ylab("Generation interval [years]") + 
  scale_color_grey("Selection path", labels=c("dam>dam", "dam>sire", "sire>dam", "sire>sire")) + ggtitle("Genomic B")
genInt$Scenario <- "Genomic B"
GI <- rbind(GI, genInt)

genInt <- read.csv("~/bin/AlphaSim1.05Linux/REAL20GenSel_GenSLO_BmGen/GenInts.txt", sep=" ")
genInt$label <- paste0(genInt$line, genInt$sex)
theme_set(theme_gray(base_size = 14))
plotBmGen <- ggplot(genInt, aes(x=Gen, y=genInt, group=label, colour=label)) + geom_path() + xlab("Generation") + ylab("Generation interval [years]") + 
  scale_color_grey("Selection path", labels=c("dam>dam", "dam>sire", "sire>dam", "sire>sire")) +  ggtitle("Genomic C")
genInt$Scenario <- "Genomic C"
GI <- rbind(GI, genInt)

genInt <- read.csv("~/bin/AlphaSim1.05Linux/REAL20GenSel_Gen/GenInts.txt", sep=" ")
genInt$label <- paste0(genInt$line, genInt$sex)
theme_set(theme_gray(base_size = 14))
plotGen <- ggplot(genInt, aes(x=Gen, y=genInt, group=label, colour=label)) + geom_path() + xlab("Generation") + ylab("Generation interval [years]") + 
  scale_color_grey("Selection path", labels=c("dam>dam", "dam>sire", "sire>dam", "sire>sire")) + ggtitle("Genomic D")
genInt$Scenario <- "Genomic D"
GI <- rbind(GI, genInt)


grid_arrange_shared_legend(plotClas, plotGenSLO, plotBmGen, plotSplosnaPop, plotGen, ncol=2, nrow=3)
GI <- GI[GI$Gen %in% 40:60,]
<<<<<<< HEAD
ggplot(GI, aes(x=Gen, y=genInt, colour=label)) + geom_path() +  
  xlab("Generation") + ylab("Generation interval [years]") + 
  scale_color_grey("", labels=c("dam>dam", "dam>sire", "sire>dam", "sire>sire")) +
  facet_grid(Scenario ~ ., scales = "free_y") + theme(legend.position = "bottom", legend.text=element_text(size=10)) + 
  theme(axis.title.x=element_text(margin=margin(15,0,0,0),size=12), axis.title.y=element_text(size=12), axis.text=element_text(size=11), legend.title=element_text(size=12))


=======
GIPlot <- ggplot(GI, aes(x=Gen, y=genInt, colour=label, linetype =label, group=label)) + geom_line(aes(linetype=label), size=1) + 
 scale_linetype_manual("", values=c("solid", "dashed", "dotted", "dotdash"), labels=c("dam>dam", "dam>sire", "sire>dam", "sire>sire")) + 
  scale_colour_manual("", values=c("palevioletred2", "purple2", "royalblue3", "darkolivegreen4"), labels=c("dam>dam", "dam>sire", "sire>dam", "sire>sire"))+ ggtitle("b") +
  xlab("Generation") + ylab("Generation interval [years]") + xlim(c(40,60)) +
  facet_grid(Scenario ~ ., scales = "free_y") + theme(legend.position = "bottom", legend.text=element_text(size=12)) + 
  theme(axis.title.x=element_text(margin=margin(15,0,0,0),size=14), axis.title.y=element_text(size=14), axis.text=element_text(size=11), legend.title=element_text(size=12))


#TO JE ZDJ COMBINED
library(gridExtra)
library(grid)
gl <- lapply(1:1, function(ii) grobTree(rectGrob(), textGrob(ii)))

grid.arrange(a, genGainPlot, GIPlot, layout_matrix = rbind(c(1,3),
                                               c(2,3),
                                               c(2,3),
                                               c(2,3),
                                               c(1,3)))
>>>>>>> 1b8d6b61038ae92f43525513dcd93532f3369678


library(ggplot2)
library(gridExtra)
library(grid)
grid_arrange_shared_legend <- function(..., nrow = 1, ncol = length(list(...)), position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
  gl <- c(gl, nrow = nrow, ncol = ncol)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
  
}
