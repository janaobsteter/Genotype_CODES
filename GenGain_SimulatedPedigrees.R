pedClas <- read.csv("~/bin/AlphaSim1.05Linux/REAL20GenSel_Class/SimulatedData/PedigreeAndGeneticValues_cat.txt", sep=" ")
gengainM <- aggregate(pedClas$gvNormUnres1 ~pedClas$Generation, pedClas, function(x) (mean = mean(x)))
colnames(gengainM) <- c("Gen", "Mean")
gengainSD <- aggregate(pedClas$gvNormUnres1 ~pedClas$Generation, pedClas, function(x) (sd = sd(x)))
colnames(gengainSD) <- c("Gen", "SD")
gengain <- merge(gengainM, gengainSD, by="Gen")
gengain <- gengain[40:60,]
gengain$stMean <- (gengain$Mean - gengain$Mean[1]) / gengain$SD[1]
gengainClass <- gengain
gengainClass$scenario <- "Conventional"
lm(gengainClass$Mean ~ gengainClass$Gen)

edClas <- read.csv("~/bin/AlphaSim1.05Linux/REAL20GenSel_GenSLO/SimulatedData/PedigreeAndGeneticValues_cat.txt", sep=" ")
gengainM <- aggregate(pedClas$gvNormUnres1 ~pedClas$Generation, pedClas, function(x) (mean = mean(x)))
colnames(gengainM) <- c("Gen", "Mean")
gengainSD <- aggregate(pedClas$gvNormUnres1 ~pedClas$Generation, pedClas, function(x) (sd = sd(x)))
colnames(gengainSD) <- c("Gen", "SD")
gengain <- merge(gengainM, gengainSD, by="Gen")
gengain <- gengain[40:60,]
gengain$stMean <- (gengain$Mean - gengain$Mean[1]) / gengain$SD[1]
gengainGenSLO <- gengain
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
gengainGenSLO_BmGen <- gengain
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
gengainGenSplosnaPop <- gengain
gengainGenSplosnaPop$scenario <- "Genomic B"
lm(gengainGenSplosnaPop$Mean ~ gengainGenSplosnaPop$Gen)

pedClas <- read.csv("~/bin/AlphaSim1.05Linux/REAL20GenSel_Gen/SimulatedData/PedigreeAndGeneticValues_cat.txt", sep=" ")
gengainM <- aggregate(pedClas$gvNormUnres1 ~pedClas$Generation, pedClas, function(x) (mean = mean(x)))
colnames(gengainM) <- c("Gen", "Mean")
gengainSD <- aggregate(pedClas$gvNormUnres1 ~pedClas$Generation, pedClas, function(x) (sd = sd(x)))
colnames(gengainSD) <- c("Gen", "SD")
gengain <- merge(gengainM, gengainSD, by="Gen")
gengain <- gengain[40:60,]
gengain$stMean <- (gengain$Mean - gengain$Mean[1]) / gengain$SD[1]
gengainGen <- gengain
gengainGen$scenario <- "Genomic D"
lm(gengainGen$Mean ~ gengainGen$Gen)


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

########################################################
########################################################
GI <- data.frame()
genInt <- read.csv("~/bin/AlphaSim1.05Linux/REAL20GenSel_Class/GenInts.txt", sep=" ")
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
ggplot(GI, aes(x=Gen, y=genInt, colour=label)) + geom_path() +  
  xlab("Generation") + ylab("Generation interval [years]") + 
  scale_color_grey("", labels=c("dam>dam", "dam>sire", "sire>dam", "sire>sire")) +
  facet_grid(Scenario ~ ., scales = "free_y") + theme(legend.position = "bottom", legend.text=element_text(size=10)) + 
  theme(axis.title.x=element_text(margin=margin(15,0,0,0),size=12), axis.title.y=element_text(size=12), axis.text=element_text(size=11), legend.title=element_text(size=12))




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
