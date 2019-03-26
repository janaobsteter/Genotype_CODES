#!/usr/bin/env Rscript

## This script defines the origin of the haplotype

library(ggplot2)
library(reshape)

## This loop needs to be adjusted with the other script (look at previous 8_HaplypeStatistics.sh)
setwd("~")

Base <- basename(getwd())

## Set core and breed working directory
BaseCore <- basename(getwd())
setwd("~")

BaseBreed <- basename(getwd())
setwd(BaseCore)

## Read in the data
generation           <- read.table("PedGen.txt", header=TRUE, stringsAsFactors=FALSE)
generation <- generation[generation$Generation %in% 40:60,]
Pedig                <- generation[,c("Indiv", "Father", "Mother")]





Bulls <- read.table("~/PedGenCat.txt", header=TRUE)
Bulls <- Bulls[Bulls$Generation %in% 35:60,]
Bulls <- Bulls[Bulls$cat %in% c("gpb", "pb"),]

## For the oldest individual plot a pedigree following the blog post from Susan Johnson
## https://susanejohnston.wordpress.com/2014/10/10/step-by-step-a-stylised-pedigree-using-reshape-and-ggplot2/
ped <- generation
pedigree           <- ped[,c("Indiv", "Father", "Mother", "Generation")]
names(pedigree)    <- c("ID", "FATHER", "MOTHER", "BirthYear")
pedigree[pedigree==0] <- NA
ped2 <- melt(pedigree, id.vars=c("ID"), measure.vars=c("MOTHER", "FATHER"))
#ped2 <- ped2[-which(is.na(ped2$value)),]
names(ped2)[which(names(ped2) == "value")] <- "Parent.ID"
ped2$Group <- 1:nrow(ped2)
ped3 <- melt(ped2, id.vars = "Group", measure.vars=c("ID", "Parent.ID"))
names(ped3)[3] <- "ID"
ped3 <- merge(ped3, pedigree[,c("ID", "BirthYear")], by="ID")
  
  generateXcoord <- function(size) {
    
    if(size %% 2 != 0 & size != 1) {   # Check if size is odd
      newsize <- size - 1
      interval <- 1/newsize
      x <- seq(0, 1, interval)
    }
    
    if(size %% 2 == 0) {    # Check if size is even
      interval <- 1/size
      x <- seq(0, 1, interval)[-size-1] + diff(seq(0, 1, interval))/2
    }
    
    if(size == 1)
      x <- 0.5
    x
  }
  
  #~~ Create an object to save the x coordinates within
  xcoords <- NULL
  
  for(j in unique(ped3$BirthYear)){
    
    # Extract the number of Unique IDs per year and generate X coords
    ids  <- unique(ped3$ID[which(ped3$BirthYear == j)])
    newx <- generateXcoord(length(ids)) # generate X coordinates
    
    # Append to xcoords
    xcoords <- rbind(xcoords,
                     data.frame(ID = ids,
                                x = sample(newx, size = length(newx), replace = F)))
    
    rm(ids, newx)
  }
  
  # Merge with ped3
  ped3 <- merge(ped3, xcoords)
  
## Extract only carrier individuals and keep only connections between carrier individuals
ped4 <- ped3[ped3$ID %in% Bulls$Indiv, ]

birthyears <- sort(unique(ped3$BirthYear[!is.na(ped3$BirthYear)]))

## Calculate number of carrier progeny for carrier parents
BullPedig  <- Pedig[Pedig$Indiv %in% Bulls$Indiv, ]
BullPedig[!(BullPedig$Father %in% Bulls$Indiv), "Father"] <- 0
BullPedig[!(BullPedig$Father %in% Bulls$Indiv), "Mother"] <- 0
BullProgeny <- data.frame(table(c(BullPedig$Father, BullPedig$Mother)))
CarrierProgeny <- merge(CarrierProgeny, Carriers[,c("V1","ITT")], by.x="Var1", by.y="V1")
colnames(BullProgeny)[1] <- "ID"
BullProgeny <- merge(BullProgeny, unique(ped4[,c("ID","x","BirthYear")]), by="ID")

try(CarrierProgeny[CarrierProgeny$x < 0.20, "x"] <- 0.20, silent=TRUE)
try(CarrierProgeny[CarrierProgeny$x > 0.80, "x"] <- 0.80, silent=TRUE)
## Extract only those with more than 50 carrier progeny
BullProgeny <- BullProgeny[BullProgeny$Freq > 50 & BullProgeny$Var1 !=0, ]
  
## Plot Complete pedigree
pdf(file="PedigrePlotComplete.pdf", width=13/2.54,height=20/2.54, pointsize=12)
print(ggplot(ped3, aes(x, -BirthYear)) +
        geom_line(aes(group = Group), alpha = 0.1) +
        geom_point(size = 0.5) +
        geom_point(size = 0.1) +
        theme(legend.position  = "none",
              axis.text.x      = element_blank(),
              axis.text.y      = element_text(colour = "darkgrey"),
              axis.ticks.y     = element_blank(),
              axis.ticks.x     = element_blank(),
              panel.grid       = element_blank(),
              plot.background  = element_rect(),
              panel.background = element_blank()) +
        theme(plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), units = "cm")) +
        scale_y_continuous(breaks = -seq(min(birthyears), max(birthyears), 1),
                           labels =  seq(min(birthyears), max(birthyears), 1)) +
        labs(x = "", y = "Birth year"))
dev.off()

## Plot Complete pedigree and carriers
pdf(file="PedigrePlotCompleteCarriers.pdf", width=13/2.54,height=20/2.54, pointsize=12)
print(ggplot(ped3, aes(x, -BirthYear)) +
        geom_line(aes(group = Group), alpha = 0.1) +
        geom_point(size = 0.1) +
        theme(plot.title       = element_text(hjust = 0.5),
              legend.position  = "none",
              axis.text.x      = element_blank(),
              axis.text.y      = element_text(colour = "darkgrey"),
              axis.ticks.y     = element_blank(),
              axis.ticks.x     = element_blank(),
              panel.grid       = element_blank(),
              plot.background  = element_rect(),
              panel.background = element_blank()) +
        theme(plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), units = "cm")) +
        scale_y_continuous(breaks = -seq(min(birthyears), max(birthyears), 1),
                           labels =  seq(min(birthyears), max(birthyears), 1)) +
        labs(x = "", y = "Birth year") +
        geom_point(data = ped4, color='red', size = 0.5) +
        geom_line(data = ped4, aes(group = Group), color='red', alpha = 0.1, size=1) +
        geom_point(data = ped4, size = 0.1))
dev.off()

## Plot pedigree for carriers only
pdf(file="PedigrePlotCarriers.pdf", width=13/2.54,height=20/2.54, pointsize=12)
print(ggplot(ped4, aes(x, -BirthYear)) +
        geom_line(aes(group = Group), color='red', alpha = 0.1) +
        geom_point(size = 0.5, color='red') +
        geom_point(size = 0.1) +
        theme(legend.position  = "none",
              axis.text.x      = element_blank(),
              axis.text.y      = element_text(colour = "darkgrey"),
              axis.ticks.y     = element_blank(),
              axis.ticks.x     = element_blank(),
              panel.grid       = element_blank(),
              plot.background  = element_rect(),
              panel.background = element_blank()) +
        theme(plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), units = "cm")) +
        scale_y_continuous(breaks = -seq(min(birthyears), max(birthyears), 1),
                           labels =  seq(min(birthyears), max(birthyears), 1)) +
        labs(x = "", y = "Birth year") +
        geom_line(aes(group = Group), color='red', alpha = 0.1))
dev.off()

