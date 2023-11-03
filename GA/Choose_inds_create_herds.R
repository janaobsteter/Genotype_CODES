# Title     : TODO
# Objective : TODO
# Created by: jana
# Created on: 13.5.2019


#ustvari črede in zapiši tabele
ped <- read.table("PedigreeAndGeneticValues_cat.txt", header=TRUE)
pedK <- ped[ped$cat=="k",c("Indiv", "Mother", "gvNormUnres1", "phenoNormUnres1")]
pedK$cluster <- NA
pedK$cluster <- kmeans(pedK[,2:4], 100)$cluster
print("DOING Herds")
herdNO <- as.data.frame(table(pedK$cluster))
colnames(herdNO) <- c("Herd", "NoAnim")
write.csv(herdNO, "HerdNo.txt", quote=FALSE, row.names=FALSE)

print("Merging and writing DFs")
pedHerd <- merge(pedK[c("Indiv", "cluster")], ped, by="Indiv", all.x=TRUE)
write.table(pedK, "PedCows_HERDS.txt", quote=FALSE, row.names=FALSE)
write.table(pedHerd, "PedCows_HERDS_Total.txt", quote=FALSE, row.names=FALSE)


print("Selecting Inds")
#izberi posameznike za optimizacijo in jih zapiši
selInd <- ped[ped$cat %in% c('k', 'bm', 'pBM', 'pb', 'potomciNP'),]
write.table(selInd$Indiv, "INDPED.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)


