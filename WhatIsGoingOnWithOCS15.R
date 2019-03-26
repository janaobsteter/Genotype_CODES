ped <- read.table("PedGen315.txt", header = TRUE)

ped <- ped[ped$Generation %in% 40:60,]
head(ped)

fathers <- ped[ped$Indiv %in% unique(ped$Father),]


Gen <- summarySE(data=ped, measurevar = "gvNormUnres1", groupvars = "Generation")
fathers$Generation <- as.factor(fathers$Generation)
Gen$Generation <- as.factor(Gen$Generation)

ggplot() + geom_bar(data=Gen, aes(x=Generation, y=gvNormUnres1), stat="identity") + 
  geom_bar(data=fathers, aes(x=Generation, y=gvNormUnres1, group=Indiv, fill=Generation), stat="identity", position="dodge", width=0.3)

#očetje generacije 46
ped[ped$Indiv %in% unique(ped$Father[ped$Generation == 43]),]
ped[ped$Indiv %in% unique(ped$Father[ped$Generation == 52]),]
ped[ped$Indiv %in% unique(ped$Father[ped$Generation == 46]),]
ped[ped$Indiv %in% unique(ped$Father[ped$Generation == 46]),]
mean(ped$gvNormUnres1[ped$Indiv %in% unique(ped$Father[ped$Generation == 46])])
ped[ped$Indiv %in% unique(ped$Mother[ped$Generation == 45]),]
mean(ped$gvNormUnres1[ped$Indiv %in% unique(ped$Mother[ped$Generation == 46])])
mean(ped$gvNormUnres1[ped$Indiv %in% unique(ped$Mother[ped$Generation == 45])])
mean(ped$gvNormUnres1[ped$Indiv %in% unique(ped$Mother[ped$Generation == 47])])
table(ped$Father[ped$Generation == 46])

ped$mst <- NA
for (row in (40*8640):nrow(ped)) {
  ped$mst[row] <- mean(ped$gvNormUnres1[ped$Indiv == ped$Mother[row]], ped$gvNormUnres1[ped$Indiv == ped$Father[row]]) - ped$gvNormUnres1[row]
}

ped$MST <- NA

fathers <-  unique(ped$Father[ped$Generation %in% 40:60])
mothers <-  unique(ped$Mother[ped$Generation %in% 40:60])
pedFather <- ped[ped$Indiv %in% fathers,c("Indiv", "gvNormUnres1")]
pedMother <- ped[ped$Indiv %in% mothers,c("Indiv", "gvNormUnres1")]
colnames(pedFather) <- c("Father", "gvF")
colnames(pedMother) <- c("Mother", "gvM")

ped1 <- ped[ped$Generation %in% 40:60,]
nrow(ped1)
ped1 <- merge(ped1, pedFather, by="Father", all.x=TRUE)
nrow(ped1)
ped1 <- merge(ped1, pedMother, by="Mother", all.x=TRUE)
nrow(ped1)

ped1$MST <- NA
#ped1$PA <- 
ped1$gvF <- as.numeric(ped1$gvF)
ped1$gvM <- as.numeric(ped1$gvM)
ped1$MST <- ((ped1$gvF +  ped1$gvM) / 2) - ped1$gvNormUnres1
tail(ped1)

ped$TGVFather <- ped$gvNormUnres1[ped$Indiv == ]


MeanParents <- data.frame(Generation=NA, Mean=NA, Parent=NA)
for (gen in 40:60) {
  MeanParents <- rbind(MeanParents, c(gen, mean(ped1$gvNormUnres1[ped1$Indiv %in% unique(ped1$Mother[ped1$Generation == gen])]), "Mother"))
  MeanParents <- rbind(MeanParents, c(gen, mean(ped1$gvNormUnres1[ped1$Indiv %in% unique(ped1$Father[ped1$Generation == gen])]), "Father"))
  MeanParents <- rbind(MeanParents, c(gen, mean(ped1$gvNormUnres1[ped1$Generation == gen]), "Population"))
  MeanParents <- rbind(MeanParents, c(gen, mean(ped1$MST[ped1$Generation == gen]), "MST"))
  
}

write.csv(MeanParents,  "~/MeanParents.csv", quote=FALSE, row.names=FALSE)

head(ped1[which(ped1$MST < 0),])

MeanParents$Mean <- as.numeric(MeanParents$Mean)
ggplot(data=MeanParents, aes(x=Generation, y=Mean, group=Parent, fill=Parent)) + geom_bar(stat="identity", position="dodge")

hist(ped$gvNormUnres1[ped$Generation == 45])
hist(ped$gvNormUnres1[ped$Generation == 46])

p45 <- qplot(gvNormUnres1, data=ped[ped$Generation == 45,], geom="histogram", bins=50) + scale_x_continuous(breaks=seq(from=-3, to=10, by=1)) + xlim(c(-3, 10)) + ggtitle("Generation 45")
p46 <- qplot(gvNormUnres1, data=ped[ped$Generation == 46,], geom="histogram", bins=50) + scale_x_continuous(breaks=seq(from=-3, to=10, by=1)) + xlim(c(-3, 10))+ ggtitle("Generation 46")

library(Rmisc)
multiplot(p45, p46)

sum(ped$Mother[ped$Generation==46] == 0)
sum(is.na(ped$Mother[ped$Generation==46]))
