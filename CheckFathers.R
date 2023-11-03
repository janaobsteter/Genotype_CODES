ped <- read.table("SimulatedData/PedigreeAndGeneticValues_cat.txt", header=TRUE, sep=" ")
ps <- read.csv("PopulationSplit.txt")
colnames(ps) <- c("Group", "Indiv")
ped <- merge(ped, ps, by="Indiv")

print("nr HOME father")
print(nrow(ped[ped$Father %in% ps$Indiv[ps$Group == "home"] & ped$cat == "nr" & ped$Group == "home",]))
print("nr IMPORT father")
print(nrow(ped[ped$Father %in% ps$Indiv[ps$Group == "import"] & ped$cat == "nr" & ped$Group == "home",]))
print("potomci HOME father")
print(nrow(ped[ped$Father %in% ps$Indiv[ps$Group == "home"] & ped$cat == "potomciNP" & ped$Group == "home",]))
print("potomci IMPORT father")
print(nrow(ped[ped$Father %in% ps$Indiv[ps$Group == "import"] & ped$cat == "potomciNP" & ped$Group == "home",]))

pedNR <- ped[ped$cat == "nr" & ped$Group == "home",]
pedNR <- ped[ped$cat == "nr" & ped$Group == "home",]
f <- ped[ped$Indiv %in% pedNR$Father,]
print("Fathers home nr")
print(table(f$cat, f$Group))

