
# ---- Requirements ----

library(package = "tidyverse")
library(package = "ggplot2")
library(package = "AlphaSimR") # for sortPed
library(package = "pedigreemm")
library(package = "igraph")

# ---- Read in data ----

dat = read_csv2(file = "podatki_2017-10-10.txt")
head(dat)

colnames(dat) = tolower(colnames(dat))
head(dat)

# ---- Filter the data ----

# Zival
dat$rod = as.factor(dat$rod)

# Datum rojstva
dat$datum_rojstva = as.Date(dat$datum_rojstva, format = "%d.%m.%Y")
sel = is.na(dat$datum_rojstva)
sum(sel)
tmp = dat[sel, ]
file.remove("Neznan_datum_rojstva.csv")
if (nrow(tmp) > 1) {
  View(tmp)
  write.csv2(x = tmp, file = "Neznan_datum_rojstva.csv", row.names = FALSE)
  dat = dat[!sel, ]
}

# Pasma
dat$pasma = as.factor(dat$pasma)
table(dat$pasma)
sel = is.na(dat$pasma)
sum(sel)
tmp = dat[sel, ]
file.remove("Neznana_pasma.csv")
if (nrow(tmp) > 1) {
  View(tmp)
  write.csv2(x = tmp, file = "Neznana_pasma.csv", row.names = FALSE)
  dat$pasma[sel] = 1
}

# Velikost gnezda
dat$velikost_gnezda = as.factor(dat$velikost_gnezda)
table(dat$velikost_gnezda)
sel = is.na(dat$velikost_gnezda)
sum(sel)
tmp = dat[sel, ]
file.remove("Neznana_velikost_gnezda.csv")
if (nrow(tmp) > 1) {
  View(tmp)
  write.csv2(x = tmp, file = "Neznana_velikost_gnezda.csv", row.names = FALSE)
  dat$velikost_gnezda[sel] = 1
}

# Zaporedno gnezdo
dat$zap_gnezdo = as.factor(dat$zap_gnezdo)
table(dat$zap_gnezdo)
sel = is.na(dat$zap_gnezdo)
sum(sel)
tmp = dat[sel, ]
file.remove("Neznano_zaporedno_gnezdo.csv")
if (nrow(tmp) > 1) {
  View(tmp)
  write.csv2(x = tmp, file = "Neznano_zaporedno_gnezdo.csv", row.names = FALSE)
  dat$zap_gnezdo[sel] = 3
}

# Izvor
dat$izvor = as.factor(dat$izvor)
tmp = table(dat$izvor)
hist(tmp)
table(tmp)
sel = is.na(dat$izvor)
sum(sel)
tmp = dat[sel, ]
file.remove("Neznan_izvor.csv")
if (nrow(tmp) > 1) {
  View(tmp)
  write.csv2(x = tmp, file = "Neznan_izvor.csv", row.names = FALSE)
  dat$izvor[sel] = 0
}

# Testna postaja
dat$testna_postaja = as.factor(dat$testna_postaja)
table(dat$testna_postaja)
sel = is.na(dat$testna_postaja)
tmp = dat[sel, ]
file.remove("Neznana_testna_postaja.csv")
if (nrow(tmp) > 1) {
  View(tmp)
  write.csv2(x = tmp, file = "Neznana_testna_postaja.csv", row.names = FALSE)
  dat = dat[!sel, ]
}

ggplot(data = dat, aes(starost)) + geom_histogram()

# Premladi
tmp = dat %>%
  filter(starost < 50) %>%
  arrange(rod, starost)
file.remove("Premladi.csv")
if (nrow(tmp) > 1) {
  View(tmp)
  write.csv2(x = tmp, file = "Premladi.csv", row.names = FALSE)
  sel = dat$rod %in% tmp$rod
  dat = dat[!sel, ]
}

# Prestari
tmp = dat %>%
  filter(starost > 400) %>%
  arrange(rod, starost)
file.remove("Prestari.csv")
if (nrow(tmp) > 1) {
  View(tmp)
  write.csv2(x = tmp, file = "Prestari.csv", row.names = FALSE)
  sel = dat$rod %in% tmp$rod
  # View(dat[sel, ] %>% arrange(datum_rojstva, rod, starost))
  # table(format(dat[sel, ]$"datum_rojstva", format = "%Y"))
  dat = dat[!sel, ]
}

ggplot(data = dat, aes(starost)) + geom_histogram()

# Mase
tmp = dat %>%
  filter(is.na(teza) | !is.finite(teza)) %>%
  arrange(rod, starost)
file.remove("Neznana_teza.csv")
if (nrow(tmp) > 1) {
  View(tmp)
  write.csv2(x = tmp, file = "Neznana_teza.csv", row.names = FALSE)
  sel = paste(dat$rod, dat$starost, sep = "-") %in% paste(tmp$rod, tmp$starost, sep = "-")
  dat = dat[!sel, ]
}

ggplot(data = dat, aes(teza)) + geom_histogram()
ggplot(data = dat, aes(x = starost, y = teza)) + geom_point()

# ---- Read in pedigree ----

ped = read.csv("~/Documents/PhD/CompBio/TestingGBLUP/PedigreeAndGeneticValues_cat.txt", sep=" ")
head(ped)
#chose genotipying candidated as the target individuals
dat <- ped[ped$cat %in% c("k", "potomciNP", "nr"),]

colnames(ped) = tolower(colnames(ped))
head(ped)

ped$datum_rojstva = as.Date(ped$datum_rojstva, format = "%d.%m.%Y")

# ---- Build pedigree for the data ----

nGenerations = 5

# Phenotyped individuals - in this case genotyping candidated
ped2 = ped[ped$cat %in% c("k", "potomciNP", "nr"),c(1,2,3,4)]
ped <- ped[,c(1,2,3,4)]
nrow(ped2)
length(unique(dat$rod))
if (nrow(ped2) != length(unique(dat$rod))) {
  stop("TODO: some phenotyped individuals not in pedigree")
}

# Add generations
for (Generation in 2:nGenerations) {
  sel = ped2$Mother %in% ped$Indiv &
    !ped2$Mother %in% ped2$Indiv &
    !is.na(ped2$Mother)
  if (sum(sel) > 0) {
    ped2 = rbind(ped2, ped[ped$Indiv %in% ped2$Mother[sel], ])
  }
  sel = ped2$Father %in% ped$Indiv &
    !ped2$Father %in% ped2$Indiv &
    !is.na(ped2$Father)
  if (sum(sel) > 0) {
    ped2 = rbind(ped2, ped[ped$Indiv %in% ped2$Father[sel], ])
  }
}

# Remaining founders
sel = !ped2$Mother %in% ped2$Indiv &
  !is.na(ped2$Mother)
if (sum(sel) > 0) {
  ped2 = rbind(ped2, data.frame(Generation = NA,
                                Indiv = ped2$Mother[sel],
                                Father = NA,
                                Mother = NA))
}
sel = !ped2$Father %in% ped2$Indiv &
  !is.na(ped2$Father)
if (sum(sel) > 0) {
  ped2 = rbind(ped2, data.frame(Generation = NA,
                                Indiv = ped2$Father[sel],
                                Father = NA,
                                Mother = NA))
}

# ---- Cleanup pedigree errors ----

sel = ped2$Mother %in% ped2$Father & !is.na(ped2$Mother)
sel2 = ped2$Father %in% ped2$Mother & !is.na(ped2$Father)
sum(sel)
sum(sel2)
ped2[sel, ]
ped2[sel2, ]

tmp = table(table(ped2$Indiv))
if (length(tmp) > 1) {
  stop("duplicates!")
}

# save.image("Data.RData")
# load("Data.RData")

# ---- Sort pedigree so that parents preceede progeny and get 1:n individual codes ----

# Make a graph and do topological sort on a graph
pedGraph = data.frame(from = c(ped2$oce, ped2$mati),
                      to   = c(ped2$rod, ped2$rod))
sel = is.na(pedGraph$from)
pedGraph = pedGraph[!sel, ]
ped2$name = ped2$rod
pedGraph = graph_from_data_frame(d = pedGraph, directed = TRUE, vertices = ped2)
ord = topo_sort(graph = pedGraph, mode = "out") # out to sort from founders to descendants

# Merge the sorted vertex names with the pedigree table
ped2 = merge(x = ped2,
             y = data.frame(rod = ord$name, rodI = 1:length(ord)))
head(ped2); tail(ped2)
# Sort by 1:n codes
ped2 = ped2[order(ped2$rodI), ]
head(ped2); tail(ped2)

if (FALSE) {
  ped2 = ped2[order(ped2$datum_rojstva), ]
  head(ped2); tail(ped2)
  tmp = new("Pedigree",
            nInd = nrow(ped2),
            ids = ped2$rod,
            mother = match(x = ped2$mati,  table = ped2$rod, nomatch = 0),
            father = match(x = ped2$oce,   table = ped2$rod, nomatch = 0))
  tmp = sortPed(x = tmp)
  tmp2 = data.frame(id = tmp@ids, idI = 1:tmp@nInd, sireI = tmp@father, damI = tmp@mother)
  head(tmp2); tail(tmp2)
  tmp2[tmp2$idI < tmp2$sireI, ]
  tmp2[tmp2$idI < tmp2$damI, ]
}

# ---- Add 1:n id codes to phenotype data ---

dat = merge(x = dat, y = ped2[, c("rod", "rodI")])

# ---- Recode parents to 1:n ----

ped2$oceI  = match(x = ped2$oce,  table = ped2$rod, nomatch = 0)
ped2$matiI = match(x = ped2$mati, table = ped2$rod, nomatch = 0)
head(ped2); tail(ped2)

sum(is.na(ped2$rodI))
sum(is.na(ped2$oceI))
sum(is.na(ped2$matiI))

summary(ped2$rodI)
summary(ped2$matiI)
summary(ped2$oceI)

tmp = ped2[ped2$rodI < ped2$oceI, ]
if (nrow(tmp) > 0) {
  stop("ordering failed")
}
tmp = ped2[ped2$rodI < ped2$matiI, ]
if (nrow(tmp) > 0) {
  stop("ordering failed")
}

# ---- Build pedigree precision matrix ----

pedX = pedigree(label = ped2$rodI,
                sire  = ped2$oceI,
                dam   = ped2$matiI)

TInv = as(pedX, "sparseMatrix")
DInv = Diagonal(x = 1 / Dmat(pedX))
QPedigree = crossprod(sqrt(DInv) %*% TInv)
