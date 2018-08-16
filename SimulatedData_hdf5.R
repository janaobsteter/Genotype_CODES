library("rhdf5")
library(readr)

setwd("/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Gen1/SimulatedData/AllIndividualsSnpChips")
gen <- as.matrix(read_table("Last20Gen_small.txt",  col_names=FALSE,  col_types = cols(X1 = col_character(), .default = col_integer())))
colnames(gen) <- c("ID", paste0("M", 1:(ncol(gen)-1)))

h5file("SimulationData.h5")
h5createFile("SimulationData.h5")
h5createGroup("SimulationData.h5", "Inbreeding")
h5createGroup("SimulationData.h5", "Genotypes")
h5ls("SimulationData.h5")

h5write(gen, file="SimulationData.h5", name="Genotypes/Class1")
h5ls("SimulationData.h5")

fid <- H5Fopen("SimulationData.h5")
id <- H5Dopen(fid, "Genotypes/Class1")
h5writeAttribute(id, name="Strategy", attr="SU55")
h5writeAttribute(id, name="Scenario", attr="Class")
h5writeAttribute(id, name="Rep", attr="0")

gid <- H5Gopen(fid, "Genotypes/")
h5writeAttribute(gid, name="chip number", attr="1")

H5Dclose(id)
H5Gclose(gid)
H5Fclose(fid)

h5readAttributes(file="SimulationData.h5", name="Genotypes")
h5readAttributes(file="SimulationData.h5", name="Genotypes/Class1")

testSubset <- as.numeric(h5read(file="SimulationData.h5", name="Genotypes/Class1", index=list(NULL,2)))
hist(testSubset)
