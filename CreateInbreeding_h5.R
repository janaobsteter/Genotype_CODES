library("rhdf5")
library(readr)

h5createFile("Inbreeding.h5")


for (strategy in c("SU55", "SU51", "SU15")) {
  h5createGroup("Inbreeding.h5", strategy)
  for (scenario in c("Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen")) {
    for (rep in 0:19) {
      inb <- as.matrix(read_table(paste0(strategy, "/Inbreeding_", scenario, rep, ".txt"),  col_names=FALSE,  col_types = cols(X1 = col_character(), .default = col_double())))
      print(head(inb))
      colnames(inb) <- c("Indiv", "F")
      #h5ls("Inbreeding.h5")
      h5write(inb, file="Inbreeding.h5", name=paste0(strategy, "/", scenario, rep))

      fid <- H5Fopen("Inbreeding.h5")
      id <- H5Dopen(fid, paste0(strategy, "/", scenario, rep))
      h5writeAttribute(id, name="Strategy", attr=strategy)
      h5writeAttribute(id, name="Scenario", attr=scenario)
      h5writeAttribute(id, name="Rep", attr=rep)
      H5Dclose(id)
      H5Fclose(fid)
    }
  }
  fid <- H5Fopen("Inbreeding.h5")
  gid <- H5Gopen(fid, paste0(strategy, "/"))
  h5writeAttribute(gid, name="Pedigree Inbreeding, strategy: ", attr=strategy)
  H5Gclose(gid)
}

fid <- H5Fopen("Inbreeding.h5")
h5writeAttribute(fid, name="Inbreeding type", attr="Pedigree")

H5Fclose(fid)
H5close()

