#general check of the pedigree
ped = pd.read_csv(AlphaSimDir +"/SimulatedData/PedigreeAndGeneticValues_cat.txt", sep=" ")
#check pb transition
pb = ped.loc[ped.cat=='pb']

set(pb.Generation)
cat40 = pd.read_table(AlphaSimDir + "Categories_gen40DF.csv", sep=",").pb.dropna().astype(int)
cat41 = pd.read_table(AlphaSimDir + "Categories_gen41DF.csv", sep=",").pb.dropna().astype(int)
cat42 = pd.read_table(AlphaSimDir + "Categories_gen42DF.csv", sep=",").pb.dropna().astype(int)