import random
import pandas as pd

chromosome =   [0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0]


#chromosome = random.sample(chromosome, 100)

ped = pd.read_csv("/home/jana/Documents/PhD/CompBio/PedigreeAndGeneticValues_Herds.txt", sep=",")
pedO = pd.read_csv("/home/jana/Documents/PhD/CompBio/TestingGBLUP/PedigreeAndGeneticValues_cat.txt", sep="\s+")
genK = [herd for (herd, gen) in zip(sorted(list(set(ped.Herd))), chromosome) if gen ==1]  
pd.DataFrame({"ID": list(ped.loc[ped.Herd.isin(genK), 'Indiv']) + list(pedO.loc[pedO.Generation.isin([40]),'Indiv']) }).to_csv('/home/jana/Documents/PhD/CompBio/TestingGBLUP/IndForGeno.txt', index=None, header=None)
print(sum(chromosome))



