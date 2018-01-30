import pandas as pd

pedigree = pd.read_table("~/Documents/PhD/Simulaton/INBREEDING.txt", sep="\s+", header=None)
pedigree.loc[:, "Rep"] = pedigree[2].str.split("_").str[0]
pedigree.loc[:, "Scenario"] = pedigree[2].str.split("_").str[1]
pedigree.columns = ["Indiv", "F", "Group", "Rep", "Scenario"]
ped = pedigree.loc[pedigree.Indiv != "dumMother"]
ped1 = pedigree.loc[pedigree.Indiv != "dumFather"]
max(pedigree.Indiv)

generation = pd.read_table("~/Documents/PhD/Simulaton/GenerationCLASS.txt", sep=" ")
generation = generation[["Indiv", "Generation"]]

len(set(generation.Indiv) & set(pedigree[0]))