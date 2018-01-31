import pandas as pd

pedigree = pd.read_table("~/Documents/PhD/Simulaton/INBREEDING_1PB.txt", sep="\s+", header=None)
pedigree.columns = ["Indiv", "F", "Group"]
groups = set(pedigree.Group)

generation = pd.read_table("~/Documents/PhD/Simulaton/GenerationCLASS.txt", sep=" ")
generation = generation[["Indiv", "Generation"]]
generation= generation.loc[generation.Indiv != "dumMother"]
generation= generation.loc[generation.Indiv != "dumFather"]
generation = generation.loc[generation.Indiv != "Indiv"]
generation.Indiv = pd.to_numeric(generation.Indiv)

for group in groups:
    ped = pedigree.loc[pedigree.Group == group]
    ped.loc[:, "Rep"] = ped.Group.str.split("_").str[0]
    ped.loc[:, "Scenario"] = ped.Group.str.split("_").str[1]
    ped.columns = ["Indiv", "F", "Group", "Rep", "Scenario"]

    ped = ped.loc[ped.Indiv != "dumMother"]
    ped = ped.loc[ped.Indiv != "dumFather"]
    ped = ped.loc[ped.Indiv != "Indiv"]

    ped.Indiv = pd.to_numeric(ped.Indiv)
    max(ped.Indiv)
    PED = pd.merge(ped, generation, on="Indiv", how="left")

    
    len(set(generation.Indiv) & set(pedigree[0]))