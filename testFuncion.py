import os
import pandas as pd

indGeno = pd.read_csv("/home/jana/IndForGeno_new.txt", header=None, names=["Indiv"])
ped = pd.read_table('/home/jana/PedigreeAndGeneticValues_cat.txt', sep='\s+')

indCatSex = pd.merge(indGeno, ped[['Indiv', 'sex', 'cat']], on = 'Indiv')
indCatSex.Indiv[indCatSex.sex == sex]


def obtainNewGenotypedInd_sex(sex):
    indGeno = pd.read_csv(AlphaSimDir + "IndForGeno_new.txt", names=['Indiv'], header=None)
    ped = pd.read_csv(AlphaSimDir + '/PedigreeAndGeneticValues_cat.txt', sep='\s+')
    indCatSex = pd.merge(indGeno, ped[['Indiv', 'sex', 'cat']], on='Indiv')
    return list(indCatSex.Indiv[indCatSex.sex == sex])


