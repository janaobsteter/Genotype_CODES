import pandas as pd
import os

os.chdir("/home/jana/Documents/PhD/CompBio/TestingGBLUP/")

blupSol = pd.read_csv('renumbered_Solutions', header=None,
                      sep='\s+', names=['renID', 'ID', 'Solution'])
AlphaPed = pd.read_table("PedigreeAndGeneticValues_cat.txt", sep=" ")
AlphaSelPed = AlphaPed.loc[:, ['Generation', 'Indiv', 'Father', 'Mother','cat', 'gvNormUnres1']]
AlphaSelPed.loc[:, 'EBV'] = blupSol.Solution
AlphaSelPed = AlphaSelPed.loc[AlphaSelPed.Generation.isin([40])]
import numpy
print numpy.corrcoef(AlphaSelPed.EBV, AlphaSelPed.gvNormUnres1)
AlphaSelPed.to_csv('GenPed_EBV.txt', index=None)
