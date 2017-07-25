import os
import pandas as pd

os.chdir(AlphaSimDir)
os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Gen/')
os.chdir("/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenSplosnaPop")
def IndCat(ind):
    indcat = []
    Cat = sorted([ i for i in os.listdir(os.getcwd()) if i.startswith('Cat')])
    for i in Cat:
        catDF = pd.read_csv(i)
        for cat in catDF.columns:
            if float(ind) in set(catDF[cat].dropna()):
                print cat
                indcat.append(cat) 
    return indcat

    
    
potomci = list(pd.read_csv('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Class/Categories_gen56DF.csv')['cak'].dropna().astype(int))
sol = pd.read_table('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Class/renumbered_Solutions_56', header=None, sep=" ", names=['rind', 'Indiv', 'EBV'])
pot = ped.loc[ped.Indiv.isin(potomci)]
a = pd.merge(pot, sol, on='Indiv')
a.loc[a.sex=='M'].sort(a.columns[17], ascending=False)[:30]

cak = list(pd.read_table('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Class/Categories_gen59DF.csv', sep=",").cak.dropna().astype(int))
sol = pd.read_table('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Class/renumbered_Solutions_59', header=None, sep=" ", names=['rind', 'Indiv', 'EBV'])
ped = pd.read_table('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Class/SimulatedData/PedigreeAndGeneticValues.txt', sep='\s+')