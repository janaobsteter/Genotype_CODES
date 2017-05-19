import os
import pandas as pd



def IndCat(ind):
    os.chdir('/home/jana/bin/AlphaSim1.05Linux/')
    indcat = []
    Cat = sorted([ i for i in os.listdir('/home/jana/bin/AlphaSim1.05Linux/') if i.startswith('Cat')])
    for i in Cat:
        catDF = pd.read_csv(i)
        for cat in catDF.columns:
            if float(ind) in set(catDF[cat].dropna()):
                indcat.append(cat) 
    return indcat

    