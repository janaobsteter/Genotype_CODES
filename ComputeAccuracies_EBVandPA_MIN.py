from scipy.stats.stats import pearsonr   
import pandas as pd
import numpy as np
from collections import defaultdict
import math
import os
import pandas as pd

#os.chdir(AlphaSimDir)
os.chdir(AlphaSimDir + 'REAL20GenSel_Class/')
TBVs_all = pd.read_table('SimulatedData/PedigreeAndGeneticValues.txt', sep="\s+").loc[:, ['Generation', 'Indiv', 'Father', 'Mother', 'gvNormUnres1']]
accuraciesEBV = pd.DataFrame(columns = ['Cycle', 'Cat', 'corEBV', 'corPA', 'meanEBV', 'meanTBV', 'meanPA'])
ped = pd.read_csv(AlphaSimDir +"/SimulatedData/PedigreeAndGeneticValues_cat.txt", sep=" ")

for i in range(41, 61):
    TBVs = TBVs_all.loc[TBVs_all.Generation.isin(range(i-11,i+1))]
    accuraciesEBV_temp = pd.DataFrame(columns = ['Cycle', 'Cat', 'corEBV', 'corPA', 'meanEBV', 'meanTBV', 'meanPA'])
    solution = pd.read_table("renumbered_Solutions_" + str(i), header=None, sep="\s+", names = ['renID', 'Indiv', 'EBV', 'Cycle']).loc[:, ['Indiv', 'EBV']]
    corDF = pd.merge(TBVs, solution, on='Indiv')
    corDF.loc[:, 'PA'] = ""

    catDFEx = 'Categories_gen' + str(i) + 'DF.csv' 
    categories = defaultdict(list)
    catDF = pd.read_csv(catDFEx)
    accuraciesEBV_temp.loc[:, 'Cycle'] = [i] * len(catDF.columns)
    accuraciesEBV_temp.loc[:,'Cat'] = list(catDF.columns)
    accuraciesEBV_temp.loc[:,'YOB'] = [""] * len(catDF.columns) 
    for cat in catDF.columns:
        values = [int(ind) for ind in catDF[cat] if not math.isnan(ind)]
        categories[cat] = values
        OneCatDF = corDF[corDF.Indiv.isin(categories[cat])]
        accuraciesEBV_temp.loc[(accuraciesEBV_temp.Cycle == i) & (accuraciesEBV_temp.Cat == cat), 'corEBV'] = pearsonr(OneCatDF.gvNormUnres1[OneCatDF.Generation == max(OneCatDF.Generation)], OneCatDF.EBV[OneCatDF.Generation == max(OneCatDF.Generation)])[0]
        accuraciesEBV_temp.loc[(accuraciesEBV_temp.Cycle == i) & (accuraciesEBV_temp.Cat == cat), 'meanEBV'] = np.mean(OneCatDF.EBV[OneCatDF.Generation == max(OneCatDF.Generation)])
        accuraciesEBV_temp.loc[(accuraciesEBV_temp.Cycle == i) & (accuraciesEBV_temp.Cat == cat), 'meanTBV'] = np.mean(OneCatDF.gvNormUnres1[OneCatDF.Generation == max(OneCatDF.Generation)])
        accuraciesEBV_temp.loc[(accuraciesEBV_temp.Cycle == i) & (accuraciesEBV_temp.Cat == cat), 'YOB'] = max(OneCatDF.Generation)


    for row in corDF.loc[corDF.Indiv.isin(categories['potomciNP'])].index:
        if (corDF.Father[row] != 0 and corDF.Mother[row] != 0):
            corDF.loc[row, 'PA'] = ((float(corDF.EBV[corDF.Indiv == corDF.Father[row]]) + float(corDF.EBV[corDF.Indiv == corDF.Mother[row]])) / 2) 
    accuraciesEBV_temp.loc[(accuraciesEBV_temp.Cycle == i) & (accuraciesEBV_temp.Cat == 'genTest'), 'corPA'] = pearsonr(corDF.gvNormUnres1[corDF.Indiv.isin(categories['genTest'])], corDF.PA[corDF.Indiv.isin(categories['genTest'])])[0]
    accuraciesEBV_temp.loc[(accuraciesEBV_temp.Cycle == i) & (accuraciesEBV_temp.Cat == 'genTest'), 'meanPA'] = np.mean(corDF.PA[corDF.Indiv.isin(categories['genTest'])])

    accuraciesEBV = pd.concat([accuraciesEBV, accuraciesEBV_temp])
    
    
accuraciesEBV.loc[:, 'Cycle'] = accuraciesEBV.Cycle.astype(int)
accuraciesEBV.to_csv('AccuraciesEBVPA.csv', index=None, sep=",")


print 'Create Data Frame: AccuraciesEBVPA.csv'

