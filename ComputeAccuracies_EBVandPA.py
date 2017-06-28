from scipy.stats.stats import pearsonr   
import pandas as pd
import numpy as np
from collections import defaultdict
import math

os.chdir(AlphaSimDir + 'REAL20GenSel_Gen/')
TBVs_all = pd.read_table('SimulatedData/PedigreeAndGeneticValues.txt', sep="\s+").loc[:, ['Generation', 'Indiv', 'Father', 'Mother', 'gvNormUnres1']]
accuraciesEBV = pd.DataFrame(columns = ['Cycle', 'Cat', 'corEBV', 'corPA', 'meanEBV', 'meanTBV', 'meanPA'])

for i in range(21, 61):
    TBVs = TBVs_all.loc[TBVs_all.Generation.isin(range(i-11,i+1))]
    accuraciesEBV_temp = pd.DataFrame(columns = ['Cycle', 'Cat', 'corEBV', 'corPA', 'meanEBV', 'meanTBV', 'meanPA'])
    solution = pd.read_table("renumbered_Solutions_" + str(i), header=None, sep=" ", names = ['renID', 'Indiv', 'EBV']).loc[:, ['Indiv', 'EBV']]
    corDF = pd.merge(TBVs, solution, on='Indiv')
    corDF.loc[:, 'PA'] = ""

    catDFEx = 'Categories_gen' + str(i) + 'DF.csv' 
    categories = defaultdict(list)
    catDF = pd.read_csv(catDFEx)
    accuraciesEBV_temp.loc[:, 'Cycle'] = [i] * len(catDF.columns)
    accuraciesEBV_temp.loc[:,'Cat'] = list(catDF.columns)
    for cat in catDF.columns:
        values = [int(ind) for ind in catDF[cat] if not math.isnan(ind)]
        categories[cat] = values
        accuraciesEBV_temp.loc[(accuraciesEBV_temp.Cycle == i) & (accuraciesEBV_temp.Cat == cat), 'corEBV'] = pearsonr(corDF.gvNormUnres1[corDF.Indiv.isin(categories[cat])], corDF.EBV[corDF.Indiv.isin(categories[cat])])[0]
        accuraciesEBV_temp.loc[(accuraciesEBV_temp.Cycle == i) & (accuraciesEBV_temp.Cat == cat), 'meanEBV'] = np.mean(corDF.EBV[corDF.Indiv.isin(categories[cat])])
        accuraciesEBV_temp.loc[(accuraciesEBV_temp.Cycle == i) & (accuraciesEBV_temp.Cat == cat), 'meanTBV'] = np.mean(corDF.gvNormUnres1[corDF.Indiv.isin(categories[cat])])



    for row in corDF.loc[corDF.Indiv.isin(categories['potomciNP'])].index:
        if (corDF.Father[row] != 0 and corDF.Mother[row] != 0):
            corDF.loc[row, 'PA'] = ((float(corDF.EBV[corDF.Indiv == corDF.Father[row]]) + float(corDF.EBV[corDF.Indiv == corDF.Mother[row]])) / 2) 
    accuraciesEBV_temp.loc[(accuraciesEBV_temp.Cycle == i) & (accuraciesEBV_temp.Cat == 'potomciNP'), 'corPA'] = pearsonr(corDF.gvNormUnres1[corDF.Indiv.isin(categories['potomciNP'])], corDF.EBV[corDF.Indiv.isin(categories['potomciNP'])])[0]
    accuraciesEBV_temp.loc[(accuraciesEBV_temp.Cycle == i) & (accuraciesEBV_temp.Cat == 'potomciNP'), 'meanPA'] = np.mean(corDF.PA[corDF.Indiv.isin(categories['potomciNP'])])

    accuraciesEBV = pd.concat([accuraciesEBV, accuraciesEBV_temp])
    
    
accuraciesEBV.loc[:, 'Cycle'] = accuraciesEBV.Cycle.astype(int)
accuraciesEBV.to_csv('AccuraciesEBVPA.csv', index=None, sep=",")


print 'Create Data Frame: AccuraciesEBVPA.csv'