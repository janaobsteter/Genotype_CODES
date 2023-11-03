from scipy.stats.stats import pearsonr   
import pandas as pd
import numpy as np
from collections import defaultdict
import math
import os
import pandas as pd
import sys

strategy = sys.argv[1]

WorkingDir = os.getcwd()


for scenario in ["Class", "GenSLO", "OtherCowsGen", "BmGen", "Gen"]:
	for rep in range(0, 21):
		print(scenario, rep)
		os.chdir(WorkingDir + '/' +  scenario + str(rep))
		#os.chdir(AlphaSimDir + 'REAL20GenSel_GenSplosnaPop/')
		TBVs_all = pd.read_table('SimulatedData/PedigreeAndGeneticValues.txt', sep="\s+").loc[:, ['Generation', 'Indiv', 'Father', 'Mother', 'gvNormUnres1']]
		accuraciesEBV = pd.DataFrame(columns = ['Strategy', 'Scenario', 'Rep', 'Cycle', 'Cat', 'mseEBV', 'msePA', 'meanEBV', 'meanTBV', 'meanPA'])
		ped = pd.read_csv("SimulatedData/PedigreeAndGeneticValues_cat.txt", sep=" ")

		files =	os.listdir(os.getcwd())
		sol = [i for i in files if i.startswith("renumbered_Solutions_")]
		number = [int(i.strip("renumbered_Solutions_")) for i in sol]

		for i in range(41,max(number)):
		    TBVs = TBVs_all.loc[TBVs_all.Generation.isin(range(i-11,i+1))]
		    accuraciesEBV_temp = pd.DataFrame(columns = ['Rep', 'Cycle', 'Cat', 'mseEBV', 'msePA', 'meanEBV', 'meanTBV', 'meanPA'])
		    solution = pd.read_table("renumbered_Solutions_" + str(i), header=None, sep="\s+", names = ['renID', 'Indiv', 'EBV', 'Cycle']).loc[:, ['Indiv', 'EBV']]
		    mseDF = pd.merge(TBVs, solution, on='Indiv')
		    mseDF.loc[:, 'PA'] = ""

		    catDFEx = 'Categories_gen' + str(i) + 'DF.csv' 
		    categories = defaultdict(list)
		    catDF = pd.read_csv(catDFEx)
		    accuraciesEBV_temp.loc[:, 'Rep'] = [rep] * len(catDF.columns)	
		    accuraciesEBV_temp.loc[:, 'Cycle'] = [i] * len(catDF.columns)
                    accuraciesEBV_temp.loc[:, 'Strategy'] = [strategy] * len(catDF.columns)
                    accuraciesEBV_temp.loc[:, 'Scenario'] = [scenario] * len(catDF.columns)


		    accuraciesEBV_temp.loc[:,'Cat'] = list(catDF.columns)
		    accuraciesEBV_temp.loc[:,'YOB'] = [""] * len(catDF.columns) 
		    for cat in catDF.columns:
			values = [int(ind) for ind in catDF[cat] if not math.isnan(ind)]
			categories[cat] = values
			OneCatDF = mseDF[mseDF.Indiv.isin(categories[cat])]
			try:
				accuraciesEBV_temp.loc[(accuraciesEBV_temp.Cycle == i) & (accuraciesEBV_temp.Cat == cat), 'mseEBV'] = np.mean((OneCatDF.gvNormUnres1[OneCatDF.Generation == max(OneCatDF.Generation)] - OneCatDF.EBV[OneCatDF.Generation == max(OneCatDF.Generation)])**2)
				accuraciesEBV_temp.loc[(accuraciesEBV_temp.Cycle == i) & (accuraciesEBV_temp.Cat == cat), 'meanEBV'] = np.mean(OneCatDF.EBV[OneCatDF.Generation == max(OneCatDF.Generation)])
				accuraciesEBV_temp.loc[(accuraciesEBV_temp.Cycle == i) & (accuraciesEBV_temp.Cat == cat), 'meanTBV'] = np.mean(OneCatDF.gvNormUnres1[OneCatDF.Generation == max(OneCatDF.Generation)])
				accuraciesEBV_temp.loc[(accuraciesEBV_temp.Cycle == i) & (accuraciesEBV_temp.Cat == cat), 'YOB'] = max(OneCatDF.Generation)
			except: 
				pass



		    accuraciesEBV = pd.concat([accuraciesEBV, accuraciesEBV_temp])
		    
		    
	accuraciesEBV.loc[:, 'Cycle'] = accuraciesEBV.Cycle.astype(int)
	accuraciesEBV.to_csv(WorkingDir + '/AccuraciesEBVPA_' + scenario + 'MSE.csv', index=None, sep=",")

	
	print 'Create Data Frame: AccuraciesEBVPA' + scenario + '.csv'

