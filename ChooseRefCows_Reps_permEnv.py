# -*- coding: utf-8 -*-
import random
import pandas as pd
import os
from itertools import chain
import sys

#rep = sys.argv[1]
rep = 1

for reference in [6000, 7000, 8000, 9000, 1000, 2000, 3000, 4000, 5000, 10000]:
	for rep in [1]:
		print(rep)
		os.chdir('/home/v1jobste/JanaO/FillInBurnIn' + str(rep) + "/") 
		print(os.getcwd())
		numAnimalperGen = reference / 3
		gen = 6
		ks = pd.DataFrame({'Indiv': sorted(list(set(chain.from_iterable([random.sample((list(pd.read_table('Categories_gen' + str(i) + 'DF.csv', sep=",")['k'].dropna().astype(int))), numAnimalperGen) for i in range(41-gen, 41)]))))})
		ped = pd.read_table('./SimulatedData/PedigreeAndGeneticValues_cat.txt', sep='\s+')
		ks = pd.concat([ks, pd.DataFrame({"Indiv": list(ped.Indiv[ped.cat == "potomciNP"])})])
		ks = pd.merge(ks, ped, on='Indiv')
		ks.loc[ks.Generation.isin(range(max(ks.Generation)-6, max(ks.Generation)+1)), 'Indiv'].to_csv('/home/v1jobste/JanaO/FillInBurnIn' + str(rep) + '/IndForGeno_' + str(reference /1000) + 'K.txt', index=None, header=None, sep='\n')


