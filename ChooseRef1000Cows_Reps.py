# -*- coding: utf-8 -*-
import random
import pandas as pd
import os
from itertools import chain
import sys

rep = sys.argv[1]

print(range(0, int(rep) + 1))

for rep in range(0, int(rep)+1):
	os.chdir('/home/v1jobste/JanaO/FillInBurnIn' + str(rep) + '/')

	numAnimalperGen = 250
	gen = 6
	ks = pd.DataFrame({'Indiv': sorted(list(set(chain.from_iterable([random.sample((list(pd.read_table('Categories_gen' + str(i) + 'DF.csv', sep=",")['k'].dropna().astype(int))), numAnimalperGen) for i in range(41-gen, 41)]))))})
	ped = pd.read_table('./SimulatedData/PedigreeAndGeneticValues_cat.txt', sep='\s+')
	ks = pd.merge(ks, ped, on='Indiv')
	ks.loc[ks.Generation.isin(range(max(ks.Generation)-6, max(ks.Generation)+1)), 'Indiv'].to_csv('/home/v1jobste/JanaO/FillInBurnIn' + str(rep) + '/IndForGeno_1000.txt', index=None, header=None, sep='\n')

#ped[ped.Indiv.isin(ks.loc[ks.Generation.isin(range(max(ks.Generation)-6, max(ks.Generation)+1)), 'Indiv'])].Generation.value_counts()

