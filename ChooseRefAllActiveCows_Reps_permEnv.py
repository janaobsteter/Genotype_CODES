# -*- coding: utf-8 -*-
import pandas as pd
import os
import sys

rep = sys.argv[1]
repEnd = sys.argv[2]

print(range(0, int(rep) + 1))

for rep in range(rep, int(repEnd)+1):
	os.chdir('/home/v1jobste/JanaO/FillInBurnIn' + str(rep) + '_permEnv/')

	ped = pd.read_table('./SimulatedData/PedigreeAndGeneticValues_cat.txt', sep='\s+')
	ped.loc[ped.cat == "k", 'Indiv'].to_csv('/home/v1jobste/JanaO/FillInBurnIn' + str(rep) + '_permEnv/IndForGeno_10K.txt', index=None, header=None, sep='\n')


