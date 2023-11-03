# -*- coding: utf-8 -*-
import pandas as pd
import os
import sys

rep = int(sys.argv[1])
repEnd = int(sys.argv[2])

print(range(rep, repEnd))

for rep in range(rep, repEnd):
	os.chdir('/home/v1jobste/JanaO/FillInBurnIn' + str(rep) + '_permEnv/')

	ped = pd.read_table('./SimulatedData/PedigreeAndGeneticValues_cat.txt', sep='\s+')
	ped.loc[ped.cat.isin(["k", "pb"]), 'Indiv'].to_csv('/home/v1jobste/JanaO/FillInBurnIn' + str(rep) + '_permEnv/IndForGeno_10K.txt', index=None, header=None, sep='\n')


