# -*- coding: utf-8 -*-
import random
import pandas as pd
import os
from itertools import chain

#os.chdir('/home/jana/bin/AlphaSim1.05Linux/REALFillIn20BurnIn20/')
numAnimalperGen = 2500
gen = 6
ks = pd.DataFrame({'Indiv': sorted(list(set(chain.from_iterable([random.sample((list(pd.read_table('Categories_gen' + str(i) + 'DF.csv', sep=",")['k'].dropna().astype(int))), numAnimalperGen) for i in range(41-gen, 41)]))))})
ped = pd.read_table('./SimulatedData/PedigreeAndGeneticValues_cat.txt', sep='\s+')
ks = pd.merge(ks, ped, on='Indiv')
#zgoraj ti je lahko pobrali tudi starejše krave, ker 5 generacij nazaj so bile še let let nazaj
ks.loc[ks.Generation.isin(range(max(ks.Generation)-6, max(ks.Generation)+1)), 'Indiv'].to_csv('IndForGeno.txt', index=None, header=None, sep='\n')

#ped[ped.Indiv.isin(ks.loc[ks.Generation.isin(range(max(ks.Generation)-6, max(ks.Generation)+1)), 'Indiv'])].Generation.value_counts()
