# -*- coding: utf-8 -*-
from selection10 import TBVPed
import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt
from pylab import legend
import os



TBV  = TBVPed(os.getcwd() + '/')
gens, means, vars = TBV.genTrend('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Gen/SimulatedData/PedigreeAndGeneticValues.txt', 41, 60)
TBV  = TBVPed(os.getcwd() + '/')
gens1, means1, vars1 = TBV.genTrend('/home/jana/compareSel/PedigreeAndGeneticValues.txt', 41, 60)



TBV  = TBVPed(os.getcwd() + '/')
gensA, meansA, varsA = TBV.genTrend('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Gen/SimulatedData/PedigreeAndGeneticValues.txt', 41, 6)


TBV  = TBVPed(os.getcwd() + '/')
gensG, meansG, varsG = TBV.genTrend('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenSLO/SimulatedData/PedigreeAndGeneticValues.txt', 41, 61)
#plt.errorbar(x = TBV.gens, y = TBV.means, yerr = TBV.vars)

TBV  = TBVPed(os.getcwd() + '/')
gensS, meansS, varsS = TBV.genTrend('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenSplosnaPop/SimulatedData/PedigreeAndGeneticValues.txt', 41, 61)


TBV  = TBVPed(os.getcwd() + '/')
gensD, meansD, varsD = TBV.genTrend('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenSLO_BmGen/SimulatedData/PedigreeAndGeneticValues.txt', 41, 61)
#plt.errorbar(x = TBV.gens, y = TBV.means, yerr = TBV.vars)
TBV  = TBVPed(os.getcwd() + '/')
gens2, means2, vars2 = TBV.genTrend('/home/jana/bin/AlphaSim1.05Linux//SimulatedData/PedigreeAndGeneticValues.txt', 41, 61)
#plt.errorbar(x = TBV.gens, y = TBV.means, yerr = TBV.vars)

TBV  = TBVPed(os.getcwd() + '/')
gensE, meansE, varsE = TBV.genTrend('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Gen/SimulatedData/PedigreeAndGeneticValues.txt', 41, 61)
#plt.errorbar(x = TBV.gens, y = TBV.means, yerr = TBV.vars)


Means = pd.DataFrame()
Means.loc[:,'Gen'] = [int(x) for x in gens]
Means.loc[:,'Mean'] = list(means)
Means.loc[:,'Var'] =  list(vars)
Means.loc[:,'Stand'] = (Means.Mean - Means.Mean[0]) / sqrt(Means.Var[0])



Means.loc[:,'Gen1'] = [int(x) for x in gens1]
Means.loc[:,'Mean1'] = list(means1)
Means.loc[:,'Var1'] =  list(vars1)
Means.loc[:,'Stand1'] = (Means.Mean1 - Means.Mean1[0]) / sqrt(Means.Var1[0])


Means.loc[:,'GenA'] = [int(x) for x in gensA]
Means.loc[:,'MeanA'] = list(meansA)
Means.loc[:,'VarA'] =  list(varsA)
Means.loc[:,'StandA'] = (Means.MeanA - Means.MeanA[0]) / sqrt(Means.VarA[0])

Means.loc[:,'GenG'] = [int(x) for x in gensG]
Means.loc[:,'MeanG'] = list(meansG)
Means.loc[:,'VarG'] =  list(varsG)
Means.loc[:,'StandG'] = (Means.MeanG - Means.MeanG[0]) / sqrt(Means.VarG[0])


Means.loc[:,'GenS'] = [int(x) for x in gensS]
Means.loc[:,'MeanS'] = list(meansS)
Means.loc[:,'VarS'] =  list(varsS)
Means.loc[:,'StandS'] = (Means.MeanS - Means.MeanS[0]) / sqrt(Means.VarS[0])

Means.loc[:,'GenD'] = [int(x) for x in gensD]
Means.loc[:,'MeanD'] = list(meansD)
Means.loc[:,'VarD'] =  list(varsD)
Means.loc[:,'StandD'] = (Means.MeanD - Means.MeanD[0]) / sqrt(Means.VarD[0])





plt.plot( Means.Gen, Means.Stand,  label = 'Conventional')
plt.plot( Means.Gen1, Means.Stand1,  label = 'genEddie')
plt.plot( Means.GenG, Means.StandG,  label = 'Conventional_SLO')
plt.plot( Means.GenS, Means.StandS,  label = 'GenBulls on Other Cows')
plt.plot( Means.GenD, Means.StandD,  label = 'GenBulls on Bull Dams')
plt.plot( Means.GenA, Means.StandA,  label = 'GenBulls on All Cows')
plt.plot( Means.GenE, Means.StandE,  label = 'GenBulls on All Cows1')
plt.xticks(Means.Gen)

plt.xlabel('Selected Generation')
plt.ylabel('Mean Generation TBV')
legend(loc='upper left')




#ZBRISALA SI EBV ZA KLASIČNO!!! PREBERI IZ RENUMBERED_SOLUTIONS
#EBVs = pd.read_table('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Class/renumbered_Solutions_60', header=None, sep=" ", names = ['rID', 'ID', 'EBV'])
#EBVs.loc[:, 'Generation'] = list(chain.from_iterable([[i] * 8640 for i in range(1, 61)]))
EBVs.to_csv('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Class/EBVs.txt', index=None, sep=" ")
gensEBV = list(set(EBVs.Generation))
meansEBV = EBVs.EBV.groupby(EBVs.Generation).aggregate(np.mean)

EBV = TBVPed(AlphaSimDir)
gensE, meansE, varsE = EBV.genTrend_EBV('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Class/EBVs.txt', 1, 60)

plt.plot(gens, meansEBV, label='EBV')
plt.xlabel('Selected Generation')
plt.ylabel('Mean Generation TBV')
legend(loc='upper left')


plt.plot(gens, means, label='Conventional')
plt.plot(gens1, means1, label='Conventional_rep1')
plt.plot(gensG, meansG, label='Conventional_SLO')
plt.plot(gensS, meansS, label='GenBulls on Other Cows')
plt.plot(gensD, meansD, label='GenBulls on Bull Dams')
plt.plot(gensE, meansE, label='GenBulls on all cows')
plt.plot(gens2, means2, label='GenBulls on Bull Dams2')

#plt.plot(gensA, meansA, label='GenBulls on all cows BitWrong')
plt.xlabel('Selected Generation')
plt.ylabel('Mean Generation TBV')
legend(loc='upper left')




