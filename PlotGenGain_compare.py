# -*- coding: utf-8 -*-
from selection10 import TBVPed
import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt
from pylab import legend
import os



TBV  = TBVPed(os.getcwd() + '/')
gens, means, vars = TBV.genTrend('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Class/SimulatedData/PedigreeAndGeneticValues.txt', 41, 61)


TBV  = TBVPed(os.getcwd() + '/')
gensA, meansA, varsA = TBV.genTrend('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenFather/SimulatedData/PedigreeAndGeneticValues.txt', 41, 61)


TBV  = TBVPed(os.getcwd() + '/')
gensG, meansG, varsG = TBV.genTrend('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenFatherReference/SimulatedData/PedigreeAndGeneticValues.txt', 41, 61)
#plt.errorbar(x = TBV.gens, y = TBV.means, yerr = TBV.vars)

Means = pd.DataFrame()
Means.loc[:,'Gen'] = [int(x) for x in gens]
Means.loc[:,'Mean'] = list(means)
Means.loc[:,'Var'] =  list(vars)
Means.loc[:,'Stand'] = (Means.Mean - Means.Mean[0]) / sqrt(Means.Var[0])
Means.loc[:,'GenG'] = [int(x) for x in gensG]
Means.loc[:,'MeanG'] = list(meansG)
Means.loc[:,'VarG'] =  list(varsG)
Means.loc[:,'StandG'] = (Means.MeanG - Means.MeanG[0]) / sqrt(Means.VarG[0])
Means.loc[:,'GenA'] = [int(x) for x in gensA]
Means.loc[:,'MeanA'] = list(meansA)
Means.loc[:,'VarA'] =  list(varsA)
Means.loc[:,'StandA'] = (Means.MeanA - Means.MeanA[0]) / sqrt(Means.VarA[0])



plt.plot( Means.Gen, Means.Stand,  label = 'Mean Gen TBV_class')

plt.plot( Means.Gen, Means.Stand,  label = 'Mean Gen TBV_genFatherRef')

plt.plot( Means.GenG, Means.StandG,  label = 'Mean Gen TBV_genFather')
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

plt.plot(gensEBV, meansEBV, label='EBV')
plt.xlabel('Selected Generation')
plt.ylabel('Mean Generation TBV')
legend(loc='upper left')



