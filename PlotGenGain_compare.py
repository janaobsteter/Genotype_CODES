# -*- coding: utf-8 -*-
from selection10 import TBVPed
import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt
from pylab import legend
import os


class TBVPed(object):  # to je tabela za grafiranje genetskih trendov čez populacije
    def __init__(self, AlphaSimDir):
        self.AlphaSimDir = AlphaSimDir

    def genTrend(self, table):
        TBVtable = pd.read_table(table, sep='\s+')
        gens = list(set(TBVtable.Generation))
        means = TBVtable.gvNormUnres1.groupby(TBVtable.Generation).aggregate(np.mean)
        vars = TBVtable.gvNormUnres1.groupby(TBVtable.Generation).aggregate(np.var)
        return gens, means, vars

    def catTrend(self):
        cat = list(
            set(pd.read_table(self.AlphaSimDir + '/SimulatedData/PedigreeAndGeneticValues_cat.txt', sep='\s+').cat))
        TBVtable = pd.read_table(self.AlphaSimDir + 'SimulatedData/PedigreeAndGeneticValues_cat.txt', sep='\s+')
        gen = max(set(TBVtable.Generation))
        means = TBVtable.gvNormUnres1.groupby(TBVtable.cat).aggregate(np.mean)
        vars = TBVtable.gvNormUnres1.groupby(TBVtable.cat).aggregate(np.var)
        return pd.DataFrame({'cat': cat, str(gen) + '_mean': means, str(gen) + '_vars': vars})

    def GenTrend(self):
        gens = list(
            set(pd.read_table(self.AlphaSimDir + '/SimulatedData/PedigreeAndGeneticValues_cat.txt',
                              sep='\s+').Generation))
        TBVtable = pd.read_table(self.AlphaSimDir + 'SimulatedData/PedigreeAndGeneticValues_cat.txt', sep='\s+')
        gen = max(set(TBVtable.Generation))
        means = TBVtable.gvNormUnres1.groupby(TBVtable.Generation).aggregate(np.mean)
        vars = TBVtable.gvNormUnres1.groupby(TBVtable.Generation).aggregate(np.var)
        return pd.DataFrame({'Gen': gens, str(gen) + '_mean': means, str(gen) + '_vars': vars})



TBV  = TBVPed(os.getcwd() + '/')
gens, means, vars = TBV.genTrend('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Class/SimulatedData/PedigreeAndGeneticValues.txt')

TBV  = TBVPed(os.getcwd() + '/')
gensG, meansG, varsG = TBV.genTrend('/home/jana/bin/AlphaSim1.05Linux/SimulatedData/PedigreeAndGeneticValues.txt')
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




plt.plot( Means.Gen, Means.Stand,  label = 'Mean Gen TBV_class')



plt.plot( Means.GenG, Means.StandG,  label = 'Mean Gen TBV_gen')
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



