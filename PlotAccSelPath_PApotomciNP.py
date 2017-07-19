from scipy.stats.stats import pearsonr   
import pandas as pd
import numpy as np
from collections import defaultdict
import math
import matplotlib.pyplot as plt
from pylab import legend
import os

os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Class/')
os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Gen/')
os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenReference/')
os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenFather/')
os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenFatherReference/')
os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenBM_UpdatedRef/')
accuraciesEBV = pd.read_table('AccuraciesEBVPA.csv', sep=",")
cats = [i for i in ['pBM', 'pb','gpb','mladi'] if i in list(accuraciesEBV.Cat)]
cats = [i for i in ['pBM', 'pb','gpb','genTest', 'k', 'pripust1', 'pripust2', 'mladi'] if i in list(accuraciesEBV.Cat)]

cats =  ['gpb']


for cat in cats:
    accCat = accuraciesEBV.loc[accuraciesEBV.Cat == cat]
    plt.plot(accCat.Cycle, accCat.corEBV, label = cat)
plt.xlabel('Selected Generation')
plt.ylabel('EBV accuracies')
legend(loc='upper left')
plt.savefig('AccuraciesSelPaths.pdf')



for cat in cats:
    trendCat = accuraciesEBV.loc[accuraciesEBV.Cat == cat]
    plt.plot(trendCat.Cycle, trendCat.meanEBV, label = cat)
plt.xlabel('Selected Generation')
plt.ylabel('Mean Generation EBV')
legend(loc='upper left')
plt.savefig('GeneticTrendEBV_SelPath.pdf')
plt.close()

for cat in cats:
    trendCat = accuraciesEBV.loc[accuraciesEBV.Cat == cat]
    plt.plot(trendCat.Cycle, trendCat.meanTBV, label = cat)
plt.xlabel('Selected Generation')
plt.ylabel('Mean Generation TBV')
legend(loc='upper left')
plt.savefig('GeneticTrendTBV_SelPath.pdf')
plt.close()

 

accCat = accuraciesEBV.loc[accuraciesEBV.Cat == 'potomciNP']
fit = np.polyfit(accCat.Cycle, accCat.corPA, deg=1)
plt.scatter(accCat.Cycle, accCat.corPA)
plt.plot(accCat.Cycle, fit[0] *accCat.Cycle + fit[1], color='red')
plt.xlabel('Selected Generation')
plt.ylabel('PA accuracies')
legend(loc='upper left')
plt.savefig('AccuraciesPA_potomciNP.pdf')
plt.close()


trendCat = accuraciesEBV.loc[accuraciesEBV.Cat == 'potomciNP']
fit = np.polyfit(trendCat.Cycle, trendCat.meanPA, deg=1)
plt.scatter(trendCat.Cycle, trendCat.meanPA)
plt.plot(trendCat.Cycle, fit[0] *trendCat.Cycle + fit[1], color='red')
plt.xlabel('Selected Generation')
plt.ylabel('Mean PA')
legend(loc='upper left')
plt.savefig('GeneticTrendPA_SelPath.pdf')
plt.close()
print 'Created plots: AccuraciesSelPaths.pdf and AccuraciesPA_potomciNP.pdf \n GeneticTrendEBV_SelPath.pdf \n GeneticTrendTBV_SelPath.pdf \n GeneticTrendBA_SelPath.pdf'