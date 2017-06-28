import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt

T = pd.read_csv('GenTrends_cat.csv')
T.index = T.cat
T = T.drop('cat', axis=1)

tT = np.transpose(T)
tT.loc[:,'Cycle'] = [i.strip('_vars').strip('_mean') for i in list(tT.index)]
tT_mean = tT.ix[0::2,:]
tT_var = tT.ix[1::2,:]

cats = [i for i in ['pBM', 'pb','gpb','genTest', 'k', 'pripust1', 'pripust2', 'mladi'] if i in tT_mean.columns]

for cat in cats:
    tT_meanP = tT_mean[[cat, 'Cycle']]
    tT_varP = tT_var[[cat, 'Cycle']]
    plt.plot(tT_meanP.Cycle, tT_meanP.loc[:,cat], label = cat)

plt.xlabel('Selected Generation')
plt.ylabel('Mean Generation TBV')
legend(loc='upper left')
plt.savefig('GenTrends_Mean_PathOfSel.pdf')

for cat in cats:
    tT_meanP = tT_mean[[cat, 'Cycle']]
    tT_varP = tT_var[[cat, 'Cycle']]
    plt.plot(tT_varP.Cycle, tT_varP.loc[:,cat], label = cat)

plt.xlabel('Selected Generation')
plt.ylabel('Mean Generation TBV')
legend(loc='upper left')
plt.savefig('GenTrends_Var_PathOfSel.pdf')

print 'Created plots: GenTrends_Mean_' + cat + '.pdf and GenTrends_Var_' + cat + '.pdf'
