import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt

genA = pd.read_csv('Accuracies_Cat.csv')
genA.index = genA.cat
genA = genA.drop('cat', axis=1)
tgenA = np.transpose(genA)
tgenA.loc[:, 'Cycle'] = list(tgenA.index)

#cats = [i for i in ['pBM', 'pb','gpb','genTest', 'k', 'pripust1', 'pripust2', 'mladi'] if i in tT_mean.columns]
cats = ['mladi']

for cat in cats:
    tgenAP = tgenA[[cat, 'Cycle']]
    plt.plot(tgenAP.Cycle, tgenAP.loc[:,cat], label = cat)

plt.xlabel('Selected Generation')
plt.ylabel('Mean Generation TBV')
legend(loc='upper left')
plt.savefig('Accuracies_' + cat + '.pdf')

print 'Created plot: Accuracies_PathsOfSel.pdf'
