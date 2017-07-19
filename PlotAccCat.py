import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt

genA = pd.read_csv('Accuracies_Cat.csv')
genA.index = genA.cat
genA = genA.drop('cat', axis=1)
tgenA = np.transpose(genA)
tgenA.loc[:, 'Cycle'] = list(tgenA.index)

cat = sys.argv[1]

tgenAP = tgenA[[cat, 'Cycle']]
plot(tgenAP.Cycle, tgenAP[[cat]])
tgenAP.set_index('Cycle').plot()
plt.savefig('Accuracies_' + cat + '.pdf')

print 'Created plot: Accuracies_' + cat+ '.pdf'




f, axarr = plt.subplots(nrows=3, ncols=2, sharex=True,sharey=True)
f.suptitle('Accuracies EBV')
cat='pb'
tgenAP = tgenA[[cat, 'Cycle']]
axarr[0,0].plot(tgenAP.Cycle, tgenAP[[cat]])
axarr[0,0].set_title(cat)
cat='genTest'
tgenAP = tgenA[[cat, 'Cycle']]
axarr[1,0].plot(tgenAP.Cycle, tgenAP[[cat]])
axarr[1,0].set_title(cat)
cat='mladi'
tgenAP = tgenA[[cat, 'Cycle']]
axarr[2,0].plot(tgenAP.Cycle, tgenAP[[cat]])
axarr[2,0].set_title(cat)
cat='pBM'
tgenAP = tgenA[[cat, 'Cycle']]
axarr[0,1].plot(tgenAP.Cycle, tgenAP[[cat]])
axarr[0,1].set_title(cat)
cat='bm'
tgenAP = tgenA[[cat, 'Cycle']]
axarr[1,1].plot(tgenAP.Cycle, tgenAP[[cat]])
axarr[1,1].set_title(cat)
cat='k'
tgenAP = tgenA[[cat, 'Cycle']]
axarr[2,1].plot(tgenAP.Cycle, tgenAP[[cat]])
axarr[2,1].set_title(cat)

axarr[0].set_title('Acuuracies EBV')
