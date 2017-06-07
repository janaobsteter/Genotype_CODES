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
tgenAP.set_index('Cycle').plot()
plt.savefig('Accuracies_' + cat + '.pdf')

