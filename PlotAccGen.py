import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt

genA = pd.read_csv('Accuracies_Gen.csv')
genA.columns = [np.count_nonzero(~np.isnan(genA[[i]])) for i in range(len(genA.columns))]
tgenA = np.transpose(genA)
tgenA.columns = ['Gen' + str(i +1) for i in list(tgenA.columns)]
tgenA.loc[:, 'Cycle'] = list(tgenA.index)

firstGen = int(sys.argv[1])
lastGen = [int(sys.argv[2]) if sys.argv[2] != 999 else len(genA.columns)]

tgenAP = tgenA.ix[:,firstGen-1:lastGen[0]]
tgenAP.loc[:, 'Cycle'] = list(tgenA.Cycle)
tgenAP.set_index('Cycle').plot()
plt.savefig('AccuraciesGen.pdf')