import pandas as pd
import sys
import numpy as np
from numpy import transpose
import matplotlib.pyplot as plt
numberGen = sys.argv[1]

genA = pd.read_csv('Accuracies_Gen.csv')
genA.columns = [np.count_nonzero(~np.isnan(genA[[i]])) for i in range(len(genA.columns))]
tgenA = transpose(genA)
tgenA.columns = ['Gen' + str(i) for i in list(tgenA.columns)]
tgenA.loc[:, 'Cycle'] = list(tgenA.index)

tgenAP = tgenA.ix[:,10:]
tgenAP.set_index('Cycle').plot()
plt.savefig('AccuraciesGen.pdf')

print "Created plot: AccuraciesGen.pdf"
