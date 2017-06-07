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


cat = sys.argv[1]
tT_meanP = tT_mean[[cat, 'Cycle']]
tT_varP = tT_var[[cat, 'Cycle']]



tT_meanP.set_index('Cycle').plot()
plt.savefig('GenTrends_Mean_' + cat + '.pdf')
tT_varP.set_index('Cycle').plot()
plt.savefig('GenTrends_Var_' + cat + '.pdf')

