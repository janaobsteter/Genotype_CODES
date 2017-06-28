import pandas as pd

a = pd.DataFrame({'g': [1,2], 'f1' : [7,8], 'm1': [4,5], 'm2':['f', 'm']})

print a

import matplotlib.pyplot as plt

for i in a.groupby('m2'):
    print i

"""
for i, group in a.groupby('m2'):
    plt.figure()
    group.plot(x='g', y='f1', title=str(i))
"""
for index, group in a.groupby(['m2']):
    group_agg = group.groupby(['f1']).aggregate(np.mean)
    group_agg.plot(y='g', label=index)