import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt

Cycle = str(sys.argv[1])

genA = pd.read_csv('Accuracies_Cat.csv')
genA.index = genA.cat
genA = genA.drop('cat', axis=1)
genAP = genA[[Cycle]]

genAP.plot(marker='o', linestyle='')
plt.savefig('AccuraciesCat_cycle' + Cycle + '.pdf')
