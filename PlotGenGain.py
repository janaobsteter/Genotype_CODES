from selection10 import TBVPed
import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt
from pylab import legend
import os

TBV  = TBVPed(os.getcwd() + '/')
TBV.genTrend()
#plt.errorbar(x = TBV.gens, y = TBV.means, yerr = TBV.vars)

plt.plot( TBV.gens, TBV.means,  label = 'Mean Gen TBV')
plt.xlabel('Selected Generation')
plt.ylabel('Mean Generation TBV')
legend(loc='upper left')
plt.savefig('GenTrendsLast_Mean.pdf')


plt.plot(TBV.gens, TBV.vars, label = 'TBV Var')
legend(loc='upper left')
plt.xlabel('Selected Generation')
plt.ylabel('Generation TBV variance')
plt.savefig('GenTrendsLast_Var.pdf')
