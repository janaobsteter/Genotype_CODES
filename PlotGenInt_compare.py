#from selection10 import genInterval
import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt
import os


f, axarr = plt.subplots(nrows=3, ncols=2, sharex=True,sharey=True)
f.suptitle('GenInt')
xs = [0,1,2,0,1,2]
ys = [0,0,0,1,1,1]

os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Class/')
GenInts = pd.read_csv('GenInts.txt', sep=" ")
GenInts.loc[:, 'label'] = GenInts['line'] + GenInts['sex']
Gens = GenInts[['Gen', 'genInt', 'label']]
plot1 = Gens.groupby(['Gen', 'label']).sum().unstack()
axarr[xs[0],ys[0]].plot(plot1)
axarr[xs[0],ys[0]].set_title("Conventional")


os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenSLO/')
GenInts = pd.read_csv('GenInts.txt', sep=" ")
GenInts.loc[:, 'label'] = GenInts['line'] + GenInts['sex']
Gens = GenInts[['Gen', 'genInt', 'label']]
plot2 = Gens.groupby(['Gen', 'label']).sum().unstack()
axarr[xs[1],ys[0]].plot(plot2)
axarr[xs[1],ys[0]].set_title("Conventional_SLO")

os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenSplosnaPop/')
GenInts = pd.read_csv('GenInts.txt', sep=" ")
GenInts.loc[:, 'label'] = GenInts['line'] + GenInts['sex']
Gens = GenInts[['Gen', 'genInt', 'label']]
plot3 = Gens.groupby(['Gen', 'label']).sum().unstack()
axarr[xs[2],ys[0]].plot(plot3)
axarr[xs[2],ys[0]].set_title("GenBulls on Other Cows")

os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenSLO_BmGen/')
GenInts = pd.read_csv('GenInts.txt', sep=" ")
GenInts.loc[:, 'label'] = GenInts['line'] + GenInts['sex']
Gens = GenInts[['Gen', 'genInt', 'label']]
plot4 = Gens.groupby(['Gen', 'label']).sum().unstack()
axarr[xs[0],ys[5]].plot(plot4)
axarr[xs[0],ys[5]].set_title("GenBulls on Bull Dams")

os.chdir('/home/jana/bin/AlphaSim1.05Linux//')
GenInts = pd.read_csv('GenInts.txt', sep=" ")
GenInts.loc[:, 'label'] = GenInts['line'] + GenInts['sex']
Gens = GenInts[['Gen', 'genInt', 'label']]
plot5 = Gens.groupby(['Gen', 'label']).sum().unstack()
axarr[xs[2],ys[5]].plot(plot5)
axarr[xs[2],ys[5]].set_title("GenBulls on All Cows")
legend(['dam>dam', 'dam>sire', 'sire>dam', 'sire>sire'], loc='upper left')



plt.xlabel('Year of birth')
plt.ylabel('Mean Generation TBV')
legend(loc='upper left')