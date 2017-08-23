#from selection10 import genInterval
import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt
import os


f, axarr = plt.subplots( sharex=True, gridspec_kw={'hspace':0}, nrows=3, ncols=2)
#f.suptitle('GenInt')
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
axarr[xs[1],ys[0]].set_title("Genomic a")

os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenSplosnaPop/')
GenInts = pd.read_csv('GenInts.txt', sep=" ")
GenInts.loc[:, 'label'] = GenInts['line'] + GenInts['sex']
Gens = GenInts[['Gen', 'genInt', 'label']]
plot3 = Gens.groupby(['Gen', 'label']).sum().unstack()
axarr[xs[2],ys[0]].plot(plot3)
axarr[xs[2],ys[0]].set_title("GenBulls on Other Cows")
axarr[xs[2],ys[0]].set_title("Genomic b")

os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenSLO_BmGen/')
GenInts = pd.read_csv('GenInts.txt', sep=" ")
GenInts.loc[:, 'label'] = GenInts['line'] + GenInts['sex']
Gens = GenInts[['Gen', 'genInt', 'label']]
plot4 = Gens.groupby(['Gen', 'label']).sum().unstack()
axarr[xs[1],ys[5]].plot(plot4)
axarr[xs[1],ys[5]].set_title("GenBulls on Bull Dams")
axarr[xs[1],ys[5]].set_title("Genomic c")

os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Gen/')
GenInts = pd.read_csv('GenInts.txt', sep=" ")
GenInts.loc[:, 'label'] = GenInts['line'] + GenInts['sex']
Gens = GenInts[['Gen', 'genInt', 'label']]
plot5 = Gens.groupby(['Gen', 'label']).sum().unstack()
axarr[xs[2],ys[5]].plot(plot5)
axarr[xs[2],ys[5]].set_title("GenBulls on All Cows")
axarr[xs[2],ys[5]].set_title("Genomic d")
legend(['dam>dam', 'dam>sire', 'sire>dam', 'sire>sire'], bbox_to_anchor=(1.04,1), loc='upper left')

"""
os.chdir('/home/jana/bin/AlphaSim1.05Linux//')
GenInts = pd.read_csv('GenInts.txt', sep=" ")
GenInts.loc[:, 'label'] = GenInts['line'] + GenInts['sex']
Gens = GenInts[['Gen', 'genInt', 'label']]
plot5 = Gens.groupby(['Gen', 'label']).sum().unstack()
axarr[xs[2],ys[5]].plot(plot5)
axarr[xs[2],ys[5]].set_title("GenBulls on BM 2")
legend(['dam>dam', 'dam>sire', 'sire>dam', 'sire>sire'], bbox_to_anchor=(1.04,1), loc='upper left')
"""


f.text(0.5, 0.04, 'Generation of birth', ha='center', va='center', size=14)
f.text(0.06, 0.5, 'Generation interval [years]', ha='center', va='center', rotation='vertical', size=14)
