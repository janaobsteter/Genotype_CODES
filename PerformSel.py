from __future__ import division
import pandas as pd
import numpy as np


#########################################
#newborn gen
#########################################
nT = 20000
nrT = 0.19
nF = 18000
nrF = int(0.15 * nF)
tel12 = 0.11 * nF
pt12 = 0.035 * nF
tel24 = 0.08 * nF
pt24 = 0.03 * nF
k = 0.625 * nF

#ped = pd.read_csv("~/Documents/PhD/Simulaton/Pedigrees/Pedigree_10burnIn_Gen1.txt", delim_whitespace=True,header='infer')
ped = pd.read_csv("~/Documents/PhD/Simulaton/Pedigrees/PedPython.txt")
#ped = pd.DataFrame(ped)
#ped = ped[ped.columns[[1,2,3,8]]]

nr = range(int((max(ped['Generation']) - 1)*(nT*nrT)),int(((max(ped['Generation']))*(nT*nrT))))

ped['cat'] = ""
#ped.loc[ped.Generation == (min(ped['Generation'])), 'cat'] = "nr"
ped['active'] = 2
ped['sex'] = ""

gen_in_ped = set(ped.Generation)


#FEMALES, Ggen = max(gen)
ped0 = ped.loc[ped.Generation == max(ped.Generation)]
ped0[0:(nF*nrF), 'sex'] = 2