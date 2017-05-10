# -*- coding: utf-8 -*-
from collections import defaultdict
import pandas as pd
import numpy as np

categories = defaultdict(list)
catDF = pd.read_csv('/home/jana/bin/AlphaSim1.05Linux//BurnInPLus19Gen_05052017/Categories_gen11DF.csv')
for cat in catDF.columns:
    values = [int(i) for i in catDF[cat] if not math.isnan(i)]
    categories[cat] = values

baba = defaultdict()  
sel = pd.read_csv('/home/jana/SelectionParamTEST.csv', header=None)
for i in range(len(sel)):
    if (sel[0][i]) not in  ['BurnInYN','EBV','gEBV']:
        try:
            baba[sel[0][i]] = int(sel[1][i])
        except:
            baba[sel[0][i]] = float(sel[1][i])  

baba['BurnInYN'] = False
baba['EBV'] = False
baba['gEBV'] = True

genPed = pd.DataFrame(columns=['Generation', 'Indiv', 'Father', 'Mother', 'EBV'])
genPed['Generation'] = list(chain.from_iterable([[i] * 6700 for i in range(10)]))
genPed['Indiv'] = range(1,67001)
genPed['Father'] = 0
genPed['Mother'] = 0
genPed['EBV'] = np.random.uniform(low=0.004, high=1.5, size = 67000)
genPed.to_csv('/home/jana/Genotipi/Genotipi_CODES/GenericPed_67000.txt', index=None)

import selection7
reload(selection7)
from selection7 import *
Ped = pedigree('/home/jana/Genotipi/Genotipi_CODES/GenericPed_67000.txt')

#

ped, c, s, a = nastavi_cat('/home/jana/Genotipi/Genotipi_CODES/GenericPed_67000.txt', **baba) #nastavi kategorije 
krogov = 5
for krog in krogov: #ponavljaj kolikor krogov selekcije hočeš
    ped, c, s, a = selekcija_total("/home/jana/PedTotal.txt", **baba)