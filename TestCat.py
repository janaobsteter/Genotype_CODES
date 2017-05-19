# -*- coding: utf-8 -*-
from collections import defaultdict
import pandas as pd
import numpy as np
import selection10
reload(selection10)
from selection10 import *

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
baba['EBV'] = True
baba['gEBV'] = False

genPed = pd.DataFrame(columns=['Generation', 'Indiv', 'Father', 'Mother', 'EBV'])
genPed['Generation'] = list(chain.from_iterable([[i] * 8640 for i in range(1,11)]))
genPed['Indiv'] = range(1,86401)
genPed['Father'] = 0
genPed['Mother'] = 0
genPed['EBV'] = np.random.uniform(low=0.004, high=1.5, size = 86400)
genPed.to_csv('/home/jana/Genotipi/Genotipi_CODES/GenericPed_67000.txt', index=None)
genPed.to_csv('/home/jana/Genotype_CODES/GenericPed_86400.txt', index=None)

genPed.to_csv('/home/jana/Genotipi/Genotipi_CODES/GenericPed_1000.txt', index=None)

import selectio8n
reload(selectio8n)
from selectio8n import *
Ped = pedigree('/home/jana/Genotipi/Genotipi_CODES/GenericPed_67000.txt')
Ped = pedigree('/home/jana/Genotipi/Genotipi_CODES/GenericPed_1000.txt')

Gender = pd.DataFrame({'Generation' : list(chain.from_iterable([[i] * 8640 for i in range(1,11)])), 'Indiv': range(1, 86401), 'Gender': [1,2] * 43200})
Gender.to_csv('/home/jana/bin/AlphaSim1.05Linux/SimulatedData/Gender.txt', index=None, sep='\t')#

IndCat = pd.DataFrame()
ped, c, s, a = nastavi_cat('/home/jana/Genotipi/Genotipi_CODES/GenericPed_67000.txt', **baba) #nastavi kategorije 
krogov = 2
ped, c, s, a = nastavi_cat('/home/jana/Genotipi/Genotipi_CODES/GenericPed_1000.txt', **baba) #nastavi kategorije 
ped, c, s, a = nastavi_cat('/home/jana/Genotype_CODES/GenericPed_86400.txt', **baba) #nastavi kategorije 



IndCat['Indiv'] = ped.ped.Indiv
IndCat['catBurnIN'] = ped.ped.cat

for krog in range(krogov): #ponavljaj kolikor krogov selekcije hočeš
    ped, c, s, a = selekcija_total("/home/jana/PedTotal.txt", **baba)
    IndCat[str('cat' + str(krog))] = ped.ped.cat