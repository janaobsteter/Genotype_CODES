# -*- coding: utf-8 -*-
import pandas as pd
from collections import defaultdict
import operator


def makeGen(row):
    return min(row[al], row[al+1]) + max(row[al], row[al+1])
#################################################################################
#kapa kazein A
##################################################################################
AB1 = pd.read_table('/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/IDBv03/BetaLGB1.ped', header=None, sep=" ")
BetaLGB = pd.DataFrame({'ID': list(AB1[1]) })
cols = len(AB1.columns)
AB1 = AB1[[1] + range(6,cols)]

#tukaj združi sosednja dva alela v genotipe
for al in [i for i in range(6,cols) if i % 2 == 0]:
    AB1.loc[:,str(al)+'g'] = AB1.apply(makeGen, axis=1)

#tukaj pridobi napogostejši genotip
def countOcc_max(alleles):   
    Count = defaultdict()    
    for geno in unique(alleles):
        if geno != '00':
            Count[geno] = alleles.count(geno)   
    return max(Count.iteritems(), key=operator.itemgetter(1))[0]
    
AB1.loc[:, 'finalGeno'] = ""
for row in AB1.index:
    AB1['finalGeno'][row] =  next(x for x in (list(AB1.ix[row][genoCol])) if x  != '00')
    genoCol.append(str(al) + 'g')
   
    
AB1.rename(columns={1:'ID', 'finalGeno':'AB1'}, inplace=True)
BetaLGB = pd.merge(BetaLGB, AB1[['ID','AB1']], on='ID') 



##########################################################################
#Kappa-CSN AB SNP2
##########################################################################
AB2 = pd.read_table('/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/IDBv03/BetaLGB2.ped', header=None, sep=" ")
cols = len(AB2.columns)
AB2 = AB2[[1] + range(6,cols)]

#tukaj združi sosednja dva alela v genotipe
genoCol = []
for al in [i for i in range(6,cols) if i % 2 == 0]:
    AB2.loc[:,str(al)+'g'] = AB2.apply(makeGen, axis=1)
    genoCol.append(str(al) + 'g')

#tukaj pridobi napogostejši genotip
AB2.loc[:, 'finalGeno'] = ""
for row in AB2.index:
    try: 
        AB2['finalGeno'][row] =  next(x for x in (list(AB2.ix[row][genoCol])) if x  != '00' ) 
    except:
        AB2['finalGeno'][row] = '00'
    
    
AB2.rename(columns={1:'ID', 'finalGeno':'AB2'}, inplace=True)
BetaLGB = pd.merge(BetaLGB, AB2[['ID','AB2']], on='ID') 