# -*- coding: utf-8 -*-
import pandas as pd
from collections import defaultdict
import operator

KapaCSN = pd.DataFrame()
def makeGen(row):
    return min(row[al], row[al+1]) + max(row[al], row[al+1])
#################################################################################
#kapa kazein A
##################################################################################
AB1 = pd.read_table('/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/IDBv03/KapaAB_1.ped', header=None, sep=" ")
AB1 = AB1[[1,6,7,8,9,10,11,12,13]]

#tukaj združi sosednja dva alela v genotipe
for al in [i for i in range(6,13) if i % 2 == 0]:
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
    
KapaCSN.loc[:, 'AB1'] = AB1['finalGeno']


##########################################################################
#Kappa-CSN AB SNP2
##########################################################################
AB2 = pd.read_table('/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/IDBv03/KapaAB_2.ped', header=None, sep=" ")
AB2 = AB2[[1] + range(6,26)]

#tukaj združi sosednja dva alela v genotipe
genoCol = []
for al in [i for i in range(6,26) if i % 2 == 0]:
    AB2.loc[:,str(al)+'g'] = AB2.apply(makeGen, axis=1)
    genoCol.append(str(al) + 'g')

#tukaj pridobi napogostejši genotip
AB2.loc[:, 'finalGeno'] = ""
for row in AB2.index:
    AB2['finalGeno'][row] =  next(x for x in (list(AB2.ix[row][genoCol])) if x  != '00')
    
    
KapaCSN.loc[:, 'AB2'] = AB2['finalGeno']

##########################################################################
#Kappa-CSN E
##########################################################################
C = pd.read_table('/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/IDBv03/KapaC.ped', header=None, sep=" ")
cols = len(C.columns)
C = C[[1] + range(6,cols)]

#tukaj združi sosednja dva alela v genotipe
genoCol = []
for al in [i for i in range(6,cols) if i % 2 == 0]:
    C.loc[:,str(al)+'g'] = C.apply(makeGen, axis=1)
    genoCol.append(str(al) + 'g')

#tukaj pridobi napogostejši genotip
C.loc[:, 'finalGeno'] = ""
for row in C.index:
    C['finalGeno'][row] =  next(x for x in (list(C.ix[row][genoCol])) if x  != '00')
    
    
KapaCSN.loc[:, 'C'] = C['finalGeno']


##########################################################################
#Kappa-CSN E
##########################################################################
E = pd.read_table('/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/IDBv03/KapaE.ped', header=None, sep=" ")
cols = len(E.columns)
E = E[[1] + range(6,cols)]

#tukaj združi sosednja dva alela v genotipe
genoCol = []
for al in [i for i in range(6,cols) if i % 2 == 0]:
    E.loc[:,str(al)+'g'] = E.apply(makeGen, axis=1)
    genoCol.append(str(al) + 'g')

#tukaj pridobi napogostejši genotip
E.loc[:, 'finalGeno'] = ""
for row in E.index:
    E['finalGeno'][row] =  next(x for x in (list(E.ix[row][genoCol])) if x  != '00')
    
    
KapaCSN.loc[:, 'E'] = E['finalGeno']

##########################################################################
#Kappa-CSN E
##########################################################################
I = pd.read_table('/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/IDBv03/KapaI.ped', header=None, sep=" ")
cols = len(I.columns)
I = I[[1] + range(6,cols)]

#tukaj združi sosednja dva alela v genotipe
genoCol = []
for al in [i for i in range(6,cols) if i % 2 == 0]:
    I.loc[:,str(al)+'g'] = I.apply(makeGen, axis=1)
    genoCol.append(str(al) + 'g')

#tukaj pridobi napogostejši genotip
I.loc[:, 'finalGeno'] = ""
for row in I.index:
    I['finalGeno'][row] = next(x for x in (list(I.ix[row][genoCol])) if x  != '00')
    
    
KapaCSN.loc[:, 'I'] = I['finalGeno']

###############################################################################
#SKUPNI
def skupniKCSN(row):
    if (row['AB1'] == 'AA') & (row['AB2'] == 'BB') & (row['C'] == 'BB') & (row['E'] == 'AA'):
        val = 'BB'
    elif (row['AB1'] == 'BB') & (row['AB2'] == 'AA') & (row['C'] == 'BB') & (row['E'] == 'AA'):
        val = 'AA'
    elif (row['AB1'] == 'AA') & (row['AB2'] == 'BB') & (row['C'] == 'AB') & (row['E'] == 'AA'):
        val = 'BC'
    elif (row['AB1'] == 'BB') & (row['AB2'] == 'AA') & (row['C'] == 'BB') & (row['E'] == 'AB'):
        val = 'AE'
    elif (row['AB1'] == 'AB') & (row['AB2'] == 'AB') & (row['C'] == 'AB') & (row['E'] == 'AA'):
        val = 'AC'
    elif (row['AB1'] == 'AB') & (row['AB2'] == 'AB') & (row['C'] == 'BB') & (row['E'] == 'AB'):
        val = 'BE'
    return val

KapaCSN.loc[:, 'SKUPEN'] = KapaCSN.apply(skupniKCSN, axis=1)

