# -*- coding: utf-8 -*-
import pandas as pd
from collections import defaultdict
import operator

BetaCSN = pd.DataFrame()
def makeGen(row):
    return min(row[al], row[al+1]) + max(row[al], row[al+1])
#################################################################################
#Beta kazein B
##################################################################################
B = pd.read_table('/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/IDBv03/BetaB.ped', header=None, sep=" ")
BetaCSN = pd.DataFrame({'ID': list(B[1]) })
cols = len(B.columns)
B = B[[1] + range(6,cols)]

#tukaj združi sosednja dva alela v genotipe
for al in [i for i in range(6,cols) if i % 2 == 0]:
    B.loc[:,str(al)+'g'] = B.apply(makeGen, axis=1)

#tukaj pridobi napogostejši genotip
def countOcc_max(alleles):   
    Count = defaultdict()    
    for geno in unique(alleles):
        if geno != '00':
            Count[geno] = alleles.count(geno)   
    return max(Count.iteritems(), key=operator.itemgetter(1))[0]
    
B.loc[:, 'finalGeno'] = ""
for row in B.index:
    B['finalGeno'][row] =  next(x for x in (list(B.ix[row][genoCol])) if x  != '00')
    genoCol.append(str(al) + 'g')

B.rename(columns={1:'ID', 'finalGeno':'B'}, inplace=True)
BetaCSN = pd.merge(BetaCSN, B[['ID','B']], on='ID')


##########################################################################
#Beta-CSN AB SNP2
##########################################################################
A2 = pd.read_table('/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/IDBv03/BetaA2.ped', header=None, sep=" ")
cols = len(A2.columns)
A2 = A2[[1] + range(6,cols)]

#tukaj združi sosednja dva alela v genotipe
genoCol = []
for al in [i for i in range(6,cols) if i % 2 == 0]:
    A2.loc[:,str(al)+'g'] = A2.apply(makeGen, axis=1)
    genoCol.append(str(al) + 'g')

#tukaj pridobi napogostejši genotip
A2.loc[:, 'finalGeno'] = ""
for row in A2.index:
    A2['finalGeno'][row] =  next(x for x in (list(A2.ix[row][genoCol])) if x  != '00')
    
    
A2.rename(columns={1:'ID', 'finalGeno':'A2'}, inplace=True)
BetaCSN = pd.merge(BetaCSN, A2[['ID','A2']], on='ID')

#tukaj sedaj preveri, kateri alel na B-ju je kateri
#B ne more iti v A2
Anzelak = pd.read_csv('/home/jana/Documents/F4F/AnzelakZivali.csv', header=None)
    
Anzelak.rename(columns={0:'ID'}, inplace=True)
AnzelakCSN = pd.merge(Anzelak, B[ ['ID','6g', '8g', '10g', '12g', '14g']], on='ID')
AnzelakCSN = pd.merge(AnzelakCSN, A2[ ['ID','6g', '8g', '10g', '12g', '14g']], on='ID')
AnzelakCSN.to_csv('/home/jana/Documents/F4F/AnzelakCipiRezultati.csv')

#zaenkrat lahko ločim le A2 alel od drugih - nadaljnje analize

###############################################################################
#SKUPNI
def skupniBCSN(row):
    if (row['B'] == 'BB') : #non-causative C & E, BB
        val = 'A2/A2'
    elif (row['B'] == 'AB') :
        val = 'A2/*'
    elif (row['B'] == 'AA'):
        val = '*/*'
    return val

BetaCSN.loc[:, 'SKUPEN'] = BetaCSN.apply(skupniBCSN, axis=1)
BetaCSN.to_csv('/home/jana/Documents/F4F/MlecniProteini/BetaCaseinGenotype_python.csv')

