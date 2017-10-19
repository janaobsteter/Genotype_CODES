# -*- coding: utf-8 -*-
import pandas as pd
from collections import defaultdict
import operator

<<<<<<< HEAD
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
=======
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
>>>>>>> 1b8d6b61038ae92f43525513dcd93532f3369678

#tukaj pridobi napogostejši genotip
def countOcc_max(alleles):   
    Count = defaultdict()    
    for geno in unique(alleles):
        if geno != '00':
            Count[geno] = alleles.count(geno)   
    return max(Count.iteritems(), key=operator.itemgetter(1))[0]
    
<<<<<<< HEAD
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
=======
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
>>>>>>> 1b8d6b61038ae92f43525513dcd93532f3369678

#tukaj združi sosednja dva alela v genotipe
genoCol = []
for al in [i for i in range(6,cols) if i % 2 == 0]:
<<<<<<< HEAD
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
    if (row['AB1'] == 'AA') & (row['AB2'] == 'BB') & (row['C'] == 'BB') & (row['E'] == 'AA'): #non-causative C & E, BB
        val = 'BB'
    elif (row['AB1'] == 'BB') & (row['AB2'] == 'AA') & (row['C'] == 'BB') & (row['E'] == 'AA'):#non-causative C & E, AA
        val = 'AA'
    elif (row['AB1'] == 'AB') & (row['AB2'] == 'AB') & (row['C'] == 'BB') & (row['E'] == 'AA'):#non-causative C & E, AB
        val = 'AB'
    elif (row['AB1'] == 'AA') & (row['AB2'] == 'BB') & (row['C'] == 'AB') & (row['E'] == 'AA'): #BB, but causative C --> BC
        val = 'BC'
    elif (row['AB1'] == 'BB') & (row['AB2'] == 'AA') & (row['C'] == 'BB') & (row['E'] == 'AB'):#AA, but causative E --> AE
        val = 'AE'
    elif (row['AB1'] == 'AB') & (row['AB2'] == 'AB') & (row['C'] == 'AB') & (row['E'] == 'AA'):#AB, but causative C --> AC
        val = 'AC'
    elif (row['AB1'] == 'AB') & (row['AB2'] == 'AB') & (row['C'] == 'BB') & (row['E'] == 'AB'):#AB, but causative E --> BE
        val = 'BE'
    return val

KapaCSN.loc[:, 'SKUPEN'] = KapaCSN.apply(skupniKCSN, axis=1)
=======
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
>>>>>>> 1b8d6b61038ae92f43525513dcd93532f3369678

