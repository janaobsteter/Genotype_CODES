# -*- coding: utf-8 -*-
import pandas as pd
from collections import defaultdict
import operator
import sys

#funkcija, ki naredi genotip iz alelov
def makeGen(row):
    return min(row[al], row[al+1]) + max(row[al], row[al+1])

#funkcija, ki najde najpogostejši genotip na podvojenih SNPih
def countOcc_max(alleles):   
    Count = defaultdict()    
    for geno in unique(alleles):
        if geno != '00':
            Count[geno] = alleles.count(geno)   
    return max(Count.iteritems(), key=operator.itemgetter(1))[0]
    
#funkcija za kreiranje skupnega genotipa (končnega) na kapa kazeinu - iz vseh SNPov
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
    else:
        val = "00"
    return val

#funkcija, ki uporabi zgornje, da skreira tabelo Kapa Kazeina
def CreateKapaGenotype(GenoDir, GenoName):
    GenoDir = sys.argv[1]#'/home/jana//Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/AB/IDBv03/Genotipi_13112017/'
    GenoName = sys.argv[2]#'AB_we_bl_03112017_IDBv03_CleanIndsMarkers_'
    KapaCSN = pd.DataFrame()

    
    #################################################################################
    #kapa kazein A
    ##################################################################################
    AB1 = pd.read_table(GenoDir + GenoName + 'KapaAB_1.ped', header=None, sep=" ")
    cols = len(AB1.columns)
    AB1 = AB1[[1] + range(6,cols)]
    
    #tukaj združi sosednja dva alela v genotipe
    genoCol = []
    for al in [i for i in range(6,cols) if i % 2 == 0]:
        AB1.loc[:,str(al)+'g'] = AB1.apply(makeGen, axis=1)
        genoCol.append(str(al) + 'g')
    
    #final genotype je PRVI neničelni genotip! (ne najpogostejši!!!)    
    AB1.loc[:, 'finalGeno'] = ""
    for row in AB1.index:
        AB1['finalGeno'][row] =  next(x for x in (list(AB1.ix[row][genoCol])) if x  != '00')
    
    AB1.rename(columns={1:'ID', 'finalGeno':'AB1'}, inplace=True)
    KapaCSN = AB1[['ID','AB1']]
    
    
    ##########################################################################
    #Kappa-CSN AB SNP2
    ##########################################################################
    AB2 = pd.read_table(GenoDir + GenoName + 'KapaAB_2.ped', header=None, sep=" ")
    cols = len(AB2.columns)
    AB2 = AB2[[1] + range(6,cols)]
    
    
    #tukaj združi sosednja dva alela v genotipe
    genoCol = []
    for al in [i for i in range(6,cols) if i % 2 == 0]:
        AB2.loc[:,str(al)+'g'] = AB2.apply(makeGen, axis=1)
        genoCol.append(str(al) + 'g')
    
    #tukaj pridobi PRVI NENULNI genotip
    AB2.loc[:, 'finalGeno'] = ""
    for row in AB2.index:
        AB2['finalGeno'][row] =  next(x for x in (list(AB2.ix[row][genoCol])) if x  != '00')
        
    
    AB2.rename(columns={1:'ID', 'finalGeno':'AB2'}, inplace=True)
    KapaCSN = pd.merge(KapaCSN, AB2[['ID','AB2']], on='ID') 
    
    ##########################################################################
    #Kappa-CSN E
    ##########################################################################
    C = pd.read_table(GenoDir + GenoName + 'KapaC.ped', header=None, sep=" ")
    cols = len(C.columns)
    C = C[[1] + range(6,cols)]
    
    #tukaj združi sosednja dva alela v genotipe
    genoCol = []
    for al in [i for i in range(6,cols) if i % 2 == 0]:
        C.loc[:,str(al)+'g'] = C.apply(makeGen, axis=1)
        genoCol.append(str(al) + 'g')
    
    #tukaj pridobi PRVI NENIČELNI genotip
    C.loc[:, 'finalGeno'] = ""
    for row in C.index:
        C['finalGeno'][row] =  next(x for x in (list(C.ix[row][genoCol])) if x  != '00')
        
    
    
    C.rename(columns={1:'ID', 'finalGeno':'C'}, inplace=True)
    KapaCSN = pd.merge(KapaCSN, C[['ID','C']], on='ID') 
    
    
    
    ##########################################################################
    #Kappa-CSN E
    ##########################################################################
    E = pd.read_table(GenoDir + GenoName + 'KapaE.ped', header=None, sep=" ")
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
        
    
    E.rename(columns={1:'ID', 'finalGeno':'E'}, inplace=True)
    KapaCSN = pd.merge(KapaCSN, E[['ID','E']], on='ID') 
    
    ##########################################################################
    #Kappa-CSN E
    ##########################################################################
    I = pd.read_table(GenoDir + GenoName + 'KapaI.ped', header=None, sep=" ")
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
        
        
    
    I.rename(columns={1:'ID', 'finalGeno':'I'}, inplace=True)
    KapACSN = pd.merge(KapaCSN, I[['ID','I']], on='ID') 

###############################################################################



    KapaCSN.loc[:, 'SKUPEN'] = KapaCSN.apply(skupniKCSN, axis=1)
    KapaCSN.to_csv(GenoDir + GenoName + 'KappaCaseinGenotype_python.csv', index = None)
