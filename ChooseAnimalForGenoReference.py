# -*- coding: utf-8 -*-
import random
import pandas as pd
import os
from itertools import chain

os.chdir('/home/jana/bin/AlphaSim1.05Linux/REALFillIn20BurnIn20/')
ks = pd.DataFrame({'Indiv': sorted(list(set(chain.from_iterable([random.sample((list(pd.read_table('Categories_gen' + str(i) + 'DF.csv', sep=",")['k'].dropna().astype(int))), 2500) for i in range(35, 41)]))))})
ped = pd.read_table(AlphaSimDir +'/SimulatedData/PedigreeAndGeneticValues_cat.txt', sep='\s+')
ks = pd.merge(ks, ped, on='Indiv')
ks.loc[ks.Generation.isin(range(max(ks.Generation)-6, max(ks.Generation)+1)), 'Indiv'].to_csv('//home/jana/bin/AlphaSim1.05Linux/IndForGeno.txt', index=None, header=None, sep='\n')

ped[ped.Indiv.isin(ks.loc[ks.Generation.isin(range(max(ks.Generation)-6, max(ks.Generation)+1)), 'Indiv'])].Generation.value_counts()

##################################
#test update of the ref
####################################
if os.path.isfile('IndForGeno.txt'):
    #first obtain new animals for genotypisation
    pd.DataFrame({0: sorted(list(set
                                    (chain.from_iterable([ped.catCurrent_indiv_sex_criteria(x, sex, xC,
                                                                                            int((len(ped.catCurrent_indiv_sex(x, sex)) * (xP / 100))))
                                                        if xP != 100 else ped.catCurrent_indiv_sex(x, sex) for (x, xP, xC, sex) in genotypedCat]))))}).to_csv(        #if ped.sel == 'class': #ZDJ TEGA NI, KER JE POSEBEJ FILE!
        'IndForGeno_new.txt', index=None, header=None)
    #then remove old animals from the previous file
    inds = pd.read_table('IndForGeno.txt', header=None, names=['Indiv'])
    ped = pd.read_table(AlphaSimDir +'/SimulatedData/PedigreeAndGeneticValues_cat.txt', sep='\s+')
    inds = pd.merge(inds, ped[['Indiv', 'Generation', 'sex']], on='Indiv')#obtain generation of the previously genotyped Individuals
#    inds = pd.merge(inds, ped[['Indiv', 'sex']], on='Indiv')#obtain sex of the previously genotyped Individuals
    sexDF = inds.loc[(inds.sex == removesex)]
    pd.concat([inds.loc[inds.sex != removesex], sexDF.loc[
        sexDF.Generation.isin(sorted(list(set(sexDF.Generation)))[rmNbGen:])]]).sort_values(by='Indiv')['Indiv'].to_csv(
        'IndForGeno.txt', index=None, header=None)
    os.system("grep -v -f IndForGeno.txt IndForGeno_new.txt > uniqNew && mv uniqNew IndForGeno_new.txt") #obtain only new ones - do you dont get duplicate cows - najbrž nepotrebno
    os.system(
        'cat IndForGeno_new.txt IndForGeno.txt | sort -n| uniq > IndGenTmp && mv IndGenTmp IndForGeno.txt')
else:
    pd.DataFrame({0: sorted(list(set
                                    (chain.from_iterable([ped.catCurrent_indiv_sex_criteria(x, sex, xC,
                                                                                            int((len(ped.catCurrent_indiv_sex(x, sex)) * (xP / 100.0))))
                                                        for (x, xP, xC, sex) in genotypedCat]))))}).to_csv(
        'IndForGeno.txt', index=None, header=None)






if os.path.isfile(AlphaSimDir + 'GenoFile.txt'): #if GenoFile.txt exists, only add the newIndividuals for genotypisation
    os.system(
        'grep -Fwf IndForGeno_new.txt ' + chipFile + ' > ChosenInd.txt')  # only individuals chosen for genotypisation - ONLY NEW - LAST GEN!
    os.system("sed 's/^ *//' ChosenInd.txt > ChipFile.txt")  # Remove blank spaces at the beginning
    os.system("cut -f1 -d ' ' ChipFile.txt > Individuals.txt")  # obtain IDs
    os.system('''awk '{$1=""; print $0}' ChipFile.txt | sed 's/ //g' > Snps.txt''')  # obtain SNP genotypes
    os.system(
        r'''paste Individuals.txt Snps.txt | awk '{printf "%- 10s %+ 15s\n",$1,$2}' > GenoFile_new.txt''')  # obtain SNP genotypes of the last generation
    os.system('grep -Fwf IndForGeno.txt GenoFile.txt  > GenoFile_Oldtmp && mv GenoFile_Oldtemp GenoFile.txt') #here obtain updated old reference - removed one generation
    os.system("cat GenoFile.txt GenoFile_new.txt > GenoFileTmp && mv GenoFileTmp GenoFile.txt")
    os.system("less GenoFile.txt | sort -n | uniq > Genotmp && mv Genotmp GenoFile.txt")
else:
    os.system(
        'grep -Fwf IndForGeno.txt ' + chipFile + ' > ChosenInd.txt')  # only individuals chosen for genotypisation - ALL
    os.system("sed 's/^ *//' ChosenInd.txt > ChipFile.txt")  # Remove blank spaces at the beginning
    os.system("cut -f1 -d ' ' ChipFile.txt > Individuals.txt")  # obtain IDs
    os.system('''awk '{$1=""; print $0}' ChipFile.txt | sed 's/ //g' > Snps.txt''')  # obtain SNP genotypes
    os.system(
        r'''paste Individuals.txt Snps.txt | awk '{printf "%- 10s %+ 15s\n",$1,$2}' > GenoFile.txt''')  # obtain SNP genotypes of the last generation
pd.read_csv(AlphaSimDir + '/SimulatedData/Chip1SnpInformation.txt', sep='\s+')[[0, 1, 2]].to_csv(
    AlphaSimDir + 'SnpMap.txt', index=None, sep=" ", header=None)



