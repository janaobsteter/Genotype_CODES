# -*- coding: utf-8 -*-
from __future__ import division
import os
import pandas as pd
import random
import sys
import math
from collections import defaultdict
import shutil
import numpy as np
import resource
from selection10 import blupf90

WorkingDir = "/home/jana/Documents/PhD/CompBio/TestingGBLUP/"
rounds = int(raw_input("Enter the number of repetitions"))
Accuracies = pd.DataFrame()
#to je skript, ki vozi GA v ponovitvah
os.chdir(WorkingDir)
for rep in range(rounds):
    #1) dobi rešitev iz GA
    os.system(" python ~/Genotipi/Genotipi_CODES/GA_genotpingHerds2.py > GAherds.txt")
    
    #ekstrahiraj rešitev
    chromosome = [int(x) for x in open(WorkingDir + "GAherds.txt").read().strip("\n")[open(WorkingDir + "GAherds.txt").read().strip("\n").find("List:"):].strip("'").strip("List:\t\t ").strip("[").strip("]").split(", ")]
    
    #ekstrahiraj živali
    ped = pd.read_csv("/home/jana/Documents/PhD/CompBio/PedigreeAndGeneticValues_Herds.txt", sep=",")
    pedO = pd.read_csv(WorkingDir + "PedigreeAndGeneticValues_cat.txt", sep="\s+")
    genK = [herd for (herd, gen) in zip(sorted(list(set(ped.Herd))), chromosome) if gen ==1]  
    pd.DataFrame({"ID": list(ped.loc[ped.Herd.isin(genK), 'Indiv']) + list(pedO.loc[pedO.cat.isin(["potomciNP", "pb"]),'Indiv']) }).to_csv(WorkingDir + '/IndForGeno.txt', index=None, header=None)
    #tukaj zapišeš IndForGeno.txt
    
    #Tukaj skreiraj GenoFile
    os.chdir(WorkingDir)
    os.system(
    'grep -Fwf IndForGeno.txt /home/jana/bin/AlphaSim1.05Linux/REALFillIn20BurnIn20/SimulatedData/AllIndividualsSnpChips/Chip1Genotype.txt > ChosenInd.txt')  # only individuals chosen for genotypisation - ALL
    os.system("sed 's/^ *//' ChosenInd.txt > ChipFile.txt")  # Remove blank spaces at the beginning
    os.system("cut -f1 -d ' ' ChipFile.txt > Individuals.txt")  # obtain IDs
    os.system('''awk '{$1=""; print $0}' ChipFile.txt | sed 's/ //g' > Snps.txt''')  # obtain SNP genotypes
    os.system(
    r'''paste Individuals.txt Snps.txt | awk '{printf "%- 10s %+ 15s\n",$1,$2}' > GenoFile.txt''')  # obtain SNP genotypes of the last generation
    pd.read_csv('/home/jana/bin/AlphaSim1.05Linux/SimulatedData/Chip1SnpInformation.txt', sep='\s+')[[0, 1, 2]].to_csv('SnpMap.txt', index=None, sep=" ", header=None)
    print "Created Geno File"
    
    #sfuraj blupf90
    os.system("./renumf90 < renumParam")  # run renumf90
    
    resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
    os.system('./blupf90 renf90.par')
    #renumber the solutions
    # copy the solution in a file that does not get overwritten
    os.system("bash Match_AFTERRenum.sh")
    
    
    #dodaj rešitve in izračunaj točnost
    blupSol = pd.read_csv('renumbered_Solutions', header=None,
                        sep='\s+', names=['renID', 'ID', 'Solution'])
    AlphaPed = pd.read_table("PedigreeAndGeneticValues_cat.txt", sep=" ")
    AlphaSelPed = AlphaPed.loc[:, ['Generation', 'Indiv', 'Father', 'Mother','cat', 'gvNormUnres1']]
    AlphaSelPed.loc[:, 'EBV'] = blupSol.Solution
    AlphaSelPed = AlphaSelPed.loc[AlphaSelPed.cat.isin(["potomciNP"])]
    Accuracies.loc[:,str(rep)] = list(np.corrcoef(AlphaSelPed.EBV, AlphaSelPed.gvNormUnres1)[0])
    AlphaSelPed.to_csv('GenPed_EBV' + str(rep) + '.txt', index=None)
    
Accuracies.to_csv(WorkingDir + "AccuraciesRep.txt")
    