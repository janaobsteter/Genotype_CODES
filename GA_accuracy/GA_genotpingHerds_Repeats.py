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
#from selection10 import blupf90

WorkingDir = "/home/jana/Documents/PhD/CompBio/TestingGBLUP/"
rounds = raw_input("Enter the number of repetitions")
Accuracies = pd.DataFrame(np.nan, index=range(rounds), columns=['Opt', 'Random', 'RandomHerd'])
#to je skript, ki vozi GA v ponovitvah


def reLu(number):
    return (0 if number < 0 else number)

for rep in range(rounds):
    #1) dobi rešitev iz GA
    os.makedirs(WorkingDir + "/Rep_" + str(rep))
    RepDir = WorkingDir + "/Rep_" + str(rep)
    os.chdir(RepDir)
    os.system("cp " + WorkingDir + "/Essentials/* .")
    os.system(" python ~/Genotipi/Genotipi_CODES/GA_genotpingHerds2.py > GAherds.txt")
    
    #ekstrahiraj rešitev
    chromosome = [int(x) for x in open(RepDir + "/GAherds.txt").read().strip("\n")[open(RepDir + "/GAherds.txt").read().strip("\n").find("List:"):].strip("'").strip("List:\t\t ").strip("[").strip("]").split(", ")]
    
    #ekstrahiraj živali
    ped = pd.read_csv("/home/jana/Documents/PhD/CompBio/TestingGBLUP/PedCows_HERDS_Total.txt", sep=" ")
    pedO = pd.read_csv(WorkingDir + "PedigreeAndGeneticValues_cat.txt", sep="\s+")
    genK = [herd for (herd, gen) in zip(sorted(list(set(ped.cluster))), chromosome) if gen ==1]  
    pd.DataFrame({"ID": list(ped.loc[ped.cluster.isin(genK), 'Indiv']) + list(pedO.loc[pedO.cat.isin(["potomciNP", "pb"]),'Indiv']) }).to_csv(RepDir + '/IndForGeno.txt', index=None, header=None)
    #tukaj zapišeš IndForGeno.txt
    
    #to je enako število random izbranih krav
    noCows = len(list(ped.loc[ped.cluster.isin(genK), 'Indiv']))
    pd.DataFrame({"ID": list(random.sample(ped.Indiv, noCows)) + list(pedO.loc[pedO.cat.isin(["potomciNP", "pb"]),'Indiv']) }).to_csv(RepDir + '/IndForGeno_Random.txt', index=None, header=None)
   
    #to je enako število random izbranih čred
    noHerds = sum(chromosome)
    randomHerds = sorted(random.sample(range(1, 101), noHerds))
    pd.DataFrame({"ID": list(ped.loc[ped.cluster.isin(randomHerds), 'Indiv']) + list(pedO.loc[pedO.cat.isin(["potomciNP", "pb"]),'Indiv']) }).to_csv(RepDir + '/IndForGeno_RandomHerds.txt', index=None, header=None)   
    
    #Tukaj skreiraj GenoFile
    os.system(
    'grep -Fwf IndForGeno.txt /home/jana/bin/AlphaSim1.05Linux/REALFillIn20BurnIn20/SimulatedData/AllIndividualsSnpChips/Chip1Genotype.txt > ChosenInd.txt')  # only individuals chosen for genotypisation - ALL
    os.system("sed 's/^ *//' ChosenInd.txt > ChipFile.txt")  # Remove blank spaces at the beginning
    os.system("cut -f1 -d ' ' ChipFile.txt > Individuals.txt")  # obtain IDs
    os.system('''awk '{$1=""; print $0}' ChipFile.txt | sed 's/ //g' > Snps.txt''')  # obtain SNP genotypes
    os.system(
    r'''paste Individuals.txt Snps.txt | awk '{printf "%- 10s %+ 15s\n",$1,$2}' > GenoFile.txt''')  # obtain SNP genotypes of the last generation
    pd.read_csv('/home/jana/bin/AlphaSim1.05Linux/REALFillIn20BurnIn20/SimulatedData/Chip1SnpInformation.txt', sep='\s+')[[0, 1, 2]].to_csv('SnpMap.txt', index=None, sep=" ", header=None)
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
    AlphaPed = pd.read_table(WorkingDir + "/PedigreeAndGeneticValues_cat.txt", sep=" ")
    AlphaSelPed = AlphaPed.loc[:, ['Generation', 'Indiv', 'Father', 'Mother','cat', 'gvNormUnres1']]
    AlphaSelPed.loc[:, 'EBV'] = blupSol.Solution
    AlphaSelPed = AlphaSelPed.loc[AlphaSelPed.cat.isin(["potomciNP"])]
    Accuracies.Opt[rep] = list(np.corrcoef(AlphaSelPed.EBV, AlphaSelPed.gvNormUnres1)[0])[1]
    AlphaSelPed.to_csv('GenPed_EBV' + str(rep) + '_Opt.txt', index=None)
  
    
      
    #potem pa naredi za vsako optimizacijo še eno random izbiro  
    #Tukaj skreiraj GenoFile
    os.system("rm GenoFile*")
    os.system(
    'grep -Fwf IndForGeno_Random.txt /home/jana/bin/AlphaSim1.05Linux/REALFillIn20BurnIn20/SimulatedData/AllIndividualsSnpChips/Chip1Genotype.txt > ChosenInd.txt')  # only individuals chosen for genotypisation - ALL
    os.system("sed 's/^ *//' ChosenInd.txt > ChipFile.txt")  # Remove blank spaces at the beginning
    os.system("cut -f1 -d ' ' ChipFile.txt > Individuals.txt")  # obtain IDs
    os.system('''awk '{$1=""; print $0}' ChipFile.txt | sed 's/ //g' > Snps.txt''')  # obtain SNP genotypes
    os.system(
    r'''paste Individuals.txt Snps.txt | awk '{printf "%- 10s %+ 15s\n",$1,$2}' > GenoFile.txt''')  # obtain SNP genotypes of the last generation
    pd.read_csv('/home/jana/bin/AlphaSim1.05Linux/REALFillIn20BurnIn20/SimulatedData/Chip1SnpInformation.txt', sep='\s+')[[0, 1, 2]].to_csv('SnpMap.txt', index=None, sep=" ", header=None)
    print "Created Geno File for Random choice"
    
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
    AlphaPed = pd.read_table(WorkingDir + "/PedigreeAndGeneticValues_cat.txt", sep=" ")
    AlphaSelPed = AlphaPed.loc[:, ['Generation', 'Indiv', 'Father', 'Mother','cat', 'gvNormUnres1']]
    AlphaSelPed.loc[:, 'EBV'] = blupSol.Solution
    AlphaSelPed = AlphaSelPed.loc[AlphaSelPed.cat.isin(["potomciNP"])]
    Accuracies.Random[rep] = list(np.corrcoef(AlphaSelPed.EBV, AlphaSelPed.gvNormUnres1)[0])[1]
    AlphaSelPed.to_csv('GenPed_EBV' + str(rep) + '_Random.txt', index=None)
    
    
    #potem pa naredi za vsako optimizacijo še eno random izbiro  ČRED
    #Tukaj skreiraj GenoFile
    os.system("rm GenoFile*")
    os.system(
    'grep -Fwf IndForGeno_RandomHerds.txt /home/jana/bin/AlphaSim1.05Linux/REALFillIn20BurnIn20/SimulatedData/AllIndividualsSnpChips/Chip1Genotype.txt > ChosenInd.txt')  # only individuals chosen for genotypisation - ALL
    os.system("sed 's/^ *//' ChosenInd.txt > ChipFile.txt")  # Remove blank spaces at the beginning
    os.system("cut -f1 -d ' ' ChipFile.txt > Individuals.txt")  # obtain IDs
    os.system('''awk '{$1=""; print $0}' ChipFile.txt | sed 's/ //g' > Snps.txt''')  # obtain SNP genotypes
    os.system(
    r'''paste Individuals.txt Snps.txt | awk '{printf "%- 10s %+ 15s\n",$1,$2}' > GenoFile.txt''')  # obtain SNP genotypes of the last generation
    pd.read_csv('/home/jana/bin/AlphaSim1.05Linux/REALFillIn20BurnIn20/SimulatedData/Chip1SnpInformation.txt', sep='\s+')[[0, 1, 2]].to_csv('SnpMap.txt', index=None, sep=" ", header=None)
    print "Created Geno File for Random HERD choice"
    
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
    AlphaPed = pd.read_table(WorkingDir + "/PedigreeAndGeneticValues_cat.txt", sep=" ")
    AlphaSelPed = AlphaPed.loc[:, ['Generation', 'Indiv', 'Father', 'Mother','cat', 'gvNormUnres1']]
    AlphaSelPed.loc[:, 'EBV'] = blupSol.Solution
    AlphaSelPed = AlphaSelPed.loc[AlphaSelPed.cat.isin(["potomciNP"])]
    Accuracies.RandomHerd[rep] = list(np.corrcoef(AlphaSelPed.EBV, AlphaSelPed.gvNormUnres1)[0])[1]
    AlphaSelPed.to_csv('GenPed_EBV' + str(rep) + '_RandomHerd.txt', index=None)

Accuracies.to_csv(WorkingDir + "AccuraciesRep.txt")


HerdsA = pd.read_csv('/home/jana/Documents/PhD/CompBio/RefADF_mean.csv')
NapA = pd.read_csv('/home/jana/Documents/PhD/CompBio/NapADF_mean.csv')
PbA = pd.read_csv('/home/jana/Documents/PhD/CompBio/PbADF_mean.csv')
HerdsAnim = pd.read_table("/home/jana/Documents/PhD/CompBio/HerdNo.txt", sep=" ")
WorkingDir = "/home/jana/Documents/PhD/CompBio/TestingGBLUP/"
cowsGen=5000
ped = pd.read_csv("/home/jana/Documents/PhD/CompBio/TestingGBLUP/PedCows_HERDS_Total.txt", sep=" ")



Relationship = pd.DataFrame(np.nan, index=range(rounds), columns=['Way', 'Rep', 'NoAnimals', 'NoHerds', 'Within', 'Between', 'Score', 'FinalScore'])
for rep in range(rounds):
    Relationship.Rep[rep] = rep
    Relationship.Way[rep] = "Opt"
    #1) dobi rešitev iz GA
    RepDir = WorkingDir + "/Rep_" + str(rep)
    os.chdir(RepDir)
    chromosome = [int(x) for x in open(RepDir + "/GAherds.txt").read().strip("\n")[open(RepDir + "/GAherds.txt").read().strip("\n").find("List:"):].strip("'").strip("List:\t\t ").strip("[").strip("]").split(", ")]
    
    NoAnimals = sum([no for (chrom, no) in zip (chromosome, HerdsAnim.NoAnim) if chrom == 1])
    chosenHerds = [herd for (chrom, herd) in zip (chromosome, HerdsAnim.Herd) if chrom == 1]
    Relationship.NoAnimals[rep] = NoAnimals
    Relationship.NoHerds[rep] = len(chosenHerds)
    
    withinA = []
    for index, vals in HerdsA.iterrows():
       if (int(vals.Herd1) in chosenHerds) and (int(vals.Herd2) in chosenHerds): 
           withinA.append(vals.A)
    
    withPb = (PbA.A[PbA.Herd.isin(chosenHerds)])
    withNap = (NapA.A[NapA.Herd.isin(chosenHerds)])
                  
    within = np.mean(list(withPb) + list(withinA))
    between = np.mean (withNap)
    
    Relationship.Within[rep] = within
    Relationship.Between[rep] = between

    #and also the number of animals 
    score = (reLu(between - within) * 10000) **2
    penalty = [-score if (NoAnimals > 1.5*cowsGen or NoAnimals < 0.85*cowsGen) else 0]
    
    Relationship.Score[rep] = score
    Relationship.FinalScore[rep] = score+penalty[0]
    RelationOpt = Relationship
    
Relationship = pd.DataFrame(np.nan, index=range(rounds), columns=['Way', 'Rep', 'NoAnimals', 'NoHerds', 'Within', 'Between', 'Score', 'FinalScore'])
for rep in range(rounds):
    Relationship.Rep[rep] = rep
    Relationship.Way[rep] = "RandomHerd"
    #1) dobi rešitev iz GA
    RepDir = WorkingDir + "/Rep_" + str(rep)
    os.chdir(RepDir)
    
    Inds = pd.read_table("IndForGeno_RandomHerds.txt", header=None)
    herds = sorted(list(set(ped.cluster[ped.Indiv.isin(list(Inds.loc[:, 0]))])))
    chromosome = [1 if herd in herds else 0 for herd in range(1, 101) ]

    NoAnimals = sum([no for (chrom, no) in zip (chromosome, HerdsAnim.NoAnim) if chrom == 1])
    chosenHerds = [herd for (chrom, herd) in zip (chromosome, HerdsAnim.Herd) if chrom == 1]
    Relationship.NoAnimals[rep] = NoAnimals
    Relationship.NoHerds[rep] = len(chosenHerds)
    
    withinA = []
    for index, vals in HerdsA.iterrows():
       if (int(vals.Herd1) in chosenHerds) and (int(vals.Herd2) in chosenHerds): 
           withinA.append(vals.A)
    
    withPb = (PbA.A[PbA.Herd.isin(chosenHerds)])
    withNap = (NapA.A[NapA.Herd.isin(chosenHerds)])
                  
    within = np.mean(list(withPb) + list(withinA))
    between = np.mean (withNap)
    
    Relationship.Within[rep] = within
    Relationship.Between[rep] = between

    #and also the number of animals 
    score = (reLu(between - within) * 10000) **2
    penalty = [-score if (NoAnimals > 1.5*cowsGen or NoAnimals < 0.85*cowsGen) else 0]
    
    Relationship.Score[rep] = score
    Relationship.FinalScore[rep] = score+penalty[0]
    
    RelationRandom = Relationship

Relationship.append(RelationOpt).to_csv(WorkingDir + "Relations.csv", index=None)