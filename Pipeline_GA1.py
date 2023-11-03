# -*- coding: utf-8 -*-
from __future__ import division
import os
from selection10 import *
import shutil
import pandas as pd
from collections import defaultdict
import numpy as np
from pyevolve import G1DList, GSimpleGA, Selectors, Statistics
from pyevolve import Initializators, Mutators, Consts, DBAdapters
from math import log, log1p, exp
from pyevolve import G1DList
from pyevolve import GSimpleGA
from pyevolve import Selectors
from pyevolve import Statistics
from pyevolve import DBAdapters
import resource
import random
import sys

rep = sys.argv[1]
rounds = sys.argv[2]


homeDir = os.getcwd()


os.chdir(homeDir + "/FillInBurnIn" + str(rep) + "/")
WorkingDir = os.getcwd()

#ustvari poddirektorij
os.mkdir("GA/")
os.chdir("GA/")

GAdir = os.getcwd()

os.system("cp " + homeDir + "/CodeDir/GA/* .")
shutil.copy(WorkingDir + "/SimulatedData/PedigreeAndGeneticValues_cat.txt",
            GAdir + "/PedigreeAndGeneticValues_cat.txt")

#select individuals for optimization and create herds for cows
os.system("/exports/cmvm/eddie/eb/groups/tier2_hickey_external/R-3.4.2/bin/Rscript"
          " Choose_inds_create_herds.R > Choose_inds_create_herds.txt")

#calculate relationship and create H matrix for the selected individuals
#copy AlphaRelate to directory
#copy AlphaRelate_GA.txt to AlphaRelateSpec.txt
#prepare the spec file + prepare the genotype and pedigree files
shutil.copy(homeDir + "/CodeDir/GA/AlphaRelateSpec_GA.txt", GAdir + "/AlphaRelateSpec.txt")
pedA = AlphaRelate(GAdir, WorkingDir)
pedA.preparePedigree()

#run AlphaRelate
pedA.runAlphaRelate()

r"""
#Calculate relatedness according to herds
herds = pd.read_table("PedCows_HERDS.txt", sep=" ")
IndGeno = pd.read_table("INDPED.txt", header=None)

#Tukaj izračunaj sorodstvo med živalmi v obema čredama
RefAmean = defaultdict()

number = 1

for herd1 in range(1, 101):
    for herd2 in range(herd1, 101):
	print(str(herd1) + "_" + str(herd2))
        ref = sorted(list(herds.Indiv[herds.cluster.isin([herd1, herd2])]))  # tukaj odberi živali v obeh čredah

        pd.DataFrame({"ID": ref}).to_csv("IndMatrix.txt", index=None, header=None)

        os.system("grep -Fwf IndMatrix.txt PedigreeNrm.txt > RefMatrix")
        a = pd.read_table("RefMatrix", sep="\s+", header=None)
        a.columns = ["Indiv"] + list(IndGeno.loc[:, 0])

        refA = a.loc[:, ref]

        meanRef = np.mean(refA).mean()

        RefAmean[number] = [herd1, herd2, meanRef]
        number = number + 1

RefDF = pd.DataFrame.from_dict(RefAmean, orient="index")
RefADF = RefDF.drop_duplicates()
RefADF.columns = ["Herd1", "Herd2", "A"]
RefADF.to_csv("RefADF_mean.csv", index=None)

#tukaj izračunaj sorodstvo med živalmi v čredi in napovedno populacijo / plemenskimi biki (referenca)
ped = pd.read_table("PedigreeAndGeneticValues_cat.txt", sep=" ")
nr = ped.Indiv[ped.cat.isin(['potomciNP'])]
pb = ped.Indiv[ped.cat == 'pb']

NapAmean = defaultdict()
PbAmean = defaultdict()
number = 1
for herd in range(1, 101):
    # odberi živali v obeh čredah
    ref = sorted(list(herds.Indiv[herds.cluster == herd]))

    #naredi tabelo krav
    pd.DataFrame({"ID": ref}).to_csv("IndHerd.txt", index=None, header=None)

    os.system("grep -Fwf IndHerd.txt PedigreeNrm.txt  > HerdMatrix")
    a = pd.read_table("HerdMatrix", sep="\s+", header=None)
    a.columns = ["Indiv"] + list(IndGeno.loc[:, 0])

    refnapA = a.loc[:, list(nr)]  # sorodstvo z napovedno populacijo
    refpbA = a.loc[:, list(pb)]  # orodstvo s plemenskimi biki
    meanRefNap = np.mean(refnapA).mean()
    meanRefPb = np.mean(refpbA).mean()

    NapAmean[number] = [herd, meanRefNap]
    PbAmean[number] = [herd, meanRefPb]
    number = number + 1

NapADF = pd.DataFrame.from_dict(NapAmean, orient="index")
NapADF.columns = ["Herd", "A"]
NapADF.to_csv("NapADF_mean.csv", index=None)

PbADF = pd.DataFrame.from_dict(PbAmean, orient="index")
PbADF.columns = ["Herd", "A"]
PbADF.to_csv("PbADF_mean.csv", index=None)

################################
################################
#spusti GA
Accuracies = pd.DataFrame(np.nan, index=range(rounds), columns=['Opt', 'Random', 'RandomHerd'])


# to je skript, ki vozi GA v ponovitvah

def reLu(number):
    return (0 if number < 0 else number)


for rep in range(rounds):
    # 1) dobi rešitev iz GA
    os.makedirs(GAdir + "/Rep_" + str(rep))
    RepDir = GAdir + "/Rep_" + str(rep)
    os.chdir(RepDir)
    os.system("cp " + homeDir + "/Essentials/* .")
    os.system("python GA_genotpingHerds2.py > GAherds.txt")

    # ekstrahiraj rešitev
    chromosome = [int(x) for x in open("GAherds.txt").read().strip("\n")[
                                  open("GAherds.txt").read().strip("\n").find("List:"):].strip("'").strip(
        "List:\t\t ").strip("[").strip("]").split(", ")]

    # ekstrahiraj živali
    ped = pd.read_csv("PedCows_HERDS_Total.txt", sep=" ")
    pedO = pd.read_csv("PedigreeAndGeneticValues_cat.txt", sep="\s+")
    #tukaj vzami chromosome in vključi izbrane črede krav + pb + potomciNP
    genK = [herd for (herd, gen) in zip(sorted(list(set(ped.cluster))), chromosome) if gen == 1]
    pd.DataFrame({"ID": list(ped.loc[ped.cluster.isin(genK), 'Indiv']) + list(
        pedO.loc[pedO.cat.isin(["potomciNP", "pb"]), 'Indiv'])}).to_csv(RepDir + '/IndForGeno.txt', index=None,
                                                                        header=None)
    # tukaj zapišeš IndForGeno.txt

    #tukaj izberi naključne krave
    # to je enako število random izbranih krav
    noCows = len(list(ped.loc[ped.cluster.isin(genK), 'Indiv']))
    pd.DataFrame({"ID": list(random.sample(ped.Indiv, noCows)) + list(
        pedO.loc[pedO.cat.isin(["potomciNP", "pb"]), 'Indiv'])}).to_csv(RepDir + '/IndForGeno_Random.txt',
                                                                        index=None, header=None)

    #tu izberi naključne črede krav
    # to je enako število random izbranih čred
    noHerds = sum(chromosome)
    randomHerds = sorted(random.sample(range(1, 101), noHerds))
    pd.DataFrame({"ID": list(ped.loc[ped.cluster.isin(randomHerds), 'Indiv']) + list(
        pedO.loc[pedO.cat.isin(["potomciNP", "pb"]), 'Indiv'])}).to_csv(RepDir + '/IndForGeno_RandomHerds.txt',
                                                                        index=None, header=None)

    # Tukaj skreiraj GenoFile
    os.system(
        'grep -Fwf IndForGeno.txt ' + WorkingDir +
        '/SimulatedData/AllIndividualsSnpChips/Chip1Genotype.txt > ChosenInd.txt')
    os.system("sed 's/^ *//' ChosenInd.txt > ChipFile.txt")
    os.system("cut -f1 -d ' ' ChipFile.txt > Individuals.txt")
    os.system('''awk '{$1=""; print $0}' ChipFile.txt | sed 's/ //g' > Snps.txt''')
    os.system(
        r'''paste Individuals.txt Snps.txt | awk '{printf "%- 10s %+ 15s\n",$1,$2}' > GenoFile.txt''')
    pd.read_csv(WorkingDir + 'SimulatedData/Chip1SnpInformation.txt',
                sep='\s+')[[0, 1, 2]].to_csv('SnpMap.txt', index=None, sep=" ", header=None)
    print("Created Geno File")

    #vstavi ime za genotipsko datoteko
    os.system("sed 's/GENOTYPEFILE/GenoFile.txt/g' renumf90_generic.par > renumf90.par")

    # sfuraj blupf90
    os.system("./renumf90 < renumParam")  # run renumf90

    resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
    os.system('./blupf90 renf90.par')
    # renumber the solutions
    # copy the solution in a file that does not get overwritten
    os.system("bash Match_AFTERRenum.sh")

    # dodaj rešitve in izračunaj točnost
    blupSol = pd.read_csv('renumbered_Solutions', header=None,
                          sep='\s+', names=['renID', 'ID', 'Solution'])
    AlphaPed = pd.read_table(WorkingDir + "/PedigreeAndGeneticValues_cat.txt", sep=" ")
    AlphaSelPed = AlphaPed.loc[:, ['Generation', 'Indiv', 'Father', 'Mother', 'cat', 'gvNormUnres1']]
    AlphaSelPed.loc[:, 'EBV'] = blupSol.Solution
    AlphaSelPed = AlphaSelPed.loc[AlphaSelPed.cat.isin(["potomciNP"])]
    Accuracies.Opt[rep] = list(np.corrcoef(AlphaSelPed.EBV, AlphaSelPed.gvNormUnres1)[0])[1]
    AlphaSelPed.to_csv('GenPed_EBV' + str(rep) + '_Opt.txt', index=None)


    # potem pa naredi za vsako optimizacijo še eno random izbiro
    # Tukaj skreiraj GenoFile
    os.system("rm GenoFile*")
    os.system(
        'grep -Fwf IndForGeno_Random.txt ' + WorkingDir +
        '/SimulatedData/AllIndividualsSnpChips/Chip1Genotype.txt > ChosenIndRandom.txt')
    os.system("sed 's/^ *//' ChosenIndRandom.txt > ChipFileRandom.txt")
    os.system("cut -f1 -d ' ' ChipFileRandom.txt > IndividualsRandom.txt")
    os.system('''awk '{$1=""; print $0}' ChipFileRandom.txt | sed 's/ //g' > SnpsRandom.txt''')
    os.system(
        r'''paste IndividualsRandom.txt SnpsRandom.txt | awk '{printf "%- 10s %+ 15s\n",$1,$2}' 
        > GenoFileRandom.txt''')
    pd.read_csv(WorkingDir + '/SimulatedData/Chip1SnpInformation.txt',
                sep='\s+')[[0, 1, 2]].to_csv('SnpMap.txt', index=None, sep=" ", header=None)
    print("Created Geno File for Random choice")

    #vstavi ime za genotipsko datoteko
    os.system("sed 's/GENOTYPEFILE/GenoFileRandom.txt/g' renumf90_generic.par > renumf90.par")

    # sfuraj blupf90
    os.system("./renumf90 < renumParam")  # run renumf90

    resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
    os.system('./blupf90 renf90.par')
    # renumber the solutions
    # copy the solution in a file that does not get overwritten
    os.system("bash Match_AFTERRenum.sh")

    # dodaj rešitve in izračunaj točnost
    blupSol = pd.read_csv('renumbered_Solutions', header=None,
                          sep='\s+', names=['renID', 'ID', 'Solution'])
    AlphaPed = pd.read_table(WorkingDir + "/PedigreeAndGeneticValues_cat.txt", sep=" ")
    AlphaSelPed = AlphaPed.loc[:, ['Generation', 'Indiv', 'Father', 'Mother', 'cat', 'gvNormUnres1']]
    AlphaSelPed.loc[:, 'EBV'] = blupSol.Solution
    AlphaSelPed = AlphaSelPed.loc[AlphaSelPed.cat.isin(["potomciNP"])]
    Accuracies.Random[rep] = list(np.corrcoef(AlphaSelPed.EBV, AlphaSelPed.gvNormUnres1)[0])[1]
    AlphaSelPed.to_csv('GenPed_EBV' + str(rep) + '_Random.txt', index=None)

    # potem pa naredi za vsako optimizacijo še eno random izbiro  ČRED
    # Tukaj skreiraj GenoFile
    os.system("rm GenoFile*")
    os.system(
        'grep -Fwf IndForGeno_RandomHerds.txt ' + WorkingDir +
        '/SimulatedData/AllIndividualsSnpChips/Chip1Genotype.txt > ChosenIndRandomHerd.txt')  # only individuals chosen for genotypisation - ALL
    os.system("sed 's/^ *//' ChosenIndRandomHerd.txt > ChipFileRandomHerd.txt")  # Remove blank spaces at the beginning
    os.system("cut -f1 -d ' ' ChipFileRandomHerd.txt > IndividualsRandomHerd.txt")  # obtain IDs
    os.system('''awk '{$1=""; print $0}' ChipFileRandomHerd.txt | sed 's/ //g' > SnpsRandomHerd.txt''')  # obtain SNP genotypes
    os.system(
        r'''paste IndividualsRandomHerd.txt SnpsRandomHerd.txt | awk '{printf "%- 10s %+ 15s\n",$1,$2}' 
        > GenoFileRandomHerd.txt''')  # obtain SNP genotypes of the last generation
    pd.read_csv(WorkingDir + '/SimulatedData/Chip1SnpInformation.txt',
                sep='\s+')[[0, 1, 2]].to_csv('SnpMap.txt', index=None, sep=" ", header=None)
    print("Created Geno File for Random HERD choice")


    #vstavi ime za genotipsko datoteko
    os.system("sed 's/GENOTYPEFILE/GenoFileRandomHerd.txt/g' renumf90_generic.par > renumf90.par")


    # sfuraj blupf90
    os.system("./renumf90 < renumParam")  # run renumf90

    resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
    os.system('./blupf90 renf90.par')
    # renumber the solutions
    # copy the solution in a file that does not get overwritten
    os.system("bash Match_AFTERRenum.sh")

    # dodaj rešitve in izračunaj točnost
    blupSol = pd.read_csv('renumbered_Solutions', header=None,
                          sep='\s+', names=['renID', 'ID', 'Solution'])
    AlphaPed = pd.read_table("PedigreeAndGeneticValues_cat.txt", sep=" ")
    AlphaSelPed = AlphaPed.loc[:, ['Generation', 'Indiv', 'Father', 'Mother', 'cat', 'gvNormUnres1']]
    AlphaSelPed.loc[:, 'EBV'] = blupSol.Solution
    AlphaSelPed = AlphaSelPed.loc[AlphaSelPed.cat.isin(["potomciNP"])]
    Accuracies.RandomHerd[rep] = list(np.corrcoef(AlphaSelPed.EBV, AlphaSelPed.gvNormUnres1)[0])[1]
    AlphaSelPed.to_csv('GenPed_EBV' + str(rep) + '_RandomHerd.txt', index=None)

Accuracies.to_csv("AccuraciesRep.txt")



os.chdir(GAdir)
#tukaj pa sedaj naredi tabelo, kjer zbereš vse podatke:
#točnost, število živali, sorodnost, scoreGA ...
HerdsA = pd.read_csv('RefADF_mean.csv')
NapA = pd.read_csv('NapADF_mean.csv')
PbA = pd.read_csv('PbADF_mean.csv')
HerdsAnim = pd.read_csv("HerdNo.txt")
cowsGen = 5000
ped = pd.read_csv("PedCows_HERDS_Total.txt", sep=" ")

Relationship = pd.DataFrame(np.nan, index=range(rounds),
                            columns=['Way', 'Rep', 'NoAnimals', 'NoHerds', 'Within', 'Between', 'Score',
                                     'FinalScore'])
for rep in range(rounds):
    Relationship.Rep[rep] = rep
    Relationship.Way[rep] = "Opt"
    # 1) dobi rešitev iz GA
    RepDir = GAdir + "/Rep_" + str(rep)
    os.chdir(RepDir)
    chromosome = [int(x) for x in open(RepDir + "/GAherds.txt").read().strip("\n")[
                                  open(RepDir + "/GAherds.txt").read().strip("\n").find("List:"):].strip("'").strip(
        "List:\t\t ").strip("[").strip("]").split(", ")]

    NoAnimals = sum([no for (chrom, no) in zip(chromosome, HerdsAnim.NoAnim) if chrom == 1])
    chosenHerds = [herd for (chrom, herd) in zip(chromosome, HerdsAnim.Herd) if chrom == 1]
    Relationship.NoAnimals[rep] = NoAnimals
    Relationship.NoHerds[rep] = len(chosenHerds)

    withinA = []
    for index, vals in HerdsA.iterrows():
        if (int(vals.Herd1) in chosenHerds) and (int(vals.Herd2) in chosenHerds):
            withinA.append(vals.A)

    withPb = (PbA.A[PbA.Herd.isin(chosenHerds)])
    withNap = (NapA.A[NapA.Herd.isin(chosenHerds)])

    within = np.mean(list(withPb) + list(withinA))
    between = np.mean(withNap)

    Relationship.Within[rep] = within
    Relationship.Between[rep] = between

    # and also the number of animals
    score = (reLu(between - within) * 10000) ** 2
    penalty = [-score if (NoAnimals > 1.5 * cowsGen or NoAnimals < 0.85 * cowsGen) else 0]

    Relationship.Score[rep] = score
    Relationship.FinalScore[rep] = score + penalty[0]
    #tukaj naredi kopijo Relationship kot RelationShipOpt - in ponovi postopek za randomherd

RelationOpt = Relationship

Relationship = pd.DataFrame(np.nan, index=range(rounds),
                            columns=['Way', 'Rep', 'NoAnimals', 'NoHerds', 'Within', 'Between', 'Score',
                                     'FinalScore'])
for rep in range(rounds):
    Relationship.Rep[rep] = rep
    Relationship.Way[rep] = "RandomHerd"
    # 1) dobi rešitev iz GA
    RepDir = "Rep_" + str(rep)
    os.chdir(RepDir)

    Inds = pd.read_table("IndForGeno_RandomHerds.txt", header=None)
    herds = sorted(list(set(ped.cluster[ped.Indiv.isin(list(Inds.loc[:, 0]))])))
    chromosome = [1 if herd in herds else 0 for herd in range(1, 101)]

    NoAnimals = sum([no for (chrom, no) in zip(chromosome, HerdsAnim.NoAnim) if chrom == 1])
    chosenHerds = [herd for (chrom, herd) in zip(chromosome, HerdsAnim.Herd) if chrom == 1]
    Relationship.NoAnimals[rep] = NoAnimals
    Relationship.NoHerds[rep] = len(chosenHerds)

    withinA = []
    for index, vals in HerdsA.iterrows():
        if (int(vals.Herd1) in chosenHerds) and (int(vals.Herd2) in chosenHerds):
            withinA.append(vals.A)

    withPb = (PbA.A[PbA.Herd.isin(chosenHerds)])
    withNap = (NapA.A[NapA.Herd.isin(chosenHerds)])

    within = np.mean(list(withPb) + list(withinA))
    between = np.mean(withNap)

    Relationship.Within[rep] = within
    Relationship.Between[rep] = between

    # and also the number of animals
    score = (reLu(between - within) * 10000) ** 2
    penalty = [-score if (NoAnimals > 1.5 * cowsGen or NoAnimals < 0.85 * cowsGen) else 0]

    Relationship.Score[rep] = score
    Relationship.FinalScore[rep] = score + penalty[0]

    RelationRandom = Relationship

Relationship.append(RelationOpt).to_csv("Relations.csv", index=None)
"""

