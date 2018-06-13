# -*- coding: utf-8 -*-
from __future__ import division
from pyevolve import G1DList, GSimpleGA, Selectors, Statistics

from pyevolve import Initializators, Mutators, Consts, DBAdapters
import os
from math import log, log1p, exp
from pyevolve import G1DList
from pyevolve import GSimpleGA
from pyevolve import Selectors
from pyevolve import Statistics
from pyevolve import DBAdapters
import pandas as pd
import numpy as np

#define reLu function
def reLu(number):
    return (0 if number < 0 else number)
    

#os.chdir("/home/jana/Documents/")
#evaluation functions for genetic algorithm
cowsGen = 5000

#the file with information about relatedness between herds
HerdsA = pd.read_csv('/home/jana/Documents/PhD/CompBio/RefADF_mean.csv')
#the file with information about relatedness of herds with test population
NapA = pd.read_csv('/home/jana/Documents/PhD/CompBio/NapADF_mean.csv')
##the file with the relatedness of test population with bulls for reference (training population)
PbA = pd.read_csv('/home/jana/Documents/PhD/CompBio/PbADF_mean.csv')
#the file containing the number of animals in each herd
HerdsAnim = pd.read_table("/home/jana/Documents/PhD/CompBio/HerdNo.txt", sep=" ")


#fitness function
def eval_func(chromosome):
    #how many animals from the herds are chosen
    NoAnimals = sum([no for (chrom, no) in zip (chromosome, HerdsAnim.NoAnim) if chrom == 1])
    #which are the chosen herds (herds 1 - 100)
    chosenHerds = [herd for (chrom, herd) in zip (chromosome, HerdsAnim.Herd) if chrom == 1]
    
    #calculate the relatedness between the animals within the chosen cow herds for reference population (training population) [list]
    withinA = []
    for index, vals in HerdsA.iterrows():
       if (int(vals.Herd1) in chosenHerds) and (int(vals.Herd2) in chosenHerds): 
           withinA.append(vals.A)
    #the relatedness of the animals in chosen cow herds with reference bulls [list]
    withPb = (PbA.A[PbA.Herd.isin(chosenHerds)])

    #the relatedness of the chosen cows with animals in the testing population
    withNap = (NapA.A[NapA.Herd.isin(chosenHerds)])
                  
    #mean relatedness within the reference (training) population
    within = np.mean(list(withPb) + list(withinA))
    #mean relatedness between the reference and testing population
    between = np.mean (withNap)

    #compute the score considering the relatedness and also the number of animals in the reference
    score = (reLu(between - within) * 10000) **2
    penalty = [-score if (NoAnimals > 1.5*cowsGen or NoAnimals < 0.85*cowsGen) else 0]
    return score+penalty[0]

    
    
# Genome instance
genome = G1DList.G1DList(100) #chromosome is a list with 100 elements (one for each herds)
genome.setParams(rangemin=0, rangemax=1) #allowed values are 0 and 1

# The evaluator function (objective function)
genome.evaluator.set(eval_func) #evaluate chromosomes with the chosen fitness function
#genome.evaluator.add(eval_func2)
genome.mutator.set(Mutators.G1DListMutatorIntegerBinary) #mutate 0 and 1


# Genetic Algorithm Instance
ga = GSimpleGA.GSimpleGA(genome)
ga.setGenerations(900) #set the number of generations
ga.selector.set(Selectors.GTournamentSelector) #set the rule for parent selection
ga.setMutationRate(0.01) #set mutation rate
ga.setCrossoverRate(0.001) #set cross-over rate
ga.setPopulationSize(50) #set population size
#ga.terminationCriteria.set(GSimpleGA.ConvergenceCriteria)


# Create DB Adapter and set as adapter
sqlite_adapter = DBAdapters.DBSQLite(identify="TestHerdsGen") #this is to write the evolution to a file
ga.setDBAdapter(sqlite_adapter)

# Using CSV Adapter
#csvfile_adapter = DBAdapters.DBFileCSV()
#ga.setDBAdapter(csvfile_adapter)

# Using the URL Post Adapter
# urlpost_adapter = DBAdapters.DBURLPost(url="http://whatismyip.oceanus.ro/server_variables.php", post=False)
# ga.setDBAdapter(urlpost_adapter)

# Do the evolution, with stats dump
# frequency of 10 generations
ga.evolve(freq_stats=10)

# Best individual
print ga.bestIndividual()
