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

def reLu(number):
    return (0 if number < 0 else number)
    

#os.chdir("/home/jana/Documents/")
#evaluation functions for genetic algorithm
cowsGen = 5000

HerdsA = pd.read_csv('/home/jana/Documents/PhD/CompBio/RefADF_mean.csv')
NapA = pd.read_csv('/home/jana/Documents/PhD/CompBio/NapADF_mean.csv')
PbA = pd.read_csv('/home/jana/Documents/PhD/CompBio/PbADF_mean.csv')
HerdsAnim = pd.read_table("/home/jana/Documents/PhD/CompBio/HerdNo.txt", sep=" ")



def eval_func(chromosome):
    NoAnimals = sum([no for (chrom, no) in zip (chromosome, HerdsAnim.NoAnim) if chrom == 1])
    chosenHerds = [herd for (chrom, herd) in zip (chromosome, HerdsAnim.Herd) if chrom == 1]
    
    withinA = []
    for index, vals in HerdsA.iterrows():
       if (int(vals.Herd1) in chosenHerds) and (int(vals.Herd2) in chosenHerds): 
           withinA.append(vals.A)
    
    withPb = (PbA.A[PbA.Herd.isin(chosenHerds)])
    withNap = (NapA.A[NapA.Herd.isin(chosenHerds)])
                  
    within = np.mean(list(withPb) + list(withinA))
    between = np.mean (withNap)

    #and also the number of animals 
    score = (reLu(between - within) * 10000) **2
    penalty = [-score if (NoAnimals > 1.5*cowsGen or NoAnimals < 0.85*cowsGen) else 0]
    return score+penalty[0]

    
    
# Genome instance
genome = G1DList.G1DList(100)
genome.setParams(rangemin=0, rangemax=1)

# The evaluator function (objective function)
genome.evaluator.set(eval_func)
#genome.evaluator.add(eval_func2)
genome.mutator.set(Mutators.G1DListMutatorIntegerRange)


# Genetic Algorithm Instance
ga = GSimpleGA.GSimpleGA(genome)
ga.setGenerations(900)
ga.selector.set(Selectors.GTournamentSelector)
ga.setMutationRate(0.01)
ga.setCrossoverRate(0.001)
ga.setPopulationSize(50)
ga.terminationCriteria.set(GSimpleGA.ConvergenceCriteria)


# Create DB Adapter and set as adapter
sqlite_adapter = DBAdapters.DBSQLite(identify="TestHerdsGen")
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
