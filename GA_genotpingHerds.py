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


os.chdir("/home/jana/Documents/")
#evaluation functions for genetic algorithm
cowsGen = 5000

HerdsA = pd.read_csv('/home/jana/Documents/PhD/CompBio/HerdsA.txt')
HerdsAnim = pd.read_csv("/home/jana/Documents/PhD/CompBio/NoByHerd.txt")



def eval_func(chromosome):
    NoAnimals = sum([no for (chrom, no) in zip (chromosome, HerdsAnim.NoAnimal)])
    score = 0
    #within
    withinA = []
    for index, vals in HerdsA.iterrows():
       if (chromosome[int(vals.num1)-1] == 1 and chromosome[int(vals.num2)-1] == 1): 
           withinA.append(vals.rel)
    betweenA = []
    for index, vals in HerdsA.iterrows():
       if (chromosome[int(vals.num1)-1] == 1 and chromosome[int(vals.num2)-1] == 0) or (chromosome[int(vals.num1)-1] == 0 and chromosome[int(vals.num2)-1] == 1): 
           betweenA.append(vals.rel)
   
    #and also the number of animals 
    deviation = 1/(abs(NoAnimals - cowsGen))
    return (sum(betweenA) - sum(withinA)) * deviation * 1000

    
    
# Genome instance
genome = G1DList.G1DList(100)
genome.setParams(rangemin=0, rangemax=1)

# The evaluator function (objective function)
genome.evaluator.set(eval_func)
#genome.evaluator.add(eval_func2)
genome.mutator.set(Mutators.G1DListMutatorIntegerRange)


# Genetic Algorithm Instance
ga = GSimpleGA.GSimpleGA(genome)
ga.setGenerations(500)
ga.selector.set(Selectors.GTournamentSelector)
ga.setMutationRate(0.09)
ga.setCrossoverRate(0.05)
ga.setPopulationSize(100)
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
