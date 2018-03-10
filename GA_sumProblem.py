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


#os.chdir("/home/jana/Documents/")
#evaluation functions for genetic algorithm






def eval_func(chromosome):
    if sum(chromosome) == 26235947428953663183191:
            score = 1000
    else:
            score = 0
    return score
    
    
# Genome instance
genome = G1DList.G1DList(89)
genome.setParams(rangemin=0, rangemax=26235947428953663183191)

# The evaluator function (objective function)
genome.evaluator.set(eval_func)
#genome.evaluator.add(eval_func2)
genome.mutator.set(Mutators.G1DListMutatorIntegerRange)


# Genetic Algorithm Instance
ga = GSimpleGA.GSimpleGA(genome)
ga.setGenerations(10000)
ga.selector.set(Selectors.GTournamentSelector)
ga.setMutationRate(0.09)
ga.setCrossoverRate(0.05)
ga.setPopulationSize(3000)
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
