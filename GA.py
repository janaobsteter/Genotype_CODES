# -*- coding: utf-8 -*-
from pyevolve import G1DList, GSimpleGA, Selectors, Statistics
from pyevolve import Initializators, Mutators, Consts, DBAdapters
import os
from math import log, log1p, exp

os.chdir("/home/jana/Documents/")
#evaluation functions for genetic algorithm

#evaluation function 1: genetic gain - MILK YIELD
#deltaG = accuracy * genetic variability * intensity / generation interval
#you have to have genetic gain for each specific category
#cat1: rodovniska napoved krave - 
    #accuracy = 0.2, intensity = 0.02, genetic variability 
#cat2: lastna preizkušnja krave
#cat3: genomska napoved krave
#cat4: rodovniška napoved biki
#cat5: progeni test biki
#cat6: genomska napoved biki

deltaG1_f = 0.2 * 0.32 * 0.05 * 10 / 2 
deltaG2_f = 0.6 * 0.32 * 0.05 * 50 / 4.25
deltaG3_f = 0.75 * 0.32 * 0.05 * 100 / 2 
deltaG1_m = 0.2 * 0.32 * 0.3 * 10 / 1.75
deltaG2_m = 0.99 * 0.32 * 0.3 * 50 / 6.5
deltaG3_m = 0.75 * 0.32 * 0.3 * 100/ 1.75

# This function is the evaluation function, we want
# to give high score to more zero'ed chromosomes
def eval_func1(chromosome):
    female3 = (100 - sum(chromosome[:2])) if sum(chromosome[:2]) < 100 else 0
    male3 = (100 - sum(chromosome[2:])) if sum(chromosome[2:]) < 100 else 0
    score_f = chromosome[0]*deltaG1_f + chromosome [1]*deltaG2_f + female3*deltaG3_f
    score_m = chromosome[2]*deltaG1_m + chromosome[3]*deltaG2_m + male3*deltaG3_m 
    score1 = 0.7*score_m + 0.3*score_f
    #print "score1 ", score1
    return exp(score1)
   
"""def eval_func2(chromosome):
    score = 0.0
    if (sum(chromosome[:3]) == 100.0) and (sum(chromosome[3:]) == 100.0):
        score += 10000.0
    return score
"""

def eval_func2(chromosome):
    female3 = (100 - sum(chromosome[:2])) if sum(chromosome[:2]) < 100 else 0
    male3 = (100 - sum(chromosome[2:])) if sum(chromosome[2:]) < 100 else 0
    score2 = chromosome[0] * 0.5 + chromosome[1] * 0.3 + female3 * 0.08 + \
    chromosome[2] * 0.5 + chromosome[3] * 0.3 + male3 * 0.08
    return log1p(score2)
    
    
# Genome instance
genome = G1DList.G1DList(4)
genome.setParams(rangemin=0.0, rangemax=100.0)

# The evaluator function (objective function)
genome.evaluator.set(eval_func1)
#genome.evaluator.add(eval_func2)
genome.mutator.set(Mutators.G1DListMutatorRealRange)


# Genetic Algorithm Instance
ga = GSimpleGA.GSimpleGA(genome)
ga.setGenerations(1000)
ga.selector.set(Selectors.GTournamentSelector)
ga.setMutationRate(0.09)
ga.setCrossoverRate(0.2)
ga.setPopulationSize(500)

# Create DB Adapter and set as adapter
sqlite_adapter = DBAdapters.DBSQLite(identify="two")
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