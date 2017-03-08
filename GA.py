# -*- coding: utf-8 -*-
from pyevolve import G1DList, GSimpleGA, Selectors, Statistics
from pyevolve import Initializators, Mutators, Consts, DBAdapters
import os
from math import log, log1p, exp
from pyevolve import G1DList
from pyevolve import GSimpleGA
from pyevolve import Selectors
from pyevolve import Statistics
from pyevolve import DBAdapters

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

deltaG1_f = 0.64 * 2.8 * 0.567 / 2 
deltaG2_f = 0.64 * 2.8 * 0.55  / 4.25
deltaG3_f = 0.64 * 2.8 * 0.8 / 2 
deltaG1_m = 3.46 * 2.8 * 0.567  / 1.75
deltaG2_m = 3.46 *2.8 * 0.99  / 6
deltaG3_m = 3.46 * 2.8 * 0.8 / 1.75

#price included
#deltaG1_f = 0.64 * 2.8 * 0.567 * 10 / 2 
#deltaG2_f = 0.64 * 2.8 * 0.55 * 50 / 4.25
#deltaG3_f = 0.64 * 2.8 * 0.8 * 100 / 2 
#deltaG1_m = 1.16 * 2.8 * 0.567 * 10 / 1.75
#deltaG2_m = 1.16 *2.8 * 0.99 * 50 / 6.5
#deltaG3_m = 1.16 * 2.8 * 0.8 * 100/ 1.75

# This function is the evaluation function, we want
# to give high score to more zero'ed chromosomes
#def eval_func1(chromosome):
#    female3 = (100 - sum(chromosome[:2])) if sum(chromosome[:2]) < 100 else 0
#    male3 = (100 - sum(chromosome[2:])) if sum(chromosome[2:]) < 100 else 0
#    score_f = chromosome[0]*deltaG1_f + chromosome [1]*deltaG2_f + female3*deltaG3_f
#    score_m = chromosome[2]*deltaG1_m + chromosome[3]*deltaG2_m + male3*deltaG3_m 
#    score1 = 0.5*score_m + 0.5*score_f
#    #print "score1 ", score1
#    return score1
#return exp(score1)

#
def eval_func1(chromosome):
    sum_f = sum(chromosome[:3])
    sum_m = sum(chromosome[3:])
    PA_f = chromosome[0]/sum_f if chromosome[0] != 0 else 0
    PT_f = chromosome[1]/sum_f if chromosome[1] != 0 else 0
    GT_f = chromosome[2]/sum_f if chromosome[2] != 0 else 0
    score_f = PA_f*deltaG1_f + PT_f*deltaG2_f + GT_f*deltaG3_f
    PA_m = chromosome[3]/sum_m if chromosome[3] != 0 else 0
    PT_m = chromosome[4]/sum_m if chromosome[4] != 0 else 0
    GT_m = chromosome[5]/sum_m if chromosome[5] != 0 else 0
    score_m = PA_m*deltaG1_m + PT_m*deltaG2_m + GT_m*deltaG3_m
    score1 = 0.5*score_m + 0.5*score_f
    return score1

#with price
#def eval_func1(chromosome):
#    sum_f = sum(chromosome[:3])
#    sum_m = sum(chromosome[3:])
#    PA_f = chromosome[0]/sum_f if chromosome[0] != 0 else 0
#    PT_f = chromosome[1]/sum_f if chromosome[1] != 0 else 0
#    GT_f = chromosome[2]/sum_f if chromosome[2] != 0 else 0
#    score_f = PA_f*deltaG1_f*4 + PT_f*deltaG2_f*3 + GT_f*deltaG3_f*2
#    PA_m = chromosome[3]/sum_m if chromosome[3] != 0 else 0
#    PT_m = chromosome[4]/sum_m if chromosome[4] != 0 else 0
#    GT_m = chromosome[5]/sum_m if chromosome[5] != 0 else 0
#    score_m = PA_m*deltaG1_m*4 + PT_m*deltaG2_m*3 + GT_m*deltaG3_m*2
#    score1 = 0.5*score_m + 0.5*score_f
#    return score1

    
   
"""def eval_func2(chromosome):
    score = 0.0
    if (sum(chromosome[:3]) == 100.0) and (sum(chromosome[3:]) == 100.0):
        score += 10000.0
    return score
"""

def eval_func2(chromosome):
    sum_f = sum(chromosome[:3])
    sum_m = sum(chromosome[3:])
    PA_f = chromosome[0]/sum_f if chromosome[0] != 0 else 0
    PT_f = chromosome[1]/sum_f if chromosome[1] != 0 else 0
    GT_f = chromosome[2]/sum_f if chromosome[2] != 0 else 0
    score_f = PA_f*5 + PT_f*3 + GT_f*1
    PA_m = chromosome[3]/sum_m if chromosome[3] != 0 else 0
    PT_m = chromosome[4]/sum_m if chromosome[4] != 0 else 0
    GT_m = chromosome[5]/sum_m if chromosome[5] != 0 else 0
    score_m = PA_m*5 + PT_m*3 + GT_m*1
    scoreP = 0.5*score_m + 0.5*score_f
    return scoreP
    
    
# Genome instance
genome = G1DList.G1DList(6)
genome.setParams(rangemin=0.0, rangemax=100)

# The evaluator function (objective function)
genome.evaluator.set(eval_func1)
#genome.evaluator.add(eval_func2)
genome.mutator.set(Mutators.G1DListMutatorIntegerRange)


# Genetic Algorithm Instance
ga = GSimpleGA.GSimpleGA(genome)
#ga.setGenerations(500)
ga.selector.set(Selectors.GTournamentSelector)
ga.setMutationRate(0.09)
ga.setCrossoverRate(0.05)
ga.setPopulationSize(5000)
ga.terminationCriteria.set(GSimpleGA.ConvergenceCriteria)


# Create DB Adapter and set as adapter
sqlite_adapter = DBAdapters.DBSQLite(identify="price")
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