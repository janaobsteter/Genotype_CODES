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

os.chdir("/home/jana/Documents/CompBioGA")
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
import rpy2.robjects as ro
import scipy.stats

#
#def eval_func1(chromosome):
#    st_F = 3400
#    st_M = 1200
#    sum_f = sum(chromosome[:3])
#    sum_m = sum(chromosome[3:])
#    PA_f_tested = chromosome[0]/sum_f if chromosome[0] != 0 else 0
#    PT_f_tested = chromosome[1]/sum_f if chromosome[1] != 0 else 0
#    GT_f_tested = chromosome[2]/sum_f if chromosome[2] != 0 else 0
#    PA_f_sel = 140 / c 
#    PT_f_sel = 
#    GT_f_sel = 
#    iC_F = scipy.stats.norm.cdf(scipy.stats.norm.ppf(PT_f)) / PT_f
#    iG_F = scipy.stats.norm.cdf(scipy.stats.norm.ppf(GT_f)) / GT_f
#
#    deltaC_f = PA_f*3400 * 2.8 * 0.55  / 4.25
#    deltaG_f = 0.64 * 2.8 * 0.8 / 2 
#    deltaG1_m = 3.46 * 2.8 * 0.567  / 1.75
#    deltaG2_m = 3.46 *2.8 * 0.99  / 6
#    deltaG3_m = 3.46 * 2.8 * 0.8 / 1.75
#
#
#
#
#    score_f = PA_f*deltaG1_f + PT_f*deltaG2_f + GT_f*deltaG3_f
#    PA_m = chromosome[3]/sum_m if chromosome[3] != 0 else 0
#    PT_m = chromosome[4]/sum_m if chromosome[4] != 0 else 0
#    GT_m = chromosome[5]/sum_m if chromosome[5] != 0 else 0
#    score_m = PA_m*deltaG1_m + PT_m*deltaG2_m + GT_m*deltaG3_m
#    score1 = 0.5*score_m + 0.5*score_f
#    return score1

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

def eval_func1(chromosome):
    st_f = 3400
    st_m = 1200
    sum_f = sum(chromosome[:2])
    sum_m = sum(chromosome[2:])
    #percentages of categories
    #############FEMALES#########################
    #PA_f = (chromosome[0]/sum_f) if sum_f != 0 else 0
    PT_f = chromosome[0]/sum_f if sum_f != 0 else 0
    GT_f = chromosome[1]/sum_f if sum_f != 0 else 0
    #intensities for FEMALES - performance testing
    PT_fp = (140 / (st_f * PT_f)) if  PT_f !=0 else 0 #percentage of females selected (selected / tested)
    PT_fi = (scipy.stats.norm.cdf(scipy.stats.norm.ppf(PT_fp)) / PT_fp) if  PT_fp !=0 else 0
    #intensities for FEMALES - genomic testing
    GT_fp = 140 / (st_f * GT_f) if  GT_f !=0 else 0 #percentage of females selected (selected / tested)
    GT_fi = (scipy.stats.norm.cdf(scipy.stats.norm.ppf(GT_fp)) / GT_fp) if  GT_fp !=0 else 0
    
    #################MALES#############
    #PA_m = chromosome[3]/sum_m if sum_m != 0 else 0
    PT_m = chromosome[2]/sum_m if sum_m != 0 else 0
    GT_m = chromosome[3]/sum_m if sum_m != 0 else 0
    #intensities for MALES - performance testing
    PT_mp = 25 / (st_m * PT_m) if  PT_m !=0 else 0 #percentage of females selected (selected / tested)
    PT_mi = (scipy.stats.norm.cdf(scipy.stats.norm.ppf(PT_mp)) / PT_mp) if  PT_mp !=0 else 0
    #intensities for MALES - genomic testing
    GT_mp = 25 / (st_m * GT_m) if  GT_m !=0 else 0 #percentage of females selected (selected / tested)
    GT_mi = (scipy.stats.norm.cdf(scipy.stats.norm.ppf(GT_mp)) / GT_mp) if  GT_mp !=0 else 0

    ###############FINAL SCORE#######################
    score = (((0.55 * PT_fi) + (0.99 * PT_mi) + (0.8 * GT_fi) + (0.8 * GT_mi))*2.8) / (4.25 + 6 + 2 + 1.75)
    return score
    
    
# Genome instance
genome = G1DList.G1DList(4)
genome.setParams(rangemin=0, rangemax=100)

# The evaluator function (objective function)
genome.evaluator.set(eval_func1)
#genome.evaluator.add(eval_func2)
genome.mutator.set(Mutators.G1DListMutatorIntegerRange)


# Genetic Algorithm Instance
ga = GSimpleGA.GSimpleGA(genome)
ga.setGenerations(500)
ga.selector.set(Selectors.GTournamentSelector)
ga.setMutationRate(0.09)
ga.setCrossoverRate(0.05)
ga.setPopulationSize(100)
#ga.terminationCriteria.set(GSimpleGA.ConvergenceCriteria)


# Create DB Adapter and set as adapter
sqlite_adapter = DBAdapters.DBSQLite(identify="noPrice")
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