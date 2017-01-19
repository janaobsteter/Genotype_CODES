import matplotlib
from pyevolve import G1DList
from pyevolve import G2DList
from pyevolve import GSimpleGA
from pyevolve import Selectors
from pyevolve import Statistics
from pyevolve import DBAdapters
from pyevolve import Mutators

#from GenomeBase import GenomeBase

#class G1DBinaryString(GenomeBase):
 #  pass

chromosome = ([0,0,0,0,0,0])

def eval_func1(chromosome):
    score = 0.0
   # iterate over the chromosome
    for F_rod, F_pro, F_gen, M_rod, M_pro, M_gen in chromosome:
        sumF = chromosome[0] + chromosome[1] + chromosome[2]
        sumM = chromosome[3] + chromosome[4] + chromosome[5]
        scoreF = (F_rod /sumF)*0.4 + (F_pro /sumF)*0.75 + (F_gen /sumF)*0.8         
        scoreM = (M_rod /sumM)*0.4 + (M_pro /sumM)*0.75 + (M_gen /sumM)*0.8                        
        score = 0.8 * scoreM + 0.2 * scoreF    
    return score
   
def eval_func2(chromosome):
    score = 0.0
    for F_rod, F_pro, F_gen, M_rod, M_pro, M_gen in chromosome:
        scoreF = F_rod*10 + F_pro*6 + F_gen*1
        scoreM = M_rod*10 + M_pro*6 + M_gen*1
        score = 0.6*scoreM + 0.4*scoreF
    return score

genome = G2DList.G2DList(1,6)
genome.evaluator.set(eval_func1)
genome.evaluator.set(eval_func2)
ga = GSimpleGA.GSimpleGA(genome)
ga.evolve(freq_stats=10)
genome.mutator.set(Mutators.G1DListMutatorRealGaussian)
genome.mutator.set(Mutators.G1DListMutatorSwap)
print ga.bestIndividual()


