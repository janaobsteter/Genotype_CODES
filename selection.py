from __future__ import division
import pandas as pd
import numpy as np

class pedigree:
    def __init__(self, pedfile):
        self.ped = pd.read_csv(pedfile)
        self.gen = set(self.ped.Generation)
        self.ngen = len(set(self.ped.Generation))
        self.ped['sex'] = ""
        self.ped['cat'] = ""
        self.ped['active'] = ""
    
    def head(self):
        return self.ped.head()
        
        
    def set_cat(self, start, stop, cat): #pregled po kategorijah
        self.ped.loc[range(start,(start + stop)), 'cat'] = cat
        
    def set_sex(self, start, stop,sex): #pregled po kategorijah
        self.ped.loc[range(start,(start + stop)), 'sex'] = sex
        
    def set_active(self, start, stop,active): #pregled po kategorijah
        self.ped.loc[range(start,(start + stop)), 'active'] = active 
    
    def cat(self):
        return self.ped.cat.value_counts()
    
    def sex(self):
        return self.ped.sex.value_counts()
        
    def active(self):
        return self.ped.active.value_counts()
    
    def row_gen(self, gen):
        return self.ped.ix[self.ped.Generation == gen].index
    
    def row_gen_F(self, gen):
        return self.ped.ix[self.ped.Generation == gen & self.ped.sex=="F"].index
        
    def row_newest_gen(self):
        return self.ped.ix[self.ped.Generation == max(self.ped.Generation)].index
        
    def row_oldest_gen(self):
        return self.ped.ix[self.ped.Generation == min(self.ped.Generation)].index
        
a = pedigree("~/Documents/PhD/Simulaton/Pedigrees/PedPython.txt")