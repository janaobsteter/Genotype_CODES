# -*- coding: utf-8 -*-
from __future__ import division
import os
import pandas as pd
import numpy as np
from collections import defaultdict

class pedigree:
    def __init__(self, pedfile):
        self.ped = pd.read_csv(pedfile)
        self.gen = set(self.ped.Generation)
        self.ngen = len(set(self.ped.Generation))
        self.ped['sex'] = ""
        self.ped['cat'] = ""
        self.ped['active'] = ""
        self.nrow = len(self.ped)
        
    
    def head(self):
        return self.ped.head()
        
    def tail(self):
        return self.ped.tail()
        
    def displayR(self, start, stop):
        return self.ped.ix[start:stop]
    
    def rows(self):
        return len(self.ped)
    
    def view_cat(self, cat):
        return self.ped[self.ped.cat==cat]
    
    def compute_age (self):
        self.ped['age'] = max(self.ped['age']) - self.ped['age']
        
    def set_cat(self, start, stop, cat): #pregled po kategorijah
        self.ped.loc[range(start,(start + stop)), 'cat'] = cat
        
    
    def set_cat_list(self, seznam, cat):
        self.ped.loc[seznam, 'cat'] = cat

    
    
    def set_cat_gen(self, gen, cat):
        self.ped.loc[self.ped.Generation == gen, 'cat'] = cat

    def set_cat_age(self, age,  cat, prevGenDict):
        self.ped.loc[(self.ped['age'] ==age) & (self.ped.Indiv.isin(prevGenDict[cat])), 'cat'] = cat
                       
    def set_cat_old(self, oldcat, cat, prevGenDict): #določi novo kategorijo glede na prejšnjo kategorijo
        self.ped.loc[self.ped.Indiv.isin(prevGenDict[oldcat]), 'cat'] = cat                           

    def set_cat_age_old(self, age,  oldcat, cat, prevGenDict):
        self.ped.loc[(self.ped['age'] ==age) & (self.ped.Indiv.isin(prevGenDict[oldcat])), 'cat'] = cat
            
    def set_sex(self, start, stop,sex): #pregled po kategorijah
        self.ped.loc[range(start, stop), 'sex'] = sex
    
    def set_sex_list(self, seznam,sex): #pregled po kategorijah
        self.ped.loc[seznam, 'sex'] = sex    
                
    def set_active(self, start, stop,active): #pregled po kategorijah
        self.ped.loc[range(start,(start + stop)), 'active'] = active 
    
    def set_active_list(self, seznam,active): #pregled po kategorijah
        self.ped.loc[seznam, 'active'] = active 
    
    def cat(self):
        return self.ped.cat.value_counts()
    
    def sex(self):
        return self.ped.sex.value_counts()
  
    def age(self):
        return self.ped.age.value_counts()

    def gens(self):
        return set(self.ped['Generation'])
               
    def active(self):
        return self.ped.active.value_counts()

    
    def row_gen(self, gen):
        return list(self.ped.ix[self.ped.Generation == gen].index)

    def row_gen_subset(self, gen, start, stop):
        return list(self.ped.ix[self.ped.Generation == gen][start:stop].index)
    
    def row_gen_F(self, gen):
        return list(self.ped.ix[self.ped.Generation == gen & self.ped.sex=="F"].index)
        
    def row_newest_gen(self):
        return list(self.ped.ix[self.ped.Generation == max(self.ped.Generation)].index)
        
    def row_oldest_gen(self):
        return list(self.ped.ix[self.ped.Generation == min(self.ped.Generation)].index)
        
    def izloci_poEBV(self, sex, stIzl, oldcat, prevGenDict):
        izlRow = list(self.ped.loc[(self.ped.sex==sex) & (self.ped.Indiv.isin(prevGenDict[oldcat])), 'EBV'].sort_values(ascending=True)[:stIzl].index) #katere izločiš
        if len(izlRow) < stIzl:
            print("Too little animals too choose from.")
        else:
            self.set_cat_list(izlRow, "izl")
            self.set_active_list(izlRow, 2)
            
    def izloci_poEBV_age(self, sex,age, stIzl, oldcat, prevGenDict):
        izlRow = list(self.ped.loc[(self.ped.age==age)&(self.ped.sex==sex) & (self.ped.Indiv.isin(prevGenDict[oldcat])), 'EBV'].sort_values(ascending=True)[:stIzl].index) #katere izločiš
        if len(izlRow) < stIzl:
            print("Too little animals too choose from.")
        else:
            self.set_cat_list(izlRow, "izl")
            self.set_active_list(izlRow, 2)
                   
    def izberi_poEBV_top(self, sex, st, oldcat, cat, prevGenDict):
        selRow = list(self.ped.loc[(self.ped.sex==sex)  & (self.ped.Indiv.isin(prevGenDict[oldcat])) , 'EBV'].sort_values(ascending=False)[:st].index) #katere izbereš
        if len(selRow) < st:
            print("Too little animals too choose from.")
        else:
            self.set_cat_list(selRow, cat)
            self.set_active_list(selRow, 1)
            
    def print_sth(self):
        return 23

    def izberi_poEBV_top_age(self, sex, age, st, oldcat, cat, prevGenDict):
        selRow = list(self.ped.loc[(self.ped.age==age) & (self.ped.sex==sex)  & (self.ped.Indiv.isin(prevGenDict[oldcat])) , 'EBV'].sort_values(ascending=False)[:st].index) #katere izbereš
        if len(selRow) < st:
            print("Too little animals too choose from.")
        else:
            self.set_cat_list(selRow, cat)
            self.set_active_list(selRow, 1)
    
    def ind_poEBV_top(self, sex, st, cat):
        selRow = list(self.ped.loc[(self.ped.sex==sex)  & (self.ped.Indiv.isin(prevGenDict[oldcat])) , 'EBV'].sort_values(ascending=False)[:st].index) #katere izbereš
        if len(selRow) < st:
            print("Too little animals too choose from.")
        else:
            self.set_cat_list(selRow, cat)
            self.set_active_list(selRow, 1)
            
    def izberi_poEBV_OdDo(self, gen, sex, start, stop,  oldcat,cat, prevGenDict):
        selRow = list(self.ped.loc[(self.ped.sex==sex) & (self.ped.Generation==gen) & (self.ped.Indiv.isin(prevGenDict[oldcat])) , 'EBV'].sort_values(ascending=False)[start:stop].index) #katere izbereš
        if len(selRow) < (stop - start):
            print("Too little animals too choose from.")
        else:
            self.set_cat_list(selRow, cat)
            self.set_active_list(selRow, 1)
        
    def izberi_random(self, gen, sex, st, oldcat, cat, prevGenDict):
        freeRow = list(self.ped.loc[(self.ped.sex==sex) & (self.ped.Generation==gen) & (self.ped.Indiv.isin(prevGenDict[oldcat])) & (self.ped.cat == "")].index) #rows with no assigned category
        if len(freeRow) < st:
            print("Too little animals to choose from")
        else:
            choice = list(np.random.choice(freeRow, st, replace=False))
            self.set_cat_list(choice, cat)
            self.set_active_list(choice, 1)
     
    def izloci_random(self, gen, sex, st, oldcat, prevGenDict):
        freeRow = list(self.ped.loc[(self.ped.sex==sex) & (self.ped.Generation==gen)  & (self.ped.Indiv.isin(prevGenDict[oldcat])) & (self.ped.cat == "")].index ) #rows with no assigned category
        if len(freeRow) < st:
            print("Too little animals to choose from")
        else:
            choice = list(np.random.choice(freeRow, st, replace=False))
            self.set_cat_list(choice, 'izl')
            self.set_active_list(choice, 2)   
    
    def izloci_cat(self, cat, prevGenDict):
        self.ped.loc[self.ped.Indiv.isin(prevGenDict[cat]), 'cat'] = 'izl'
        self.ped.loc[self.ped.Indiv.isin(prevGenDict[cat]), 'active'] = 2
 
    def izloci_age_cat(self, age, cat, prevGenDict):
        self.ped.loc[(self.ped['age'] ==age) & (self.ped.Indiv.isin(prevGenDict[cat])), 'cat'] = 'izl'
        self.ped.loc[(self.ped['age'] ==age) & (self.ped.Indiv.isin(prevGenDict[cat])), 'active'] = 2      
  
    def save_cat(self):
        cat = defaultdict(list)
        for i in set(self.ped['cat']):
            cat[i] = list(self.ped.loc[self.ped.cat==i, 'Indiv'])
        return cat        
        
    def save_active(self):
        active = defaultdict(list)
        for i in set(self.ped['active']):
            active[i] = list(self.ped.loc[self.ped.active==i, 'Indiv'])
        return active

    def save_sex(self):
        sex = defaultdict(list)
        for i in set(self.ped['sex']):
            sex[i] = list(self.ped.loc[self.ped.sex==i, 'Indiv'])
        return sex
    
    def save_age(self):
        age = defaultdict(list)
        self.ped['age'] = self.ped['Generation'] - max(self.ped['Generation'])
        for i in set(self.ped['age']):
            age[i] = list(self.ped.loc[self.ped.age==i, 'Indiv'])
        return age
         
    def add_new_gen_blank(self, stF, stM):
        nr = pd.DataFrame({'Indiv': range((max(self.ped.Indiv)+1), (max(self.ped.Indiv) + 1 + stM + stF)), \
        'cat': ['nr']*(stF + stM), 'sex': (['F']*stF + ['M'] * stM), \
        'Generation' : [(max(self.ped.Generation) + 1)]*(stF + stM), \
        'EBV' : list(np.random.uniform(0.004, 1.5, (stF + stM))), \
        'active' : [1]*(stM + stF)}, \
        index=range(max(self.ped.index)+1, max(self.ped.index)+1 + stM + stF ))
        self.ped = pd.concat([self.ped,nr])
        