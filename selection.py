# -*- coding: utf-8 -*-
from __future__ import division
import os
import pandas as pd
import numpy as np
from collections import defaultdict
import random
import sys
import shutil
from itertools import chain
from subprocess import call

class pedigree:
    def __init__(self, pedfile):
        self.ped = pd.read_csv(pedfile)
        self.gen = set(self.ped.Generation)
        self.ngen = len(set(self.ped.Generation))
        self.ped['sex'] = ""
        self.ped['cat'] = ""
        self.ped['active'] = ""
        self.nrow = len(self.ped)
    
    def printDate(self):
        print 20042017 
    
    def head(self):
        return self.ped.head()
        
    def tail(self):
        return self.ped.tail()
        
    def displayR(self, start, stop):
        return self.ped.ix[start:stop]
    
    def displayCat (self, cat):
        return self.ped[self.ped.cat == cat]
    
    def rows(self):
        return len(self.ped)
    
    def view_cat(self, cat):
        return self.ped[self.ped.cat==cat]
    
    def gens(self):
        return set(self.ped['Generation'])
    
    def compute_age (self):
        self.ped['age'] = max(self.gens()) - self.ped.Generation 

        
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
   
    def set_cat_sex_old(self, sex,  oldcat, cat, prevGenDict):
        self.ped.loc[(self.ped['sex'] ==sex) & (self.ped.Indiv.isin(prevGenDict[oldcat])), 'cat'] = cat
   
    def set_cat_mother_catCurrent(self, mothercat,  cat):
        self.ped.loc[(self.ped.Mother.isin(self.catCurrent_indiv(mothercat))),  'cat'] = cat   
        
    def set_sex(self, start, stop,sex): #pregled po kategorijah
        self.ped.loc[range(start, stop), 'sex'] = sex
    
    def set_sex_list(self, seznam,sex): 
        self.ped.loc[seznam, 'sex'] = sex    

    def set_sex_prevGen(self, prevGenDict): 
        self.ped.loc[self.ped.Indiv.isin(prevGenDict['F']), 'sex'] = "F"
        self.ped.loc[self.ped.Indiv.isin(prevGenDict['M']), 'sex'] = "M"  
        
                                                
    def set_active(self, start, stop,active): #pregled po kategorijah
        self.ped.loc[range(start,(start + stop)), 'active'] = active 
    
    def set_active_list(self, seznam,active): #pregled po kategorijah
        self.ped.loc[seznam, 'active'] = active 
    
    def set_active_prevGen(self, prevGenDict): 
        self.ped.loc[self.ped.Indiv.isin(prevGenDict[1]), 'active'] = 1
        self.ped.loc[self.ped.Indiv.isin(prevGenDict[2]), 'active'] = 2 
    
    def set_active_age_cat(self, cat, age, active, prevGenDict):
        self.ped.loc[(self.ped['age'] ==age) & (self.ped.Indiv.isin(prevGenDict[cat])), 'active'] = 2 
    
    def set_active_cat(self, cat, active, prevGenDict):
        self.ped.loc[(self.ped.Indiv.isin(prevGenDict[cat])), 'active'] = 2 
    
    def cat(self):
        return self.ped.cat.value_counts()
    
    def sex(self):
        return self.ped.sex.value_counts()
  
    def age(self):
        return self.ped.age.value_counts()

    def mother_blank(self):
        return sum(self.ped.Mother == 0)

    def father_blank(self):
        return sum(self.ped.Father == 0)
        
    def mother_nr_blank(self):
        return len(self.ped[(self.ped.cat=='nr') & (self.ped.Mother == 0)])

    def father_nr_blank(self):
        return len(self.ped[(self.ped.cat=='nr') & (self.ped.Father == 0)])
               
    def active(self):
        return self.ped.active.value_counts()

    def cat_age(self, cat):
        return self.ped.loc[self.ped.cat==cat, 'age'].value_counts()
     
    def catCurrent_indiv(self, cat):
        return list(self.ped.loc[self.ped.cat==cat, 'Indiv'])
    
    def catCurrent_indiv_age(self, cat, age):
        return list(self.ped.loc[(self.ped.cat==cat) & (self.ped.age==age), 'Indiv'])
    
    def cat_active(self, cat):
        return self.ped.loc[self.ped.cat==cat, 'active'].value_counts()
           
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
  
    def row_cat(self, cat):
        return list(self.ped.ix[self.ped.cat == cat].index)
        
        
    def izloci_poEBV(self, sex, stIzl, oldcat, prevGenDict):
        izlRow = list(self.ped.loc[(self.ped.sex==sex) & (self.ped.Indiv.isin(prevGenDict[oldcat])), 'EBV'].sort_values(ascending=True)[:stIzl].index) #katere izločiš
        if len(izlRow) < stIzl:
            print("Too little animals to choose from. <{} {}>".format("izloci po EBV", oldcat))
        else:
            self.set_cat_list(izlRow, "izl")
            self.set_active_list(izlRow, 2)
            
    def izloci_poEBV_age(self, sex,age, stIzl, oldcat, prevGenDict):
        izlRow = list(self.ped.loc[(self.ped.age==age)&(self.ped.sex==sex) & (self.ped.Indiv.isin(prevGenDict[oldcat])), 'EBV'].sort_values(ascending=True)[:stIzl].index) #katere izločiš
        if len(izlRow) < stIzl:
            print("Too little animals to choose from. <{} {}>".format("izloci po EBV_age", oldcat))
        else:
            self.set_cat_list(izlRow, "izl")
            self.set_active_list(izlRow, 2)
                   
    def izberi_poEBV_top(self, sex, st, oldcat, cat, prevGenDict):
        selRow = list(self.ped.loc[(self.ped.sex==sex)  & (self.ped.Indiv.isin(prevGenDict[oldcat])) , 'EBV'].sort_values(ascending=False)[:st].index) #katere izbereš
        if len(selRow) < st:
            print("Too little animals to choose from. <{} {} > {}>".format("izberi po EBV", oldcat, cat))
        else:
            self.set_cat_list(selRow, cat)
            self.set_active_list(selRow, 1)
 
    def izberi_poEBV_top_catCurrent(self, sex, st, cat, newcat): #to je bolj kot ne samo za prvi krog, ko nimaš slovarja od prejšnje generacije
        selRow = list(self.ped.loc[(self.ped.sex==sex)  & (self.ped.cat == cat) , 'EBV'].sort_values(ascending=False)[:st].index) #katere izbereš
        if len(selRow) < st:
            print("Too little animals to choose from. <{} {} > {}>".format("izberi po EBV_catCurrent", cat, newcat))
        else:
            self.set_cat_list(selRow, newcat)
            self.set_active_list(selRow, 1)
            
            
    def izberi_poEBV_top_age_naive(self, age, st, newcat): #to je bolj kot ne samo za prvi krog, ko nimaš slovarja od prejšnje generacije
        selRow = list(self.ped.loc[(self.ped.age==age) & (self.ped.cat == ""), 'EBV'].sort_values(ascending=False)[:st].index) #katere izbereš
        if len(selRow) < st:
            print("Too little animals to choose from. <{} {} > {}>".format("izberi po EBV_catCurrent", newcat))
        else:
            self.set_cat_list(selRow, newcat)
            self.set_active_list(selRow, 1)

    def izberi_poEBV_top_age(self, sex, age, st, oldcat, cat, prevGenDict):
        selRow = list(self.ped.loc[(self.ped.age==age) & (self.ped.sex==sex)  & (self.ped.Indiv.isin(prevGenDict[oldcat])) , 'EBV'].sort_values(ascending=False)[:st].index) #katere izbereš
        if len(selRow) < st:
            print("Too little animals to choose from. <{} {} > {}>".format("izberi po EBV_topAge", oldcat, cat))
        else:
            self.set_cat_list(selRow, cat)
            self.set_active_list(selRow, 1)
    
    def ind_poEBV_top(self, sex, st, cat, prevGenDict):
        selRow = list(self.ped.loc[(self.ped.sex==sex)  & (self.ped.Indiv.isin(prevGenDict[cat])) , 'EBV'].sort_values(ascending=False)[:st].index) #katere izbereš
        if len(selRow) < st:
            print("Too little animals to choose from. <{} {}>".format("Indiv po EBV_top", cat))
        else:
            self.set_cat_list(selRow, cat)
            self.set_active_list(selRow, 1)
            
    def izberi_poEBV_OdDo(self, sex, start, stop,  oldcat,cat, prevGenDict):
        selRow = list(self.ped.loc[(self.ped.sex==sex) &  (self.ped.Indiv.isin(prevGenDict[oldcat])) , 'EBV'].sort_values(ascending=False)[start:stop].index) #katere izbereš
        if len(selRow) < (stop - start):
            print("Too little animals to choose from. <{} {} > {}>".format("izberi po EBV_OdDo", oldcat, cat))
        else:
            self.set_cat_list(selRow, cat)
            self.set_active_list(selRow, 1)
    
    def izberi_poEBV_OdDo_age_naive(self, age, start, stop,  cat):
        selRow = list(self.ped.loc[(self.ped.age==age) & (self.ped.cat == ""), 'EBV'].sort_values(ascending=False)[start:stop].index) #katere izbereš
        if len(selRow) < (stop - start):
            print("Too little animals to choose from. <{} {} > {}>".format("izberi po EBV_OdDo", cat))
        else:
            self.set_cat_list(selRow, cat)
            self.set_active_list(selRow, 1)    
            
                    
    def izberi_random(self,  sex, st, oldcat, cat, prevGenDict):
        freeRow = list(self.ped.loc[(self.ped.sex==sex) & (self.ped.Indiv.isin(prevGenDict[oldcat])) & (self.ped.cat == "")].index) #rows with no assigned category
        if len(freeRow) < st:
            print("Too little animals to choose from. <{} {} > {}>".format("izberi random", oldcat, cat))
        else:
            choice = list(np.random.choice(freeRow, st, replace=False))
            self.set_cat_list(choice, cat)
            self.set_active_list(choice, 1)
    
    def izberi_random_age_naive(self,  age, st, cat):
        freeRow = list(self.ped.loc[(self.ped.age==age)  & (self.ped.cat == "")].index) #rows with no assigned category
        if len(freeRow) < st:
            print("Too little animals to choose from. <{} {} > {}>".format("izberi random", cat))
        else:
            choice = list(np.random.choice(freeRow, st, replace=False))
            self.set_cat_list(choice, cat)
            self.set_active_list(choice, 1) 
      
        
    def izloci_random(self,  sex, st, oldcat, prevGenDict):
        freeRow = list(self.ped.loc[(self.ped.sex==sex) & (self.ped.Indiv.isin(prevGenDict[oldcat])) & (self.ped.cat == "")].index ) #rows with no assigned category
        if len(freeRow) < st:
            print("Too little animals to choose from. <{} {}>".format("izloči random", oldcat))
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
         
    def add_new_gen_naive(self, stNR, potomciNPn):
        nr = pd.DataFrame({'Indiv': range((max(self.ped.Indiv)+1), (max(self.ped.Indiv) + 1 + stNR)), \
        'cat': ['nr']*(stNR - potomciNPn) + ['potomciNP']*potomciNPn, 'sex': ['F', 'M'] * int(stNR / 2), \
        'Generation' : [(max(self.ped.Generation) + 1)]*stNR, \
        'EBV' : list(np.random.uniform(0.004, 1.5, stNR)), \
        'Mother' : [0]*stNR,
        'Father' : [0]*stNR,
        'active' : [1]*stNR}, \
        index=range(max(self.ped.index)+1, max(self.ped.index)+1 + stNR ))
        self.ped = pd.concat([self.ped,nr])
    
    def set_mother(self,  MotherList):
        first_blank_row = min(self.ped[(self.ped['Mother'] == 0) & ((self.ped.cat=='nr')) ].index)
        self.ped.loc[(first_blank_row):(first_blank_row + len(MotherList)-1), 'Mother'] = MotherList
    
    def set_mother_catPotomca(self,  MotherList, cat_potomca):
        freeRow = list(self.ped[(self.ped['Mother'] == 0) & (self.ped.cat==cat_potomca) ].index)[0:len(MotherList)]
        self.ped.loc[freeRow, 'Mother'] = MotherList
     
    def set_father(self,  FatherList, mothercat, prevGenDict):
        mother_cat_rows = list(self.ped[(self.ped.cat == 'nr') & (self.ped.Mother.isin(self.catCurrent_indiv(mothercat))) ].index)[0:len(FatherList)]
        self.ped.loc[mother_cat_rows, 'Father'] = FatherList
    
    def set_father_catPotomca(self,  FatherList, cat_potomca):
        freeRow = list(self.ped[(self.ped['Father'] == 0) & (self.ped.cat==cat_potomca) ].index)[0:len(FatherList)]
        self.ped.loc[freeRow, 'Father'] = FatherList       
        
    def select_mother_random(self, cat, st):
        pot_Mother = list(self.ped.loc[(self.ped.cat==cat) & (self.ped.active==1), 'Indiv'])
        return list(random.sample(pot_Mother, st))
        
    def select_mother_EBV_top(self, cat, st): #tukaj je trenutna categorija, saj jih doloačs koec generacije
        selRow = list(self.ped.loc[(self.ped.cat==cat) , 'EBV'].sort_values(ascending=False)[:st].index) #katere izbereš
        return list(self.ped.loc[selRow, 'Indiv'])
        
    def write_ped(self, path):
        pedNB = self.ped[self.ped.Generation == max(self.gens())]
        pedNB.to_csv(path, columns = ['Indiv', 'Father', 'Mother'], quoting=None, index=False, header=False)
        
    def select_age_0_1(self, categories, nrFn, nrMn, telFn, vhlevljenin, potomciNPn, telMn ): #tukaj odbereš iz novorojenih živali tel, ptel in mlade bike, pripust1
    #FEMALES
        self.set_cat_sex_old("F", "potomciNP", "telF", categories)
        izlF = nrFn - telFn#koliko jih izločiš
        self.izloci_poEBV("F", izlF, "nr", categories) #tukaj jih izloči, funkcija v modulu
    
        self.izberi_poEBV_top("F", (nrFn - izlF), "nr", "telF", categories) #izberi telice, ki jih osemeniš --> krave
        
        
        #MALES
        self.izberi_poEBV_top( "M", vhlevljenin, "potomciNP", "vhlevljeni", categories) #odberi mlade TO SAMO NA ZAČETKU; POTEM POTOMCI BM IN ELITE!
        self.izloci_poEBV("M", int(potomciNPn - vhlevljenin), 'potomciNP', categories)
        self.izberi_random("M", telMn, "nr", "telM", categories)
        self.izloci_random("M", int(nrMn - telMn), "nr", categories)

    def select_age_1_2(self, categories, ptn, mladin, vhlevljenin, bik12n): # tukaj odbereš nič pri kravah - razen, če so že bikovske matere, pripust 2, bike24
        #FEMALES
        self.izberi_poEBV_top("F", ptn, 'telF', 'pt', categories)
        self.izloci_poEBV("F", (len(categories['telF']) - ptn), 'telF', categories) #terlice postanejo
    
        
        #MALES
        self.izberi_poEBV_top( "M", mladin, "vhlevljeni", "mladi", categories) #odberi mlade
        self.izberi_poEBV_OdDo( "M", mladin, vhlevljenin, "vhlevljeni", "pripust1", categories) #odberi v pripustu
        self.izberi_random( "M", bik12n, 'telM', 'bik12', categories)
        self.izloci_random( "M", (len(categories['telM']) - bik12n), 'telM', categories)


    def select_age_2_3(self, categories, kraveUp, bmOdbira, bmn, bmUp, cak, pripust1n, pripust2n, mladin):
        #FEMALES
        #najprej dodaj nove krave
        self.set_cat_old('pt', 'k', categories) #osemenjene telice postanejo krave - predpostavimo, da vse
        #potem izloči najstarejše krave - po 4. laktaciji
        if ('k' in categories.keys()) and ((kraveUp+2) in self.age()): #izloči koliko laktacij + 2 leti
            self.izloci_age_cat((kraveUp+2), 'k', categories)
        #ostale krave prestavi naprej v krave - OZIROMA PODALJŠAJ STATUS!
        self.set_cat_age_old(3, 'k', 'k', categories)
        self.set_cat_age_old(4, 'k', 'k', categories) 
        self.set_cat_age_old(5, 'k', 'k', categories)    
                
        #če imaš že dovolj stare krave, potem odberi BM
        #BM se odbira po drugi laktaciji - to je starost 3 - 4 (starost v pedigreju = 3, ker imaš tudi 0)
        if ('k' in categories.keys()) and ((1 + bmOdbira) in self.age()):
            self.izberi_poEBV_top_age("F",3, int(bmn /bmUp), 'k', 'pBM', categories) #izberi bikovske matere
        #in izloči najastarejše BM, če jih imaš
        if ('bm' in categories.keys()):
            self.izloci_cat('bm', categories)
        #ostale BM prestavi naprej
        for i in range((1 + bmOdbira + 1), (1 + bmOdbira + bmUp)): #1 leto prva osemenitev, bm odbrane po 2. laktaciji, +1 da začneš prestavljat
            self.set_cat_age_old(i, 'pBM', 'pBM', categories)
        self.set_cat_age_old((1 + bmOdbira + bmUp), 'pBM', 'bm', categories) #spremeni kategorijo iz plemenskih BM v bm v zadnji laktaciji
        
        #MALES
        #mladi biki postanejo čakajoči (~1 leto, da se osemeni krave s semenom oz. malo po 2. letu)
        self.set_cat_old('mladi', 'cak', categories) 
        self.set_active_cat('mladi', 2, categories)
        
        #čakajočim bikov podaljšaj status (do starosti 5 let)
        #hkrati jim tudi nastavi status izl
        #ped.set_cat_age_old(2, 'cak', 'cak', categories)
        for i in range((2 + 1), (2 + cak)): #1 leto, ko začnejo semenit in so mladi biki, 3 so čakajoči, +1 da začneš prestavlajt
            self.set_cat_age_old(i, 'cak', 'cak', categories)
        
        #povprečna doba v pripustu - glede na to odberi bike, ki preživijo še eno leto
        if 'pripust1' in categories.keys():
            self.izberi_random( "M", pripust2n, 'pripust1', 'pripust2', categories)
            self.izloci_random( "M", (pripust1n - pripust2n), 'pripust1', categories)
    
        #plemenske bike prestavljaj naprej
        self.set_cat_old('pb', 'pb', categories)
        self.izloci_cat('bik12', categories)
        self.izloci_cat('pripust2', categories)
        if ('cak' in categories.keys()) and ((cak+2) in self.age()): #+2 - eno leto so teleta, eno leto mladi biki
            self.izberi_poEBV_top_age("M", (cak +2), int(mladin * 0.5), 'cak', 'pb', categories)
            self.set_active_cat('cak', 2, categories) #tukaj moraš to nastaviti, zato ker fja izberi avtomatsko nastavi na active=1
            self.izloci_poEBV_age("M",(cak+2), int(mladin * 0.5), 'cak', categories) #TUKAJ MORA BITI ŠE STAROST!!!!!!!!!!!
    

    def doloci_matere(self, stNB, ptn, kraveUp):
        #MATERE
        sTbmMother = 90 if len(self.catCurrent_indiv('pBM')) >= 90 else len(self.catCurrent_indiv('pBM'))
        if sTbmMother != 0:
            bmMother = self.select_mother_random('pBM', sTbmMother)
            self.set_mother_catPotomca(bmMother, 'potomciNP')
        #
        
        if 'k' in self.cat():#TUKAJ SO DOLOČENE SEDAJ VSE MATERE!!!
            mother = self.select_mother_EBV_top('k', int(round(ptn*kraveUp*0.7))) #tukaj odberi brez tistih, ki so za gospodarsko križanje
            if len(mother) >= (stNB - sTbmMother): # če že imaš dovolj krav, določi matere vsem novorojenim oz. odbiraš matere, saj jih imaš preveč!
                motherOther = random.sample(mother, (stNB - sTbmMother))
                self.set_mother_catPotomca(motherOther, 'nr') #TUKAJ SO DOLOČENE SEDAJ VSE MATERE!!!
            elif len(mother) < (stNB - sTbmMother): # če jih še ni dovolj, ne odbiraš mater, ampak uporabiš vse MINUS gosp. križanmja
                self.set_mother_catPotomca(mother, 'nr') 

    def doloci_ocete(self, stNB, potomciNPn, cak, pbUp, pripustDoz, mladiDoz, pozitivnoTestDoz):
    #OČETJE
        mladiOce = self.catCurrent_indiv('mladi')
        pripustOce = self.catCurrent_indiv('pripust1') + self.catCurrent_indiv('pripust2') 
        testiraniOce = list(chain.from_iterable([self.catCurrent_indiv_age('pb', (2 + cak + x)) for x in range(1, pbUp+1)])) # v času, ko določaš potomce, so že eno leto starjši!!!
        bmMother = 90 if len(self.catCurrent_indiv('pBM')) >= 90 else len(self.catCurrent_indiv('pBM'))
        if 'pb' in self.cat():
            elita = np.random.choice(self.catCurrent_indiv('pb'), bmMother, replace=True) #navidezna elita
    #        pd.Series(elita).value_counts()#preveri zastopanost po bikih
            #naštimaj očete elite --> BM
            self.set_father_catPotomca(elita, 'potomciNP')    
    
        ocetje = pripustOce*pripustDoz + testiraniOce*pozitivnoTestDoz + mladiOce*mladiDoz
        if len(ocetje) >= (stNB - potomciNPn*2): #če imaš dovolj DOZ za vse NB
            ocetjeNB = random.sample(ocetje, (stNB - potomciNPn*2)) #tukaj izbereš očete za vse krave  - razen BM!
            self.set_father_catPotomca(ocetjeNB, 'nr')
        if len(ocetje) < (stNB - potomciNPn*2):
            self.set_father_catPotomca(ocetje, 'nr')

    
def selekcija_ena_gen(pedFile, categories = None, sex = None, active = None, stNB=None, nrFn =None, \
nrMn=None, telFn=None, telMn=None, potomciNPn=None, vhlevljenin=None, ptn=None, mladin=None, bik12n=None, \
pripust1n=None, pripust2n=None, cak=None, kraveUp=None, bmOdbira=None, bmn=None, bmUp=None, pripustDoz=None, mladiDoz=None, \
pozitivnoTestDoz=None, pbUp=None):
    ped = pedigree(pedFile) 
    
    if max(ped.gen) == 1:
        ped.set_cat_gen(max(ped.gen), "nr")  # to je samo na prvem loopu
        ped.set_sex_list([x for x in range(0, ped.rows()) if x % 2 == 0], "F")
        ped.set_sex_list([x for x in range(0, ped.rows()) if x % 2 != 0], "M")
        ped.izberi_poEBV_top_catCurrent("F", int(potomciNPn), 'nr', 'potomciNP')
        ped.izberi_poEBV_top_catCurrent("M", int(potomciNPn), 'nr', 'potomciNP')
    
        #global categories #to moraš dat global samo v prvenm loopu, drugje dobiš return
        categories = ped.save_cat()
        #global sex
        sex = ped.save_sex()
        
        ped = pedigree(pedFile)
        
        ped.set_sex_prevGen(sex)# prva odbira
        ped.compute_age()
        ped.select_age_0_1(categories, nrFn, nrMn, telFn, vhlevljenin, potomciNPn, telMn)
        ped.add_new_gen_naive(stNB, potomciNPn*2)
        
        ped.compute_age()
        #dodaj matere
        ped.doloci_matere(stNB, ptn, kraveUp)
        #preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
        ped.mother_nr_blank()
        #dodaj očete
        ped.doloci_ocete(stNB, potomciNPn, cak, pripustDoz, mladiDoz, pozitivnoTestDoz)
        #preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
        ped.mother_nr_blank()   
    
        categories.clear()
        
        
    if max(ped.gen) == 2:
        # SETSEX!!!
        ped.set_sex_prevGen(sex)    
        ped.set_active_prevGen(active) 
                                
        # druga odbira
        ped.set_cat_gen(1, "")
        ped.set_cat_gen(2, "")
        ped.set_cat_old('izl', 'izl', categories)

        ped.compute_age()
        ped.select_age_0_1(categories, nrFn, nrMn, telFn, vhlevljenin, potomciNPn, telMn )
        ped.select_age_1_2(categories, ptn, mladin, vhlevljenin, bik12n)        
        ped.add_new_gen_naive(stNB, potomciNPn*2)
        ped.compute_age()
        
        #dodaj matere
        ped.doloci_matere(stNB, ptn, kraveUp)
        #preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
        ped.mother_nr_blank()
        #dodaj očete
        ped.doloci_ocete(stNB, potomciNPn, cak, pripustDoz, mladiDoz, pozitivnoTestDoz)
        #preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
        ped.mother_nr_blank()
                
        categories.clear() #sprazni slovar od prejšnjega leta
        
#dodaj starše novorojenim - VEDNO PRVO MAME, KER JE FUNKCIJA ZA OČETE NAŠTIMANA, DA SE RAVNA PO MAMAH!
    
    
    if max(ped.gen) >= 3:
        ped.set_sex_prevGen(sex)
        ped.set_active_prevGen(active)
        
        
        for i in ped.gens():
            ped.set_cat_gen(i, "")

        ped.set_cat_old('izl', 'izl', categories)

        ped.compute_age()
        ped.select_age_0_1(categories, nrFn, nrMn, telFn, vhlevljenin, potomciNPn, telMn )
        ped.select_age_1_2(categories, ptn, mladin, vhlevljenin, bik12n)        
        ped.select_age_2_3(categories, kraveUp, bmOdbira, bmUp, cak, pripust1n, pripust2n, mladin)
        ped.add_new_gen_naive(stNB, potomciNPn*2)  
        ped.compute_age()
        
        #dodaj matere
        ped.doloci_matere(stNB, ptn, kraveUp)
        #preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
        ped.mother_nr_blank()
        #dodaj očete
        ped.doloci_ocete(stNB, potomciNPn, cak, pripustDoz, mladiDoz, pozitivnoTestDoz)
        #preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
        ped.mother_nr_blank()
                
        categories.clear() #sprazni slovar od prejšnjega leta
        ped.write_ped("/home/jana/bin/AlphaSim1.05Linux/ExternalPedigree.txt")

    return ped, ped.save_cat(), ped.save_sex(), ped.save_active()




def nastavi_cat (PedFile, categories = None, sex = None, active = None, stNB=None, nrFn =None, \
nrMn=None, telFn=None, telMn=None, potomciNPn=None, vhlevljenin=None, ptn=None, bmn=None, mladin=None, bik12n=None, \
pripust1n=None, pripust2n=None,cak=None, kraveUp=None, bmOdbira=None, bmUp=None, pripustDoz=None, \
mladiDoz=None, pozitivnoTestDoz=None, pbUp=None):
    ped = pedigree(PedFile)
    ped.compute_age()
    
    
    #MALES FIRST
    #age 0
    #določi vhlevljene
    ped.izberi_poEBV_top_age_naive(0, vhlevljenin, 'vhlevljeni')
    #določi moška teleta pod 12
    ped.izberi_random_age_naive(0, telMn, 'telM')
    
    #age1
    #določi mlade
    ped.izberi_poEBV_top_age_naive(1, mladin, 'mladi')
    #določi pripust - 1. leto
    ped.izberi_poEBV_OdDo_age_naive(1, mladin, vhlevljenin, 'pripust1')
    #določi bike nad 12 m
    ped.izberi_random_age_naive(1, bik12n, 'bik12')
    
    
    #age2
    ped.izberi_poEBV_top_age_naive(2, mladin, 'cak')
    ped.izberi_poEBV_OdDo_age_naive(1, mladin, (mladin + pripust2n), 'pripust2')
    
    #age3,4
    for i in [3,4]:
        ped.izberi_poEBV_top_age_naive(i, mladin, 'cak')
    
    #age 5 - 10: pb
    pbAge = range((2 + cak), (2 + cak + pbUp)) if (2 + cak + pbUp) <= max(ped.gens()) else range((2 + cak), max(ped.gens()))
    for i in pbAge:
        ped.izberi_poEBV_top_age_naive(i, 4, 'pb')
    
    
    #FEMALES
    #age 0
    #določi ženska teleta pod 12
    ped.izberi_poEBV_top_age_naive(0, telFn, 'telF')
    
    #age1
    #določi plemenske telice
    ped.izberi_poEBV_top_age_naive(1, ptn, 'pt')
    
    #age2
    for i in range(2, (1 + bmOdbira)):
        ped.izberi_poEBV_top_age_naive(i, ptn, 'k')
    
    #age3,4,5
    #odberi plemenske bm najprej
    for i in range((1 + bmOdbira), (1 + bmOdbira + bmUp)):
        ped.izberi_poEBV_top_age_naive(i, int(bmn / bmUp), 'pBM')
        ped.izberi_poEBV_top_age_naive(i, (ptn - int(bmn / bmUp)), 'k')
    
    #age 6
    #izberi odslužene bm
    ped.izberi_poEBV_top_age_naive((1 + bmOdbira + bmUp), int(bmn / bmUp), 'bm')
    
    
    #ostali so izločeni
    #določi spol ženskim živalim
    ped.set_sex_list(ped.row_cat('telF'), "F")
    ped.set_sex_list(ped.row_cat('pt'), "F")
    ped.set_sex_list(ped.row_cat('k'), "F")
    ped.set_sex_list(ped.row_cat('pBM'), "F")
    ped.set_sex_list(ped.row_cat('bm'), "F")
    
    
    #določi spol moškim živalim
    ped.set_sex_list(ped.row_cat('vhlevljeni'), "M")
    ped.set_sex_list(ped.row_cat('telM'), "M")
    ped.set_sex_list(ped.row_cat('bik12'), "M")
    ped.set_sex_list(ped.row_cat('mladi'), "M")
    ped.set_sex_list(ped.row_cat('cak'), "M")
    ped.set_sex_list(ped.row_cat('pb'), "M")
    ped.set_sex_list(ped.row_cat('pripust1'), "M")
    ped.set_sex_list(ped.row_cat('pripust2'), "M")
    
    #določi še izločene
    ped.set_cat_list(ped.row_cat(""), 'izl')
    
    ped.add_new_gen_naive(stNB, potomciNPn*2)
    
    ped.compute_age()
    #dodaj matere
    ped.doloci_matere(stNB, ptn, kraveUp)
    #preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
    ped.mother_nr_blank()
    #dodaj očete
    ped.doloci_ocete(stNB, potomciNPn, cak, pripustDoz, mladiDoz, pozitivnoTestDoz)
 
    
    ped.write_ped("/home/jana/bin/AlphaSim1.05Linux/ExternalPedigree.txt")
    
    return ped, ped.save_cat(), ped.save_sex(), ped.save_active()

class TBVGenTable:
    def __init__(self, TBVTable):
        self.TBVtable = pd.read_table(TBVTable, header=None, sep='\s+', names=['Indiv', 'TBV'])
        self.TBVmean = np.mean(self.TBVtable.TBV)
        self.TBVsd = np.std(self.TBVtable.TBV)
        self.TBVvar = np.var(self.TBVtable.TBV)
        self.TBVse = stats.sem(self.TBVtable.TBV)
    
        
        
        
