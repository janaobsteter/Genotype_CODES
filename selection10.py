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
import math
from pylab import legend
from scipy import stats
import matplotlib.pyplot as plt

plt.style.use('ggplot')


class classPed(object):
    def __init__(self, pedfile):
        self.name = pedfile
        self.ped = pd.read_csv(pedfile)
        self.gen = set(self.ped.Generation)
        self.ngen = len(set(self.ped.Generation))
        self.ped['sex'] = ""
        self.ped['cat'] = ""
        self.ped['active'] = ""
        self.ped['PA'] = ""
        self.nrow = len(self.ped)


class pedigree(classPed):
    def __init__(self, pedfile):
        super(pedigree, self).__init__(pedfile)

    def printDate(self):
        print 20042017

    def head(self):
        return self.ped.head()

    def tail(self):
        return self.ped.tail()

    def displayR(self, start, stop):
        return self.ped.ix[start:stop]

    def displayCat(self, cat):
        return self.ped[self.ped.cat == cat]

    def displayInd(self, ind):
        return self.ped[self.ped.Indiv == ind]

    def displayCat_prevGen(self, cat, prevGenDict):
        return self.ped.loc[(self.ped.Indiv.isin(prevGenDict[cat]))]

    def displayInd_prevGen(self, cat, prevGenDict):
        catInd = self.catCurrent_indiv(cat)
        categories = []
        for ind in catInd:
            for (key, value) in prevGenDict.iteritems():
                if ind in value:
                    categories.append(key)
        return categories

    def rows(self):
        return len(self.ped)

    def view_cat(self, cat):
        return self.ped[self.ped.cat == cat]

    def gens(self):
        return set(self.ped['Generation'])

    def compute_age(self):
        self.ped['age'] = max(self.gens()) - self.ped.Generation

    def set_cat(self, start, stop, cat):  # pregled po kategorijah
        self.ped.loc[range(start, (start + stop)), 'cat'] = cat

    def set_cat_list(self, seznam, cat):
        self.ped.loc[seznam, 'cat'] = cat

    def set_cat_gen(self, gen, cat):
        self.ped.loc[self.ped.Generation == gen, 'cat'] = cat

    def set_cat_age(self, age, cat, prevGenDict):
        self.ped.loc[(self.ped['age'] == age) & (self.ped.Indiv.isin(prevGenDict[cat])), 'cat'] = cat

    def set_cat_old(self, oldcat, cat, prevGenDict):  # določi novo kategorijo glede na prejšnjo kategorijo
        self.ped.loc[self.ped.Indiv.isin(prevGenDict[oldcat]), 'cat'] = cat

    def set_cat_age_old(self, age, oldcat, cat, prevGenDict):
        self.ped.loc[(self.ped['age'] == age) & (self.ped.Indiv.isin(prevGenDict[oldcat])), 'cat'] = cat

    def set_cat_sex_old(self, sex, oldcat, cat, prevGenDict):
        self.ped.loc[(self.ped['sex'] == sex) & (self.ped.Indiv.isin(prevGenDict[oldcat])), 'cat'] = cat

    def set_cat_mother_catCurrent(self, mothercat, cat):
        self.ped.loc[(self.ped.Mother.isin(self.catCurrent_indiv(mothercat))), 'cat'] = cat

    def set_sex(self, start, stop, sex):  # pregled po kategorijah
        self.ped.loc[range(start, stop), 'sex'] = sex

    def set_sex_list(self, seznam, sex):
        self.ped.loc[seznam, 'sex'] = sex

    def set_sex_prevGen(self, prevGenDict):
        self.ped.loc[self.ped.Indiv.isin(prevGenDict['F']), 'sex'] = "F"
        self.ped.loc[self.ped.Indiv.isin(prevGenDict['M']), 'sex'] = "M"

    def set_sex_AlphaSim(self, AlphaSimDir):
        # določi spol
        #       gender = pd.read_table(AlphaSimDir + '/Gender_BURNIN.txt', sep='\s+')
        gender = pd.read_table(AlphaSimDir + '/SimulatedData/Gender.txt', sep='\s+')
        females = list(gender[gender.Gender == 2]['Indiv'])
        males = list(gender[gender.Gender == 1]['Indiv'])
        self.set_sex_list(self.ped[self.ped.Indiv.isin(females)].index.tolist(), "F")
        self.set_sex_list(self.ped[self.ped.Indiv.isin(males)].index.tolist(), "M")

    def set_active(self, start, stop, active):  # pregled po kategorijah
        self.ped.loc[range(start, (start + stop)), 'active'] = active

    def set_active_list(self, seznam, active):  # pregled po kategorijah
        self.ped.loc[seznam, 'active'] = active

    def set_active_prevGen(self, prevGenDict):
        self.ped.loc[self.ped.Indiv.isin(prevGenDict['1']), 'active'] = 1
        self.ped.loc[self.ped.Indiv.isin(prevGenDict['2']), 'active'] = 2

    def set_active_age_cat(self, cat, age, active, prevGenDict):
        self.ped.loc[(self.ped['age'] == age) & (self.ped.Indiv.isin(prevGenDict[cat])), 'active'] = active

    def set_active_cat(self, cat, active, prevGenDict):
        self.ped.loc[(self.ped.Indiv.isin(prevGenDict[cat])), 'active'] = active

    def set_active_catCurrent(self, cat, active):
        self.ped.loc[(self.ped['cat'] == cat), 'active'] = active

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
        return len(self.ped[(self.ped.cat == 'nr') & (self.ped.Mother == 0)])

    def father_nr_blank(self):
        return len(self.ped[(self.ped.cat == 'nr') & (self.ped.Father == 0)])

    def active(self):
        return self.ped.active.value_counts()

    def cat_age(self, cat):
        return self.ped.loc[self.ped.cat == cat, 'age'].value_counts()

    def cat_sex(self, cat):
        return self.ped.loc[self.ped.cat == cat, 'sex'].value_counts()

    def catCurrent_indiv(self, cat):
        return list(self.ped.loc[self.ped.cat == cat, 'Indiv'])

    def catCurrent_indiv_age(self, cat, age):
        return list(self.ped.loc[(self.ped.cat == cat) & (self.ped.age == age), 'Indiv'])

    def catCurrent_indiv_sex(self, cat, sex):
        return list(self.ped.loc[(self.ped.cat == cat) & (self.ped.sex == sex), 'Indiv'])

    def catCurrent_indiv_sex_CriteriaRandom(self, cat, sex, number):
        return random.sample(list(self.ped.loc[(self.ped.cat == cat) & (self.ped.sex == sex), 'Indiv']), number)

    def catCurrent_indiv_sex_CriteriaEBV(self, cat, sex, number):
        selRow = list(
            self.ped.loc[(self.ped.cat == cat) & (self.ped.sex == sex), 'EBV'].sort_values(
                ascending=False)[:number].index)  # katere izbereš
        return list(self.ped.Indiv[selRow])

    def catCurrent_indiv_age_CriteriaEBV(self, cat, age, number):
        selRow = list(
            self.ped.loc[(self.ped.cat == cat) & (self.ped.age == age), 'EBV'].sort_values(
                ascending=False)[:number].index)  # katere izbereš
        return list(self.ped.Indiv[selRow])

    def catCurrent_indiv_sex_CriteriaPA(self, cat, sex, number):
        self.computePA_currentCat(cat)
        selRow = list(
            self.ped.loc[(self.ped.cat == cat) & (self.ped.sex == sex), 'PA'].sort_values(
                ascending=False)[:number].index)  # katere izbereš
        return list(self.ped.Indiv[selRow])

    def catCurrent_indiv_sex_criteria(self, cat, sex, criteria, number):
        if criteria == 'random':
            return self.catCurrent_indiv_sex_CriteriaRandom(cat, sex, number)
        if criteria == 'EBV':
            return self.catCurrent_indiv_sex_CriteriaEBV(cat, sex, number)
        if criteria == 'PA':
            return self.catCurrent_indiv_sex_CriteriaPA(cat, sex, number)

    def cat_active(self, cat):
        return self.ped.loc[self.ped.cat == cat, 'active'].value_counts()

    def sex_active(self, sex):
        return self.ped.loc[self.ped.sex == sex, 'active'].value_counts()

    def row_gen(self, gen):
        return list(self.ped.ix[self.ped.Generation == gen].index)

    def row_gen_subset(self, gen, start, stop):
        return list(self.ped.ix[self.ped.Generation == gen][start:stop].index)

    def row_gen_F(self, gen):
        return list(self.ped.ix[self.ped.Generation == gen & self.ped.sex == "F"].index)

    def row_newest_gen(self):
        return list(self.ped.ix[self.ped.Generation == max(self.ped.Generation)].index)

    def row_oldest_gen(self):
        return list(self.ped.ix[self.ped.Generation == min(self.ped.Generation)].index)

    def row_cat(self, cat):
        return list(self.ped.ix[self.ped.cat == cat].index)

    def izloci_poEBV(self, sex, stIzl, oldcat, prevGenDict):
        izlRow = list(
            self.ped.loc[(self.ped.sex == sex) & (self.ped.Indiv.isin(prevGenDict[oldcat])), 'EBV'].sort_values(
                ascending=True)[:stIzl].index)  # katere izločiš
        if len(izlRow) < stIzl:
            print("Too little animals to choose from. <{} {}>".format("izloci po EBV", oldcat))
        else:
            self.set_cat_list(izlRow, "izl")
            self.set_active_list(izlRow, 2)

    def izloci_poEBV_age(self, sex, age, stIzl, oldcat, prevGenDict):
        izlRow = list(self.ped.loc[(self.ped.age == age) & (self.ped.sex == sex) & (
            self.ped.Indiv.isin(prevGenDict[oldcat])), 'EBV'].sort_values(ascending=True)[
                      :stIzl].index)  # katere izločiš
        if len(izlRow) < stIzl:
            print("Too little animals to choose from. <{} {}>".format("izloci po EBV_age", oldcat))
        else:
            self.set_cat_list(izlRow, "izl")
            self.set_active_list(izlRow, 2)

    def izberi_poEBV_top(self, sex, st, oldcat, cat, prevGenDict):
        selRow = list(
            self.ped.loc[(self.ped.sex == sex) & (self.ped.Indiv.isin(prevGenDict[oldcat])), 'EBV'].sort_values(
                ascending=False)[:st].index)  # katere izbereš
        if len(selRow) < st:
            print("Too little animals to choose from. <{} {} > {}>".format("izberi po EBV", oldcat, cat))
        else:
            self.set_cat_list(selRow, cat)
            self.set_active_list(selRow, 1)

    def izberi_poEBV_top_catCurrent(self, sex, st, cat,
                                    newcat):  # to je bolj kot ne samo za prvi krog, ko nimaš slovarja od prejšnje generacije
        selRow = list(self.ped.loc[(self.ped.sex == sex) & (self.ped.cat == cat), 'EBV'].sort_values(ascending=False)[
                      :st].index)  # katere izbereš
        if len(selRow) < st:
            print("Too little animals to choose from. <{} {} > {}>".format("izberi po EBV_catCurrent", cat, newcat))
        else:
            self.set_cat_list(selRow, newcat)
            self.set_active_list(selRow, 1)

    def izberi_poEBV_top_age_naive(self, sex, age, st,
                                   newcat):  # to je bolj kot ne samo za prvi krog, ko nimaš slovarja od prejšnje generacije
        selRow = list(
            self.ped.loc[(self.ped.sex == sex) & (self.ped.age == age) & (self.ped.cat == ""), 'EBV'].sort_values(
                ascending=False)[:st].index)  # katere izbereš
        if len(selRow) < st:
            print("Too little animals to choose from. <{} {} > {}>".format("izberi po EBV_catCurrent", newcat))
        else:
            self.set_cat_list(selRow, newcat)
            self.set_active_list(selRow, 1)

    def izberi_poEBV_top_age(self, sex, age, st, oldcat, cat, prevGenDict):
        selRow = list(self.ped.loc[(self.ped.age == age) & (self.ped.sex == sex) & (
            self.ped.Indiv.isin(prevGenDict[oldcat])), 'EBV'].sort_values(ascending=False)[:st].index)  # katere izbereš
        if len(selRow) < st:
            print("Too little animals to choose from. <{} {} > {}>".format("izberi po EBV_topAge", oldcat, cat))
        else:
            self.set_cat_list(selRow, cat)
            self.set_active_list(selRow, 1)

    def ind_poEBV_top(self, sex, st, cat, prevGenDict):
        selRow = list(self.ped.loc[(self.ped.sex == sex) & (self.ped.Indiv.isin(prevGenDict[cat])), 'EBV'].sort_values(
            ascending=False)[:st].index)  # katere izbereš
        if len(selRow) < st:
            print("Too little animals to choose from. <{} {}>".format("Indiv po EBV_top", cat))
        else:
            self.set_cat_list(selRow, cat)
            self.set_active_list(selRow, 1)

    def izberi_poEBV_OdDo(self, sex, start, stop, oldcat, cat, prevGenDict):
        selRow = list(
            self.ped.loc[(self.ped.sex == sex) & (self.ped.Indiv.isin(prevGenDict[oldcat])), 'EBV'].sort_values(
                ascending=False)[start:stop].index)  # katere izbereš
        if len(selRow) < (stop - start):
            print("Too little animals to choose from. <{} {} > {}>".format("izberi po EBV_OdDo", oldcat, cat))
        else:
            self.set_cat_list(selRow, cat)
            self.set_active_list(selRow, 1)

    def izberi_poEBV_OdDo_age_naive(self, sex, age, start, stop, cat):
        selRow = list(
            self.ped.loc[(self.ped.sex == sex) & (self.ped.age == age) & (self.ped.cat == ""), 'EBV'].sort_values(
                ascending=False)[start:stop].index)  # katere izbereš
        if len(selRow) < (stop - start):
            print("Too little animals to choose from. <{} {} > {}>".format("izberi po EBV_OdDo", cat))
        else:
            self.set_cat_list(selRow, cat)
            self.set_active_list(selRow, 1)

    def izberi_random(self, sex, st, oldcat, cat, prevGenDict):
        freeRow = list(self.ped.loc[(self.ped.sex == sex) & (self.ped.Indiv.isin(prevGenDict[oldcat])) & (
            self.ped.cat == "")].index)  # rows with no assigned category
        if len(freeRow) < st:
            print("Too little animals to choose from. <{} {} > {}>".format("izberi random", oldcat, cat))
        else:
            choice = list(np.random.choice(freeRow, st, replace=False))
            self.set_cat_list(choice, cat)
            self.set_active_list(choice, 1)

    def izberi_random_age(self, sex, age, st, oldcat, cat, prevGenDict):
        freeRow = list(self.ped.loc[(self.ped.age == age) & (self.ped.sex == sex) & (
            self.ped.Indiv.isin(prevGenDict[oldcat])) & (self.ped.cat == "")].index)  # rows with no assigned category
        if len(freeRow) < st:
            print("Too little animals to choose from. <{} {} > {}>".format("izberi random", oldcat, cat))
        else:
            choice = list(np.random.choice(freeRow, st, replace=False))
            self.set_cat_list(choice, cat)
            self.set_active_list(choice, 1)

    def izberi_random_age_naive(self, sex, age, st, cat):
        freeRow = list(self.ped.loc[(self.ped.sex == sex) & (self.ped.age == age) & (
            self.ped.cat == "")].index)  # rows with no assigned category
        if len(freeRow) < st:
            print("Too little animals to choose from. <{} {} > {}>".format("izberi random", cat))
        else:
            choice = list(np.random.choice(freeRow, st, replace=False))
            self.set_cat_list(choice, cat)
            self.set_active_list(choice, 1)

    def izloci_random(self, sex, st, oldcat, prevGenDict):
        freeRow = list(self.ped.loc[(self.ped.sex == sex) & (self.ped.Indiv.isin(prevGenDict[oldcat])) & (
            self.ped.cat == "")].index)  # rows with no assigned category
        if len(freeRow) < st:
            print("Too little animals to choose from. <{} {}>".format("izloči random", oldcat))
        else:
            choice = list(np.random.choice(freeRow, st, replace=False))
            self.set_cat_list(choice, 'izl')
            self.set_active_list(choice, 2)

    def izloci_random_age(self, sex, age, st, oldcat, prevGenDict):
        freeRow = list(self.ped.loc[(self.ped.age == age) & (self.ped.sex == sex) & (
            self.ped.Indiv.isin(prevGenDict[oldcat])) & (self.ped.cat == "")].index)  # rows with no assigned category
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
        self.ped.loc[(self.ped['age'] == age) & (self.ped.Indiv.isin(prevGenDict[cat])), 'cat'] = 'izl'
        self.ped.loc[(self.ped['age'] == age) & (self.ped.Indiv.isin(prevGenDict[cat])), 'active'] = 2

    def save_cat(self):
        cat = defaultdict(list)
        for i in set(self.ped['cat']):
            cat[i] = list(self.ped.loc[self.ped.cat == i, 'Indiv'])
        return cat

    def save_cat(self):
        cat = defaultdict(list)
        for i in set(self.ped['cat']):
            cat[i] = list(self.ped.loc[self.ped.cat == i, 'Indiv'])
        return cat

    def save_active(self):
        active = defaultdict(list)
        for i in set(self.ped['active']):
            active[i] = list(self.ped.loc[self.ped.active == i, 'Indiv'])
        return active

    def save_sex(self):
        sex = defaultdict(list)
        for i in set(self.ped['sex']):
            sex[i] = list(self.ped.loc[self.ped.sex == i, 'Indiv'])
        return sex

    def save_age(self):
        age = defaultdict(list)
        self.ped['age'] = self.ped['Generation'] - max(self.ped['Generation'])
        for i in set(self.ped['age']):
            age[i] = list(self.ped.loc[self.ped.age == i, 'Indiv'])
        return age

    def add_new_gen_naive(self, stNR, potomciNPn):
        nr = pd.DataFrame({'Indiv': range((max(self.ped.Indiv) + 1), (max(self.ped.Indiv) + 1 + stNR)), \
                           'cat': ['nr'] * (stNR - potomciNPn) + ['potomciNP'] * potomciNPn,
                           'sex': ['F', 'M'] * int(stNR / 2), \
                           'Generation': [(max(self.ped.Generation) + 1)] * stNR, \
                           'EBV': list(np.random.uniform(0, 0, stNR)), \
                           'Mother': [0] * stNR,
                           'Father': [0] * stNR,
                           'active': [1] * stNR}, \
                          index=range(max(self.ped.index) + 1, max(self.ped.index) + 1 + stNR))
        self.ped = pd.concat([self.ped, nr])

    def set_mother(self, MotherList):
        first_blank_row = min(self.ped[(self.ped['Mother'] == 0) & ((self.ped.cat == 'nr'))].index)
        self.ped.loc[(first_blank_row):(first_blank_row + len(MotherList) - 1), 'Mother'] = MotherList

    def set_mother_catPotomca(self, MotherList, cat_potomca):
        freeRow = list(self.ped[(self.ped['Mother'] == 0) & (self.ped.cat == cat_potomca)].index)[0:len(MotherList)]
        self.ped.loc[freeRow, 'Mother'] = MotherList

    def set_father(self, FatherList, mothercat, prevGenDict):
        mother_cat_rows = list(
            self.ped[(self.ped.cat == 'nr') & (self.ped.Mother.isin(self.catCurrent_indiv(mothercat)))].index)[
                          0:len(FatherList)]
        self.ped.loc[mother_cat_rows, 'Father'] = FatherList

    def set_father_catPotomca(self, FatherList, cat_potomca):
        freeRow = list(self.ped[(self.ped['Father'] == 0) & (self.ped.cat == cat_potomca)].index)[0:len(FatherList)]
        self.ped.loc[freeRow, 'Father'] = FatherList

    def select_mother_random(self, cat, st):
        pot_Mother = list(self.ped.loc[(self.ped.cat == cat), 'Indiv'])
        return list(random.sample(pot_Mother, st))

    def select_mother_EBV_top(self, cat, st):  # tukaj je trenutna categorija, saj jih doloačs koec generacije
        selRow = list(
            self.ped.loc[(self.ped.cat == cat), 'EBV'].sort_values(ascending=False)[:st].index)  # katere izbereš
        return list(self.ped.loc[selRow, 'Indiv'])

    def write_ped(self, path):
        pedNB = self.ped[self.ped.Generation == max(self.gens())]
        pedNB.to_csv(path, columns=['Indiv', 'Father', 'Mother'], quoting=None, index=False, header=False)

    def write_pedTotal(self, path):
        self.ped.to_csv(path, quoting=None, index=False, header=True)

    def computePA_currentCat(self, cat):  # therefore you have to use this at the end of the selection cycle
        for i in self.ped.loc[self.ped.cat == cat].index:
            self.ped.loc[i, 'PA'] = (float(self.ped.EBV[self.ped.Indiv == self.ped.Father[i]]) + float(
                self.ped.EBV[self.ped.Indiv == self.ped.Mother[i]])) / 2

    def computePA_previousCat(self, oldcat, sex,
                              prevGenDict):  # therefore you have to use this at the end of the selection cycle
        for i in self.ped.loc[(self.ped.sex == sex) & (self.ped.Indiv.isin(prevGenDict[oldcat]))].index:
            self.ped.loc[i, 'PA'] = (float(self.ped.EBV[self.ped.Indiv == self.ped.Father[i]]) + float(
                self.ped.EBV[self.ped.Indiv == self.ped.Mother[i]])) / 2

    def izberi_poEBV_top(self, sex, st, oldcat, cat, prevGenDict):
        selRow = list(
            self.ped.loc[(self.ped.sex == sex) & (self.ped.Indiv.isin(prevGenDict[oldcat])), 'EBV'].sort_values(
                ascending=False)[:st].index)  # katere izbereš
        if len(selRow) < st:
            print("Too little animals to choose from. <{} {} > {}>".format("izberi po EBV", oldcat, cat))
        else:
            self.set_cat_list(selRow, cat)
            self.set_active_list(selRow, 1)

    def UpdateIndCat(self, Dir):
        if not os.path.isfile(Dir + '/IndCat.csv'):
            self.IndCat = pd.DataFrame()
            self.IndCat['Indiv'] = self.ped.Indiv
            self.IndCat['cat0'] = self.ped.cat
        else:
            self.OldIndCat = pd.read_csv(Dir + '/IndCat.csv')
            self.NewIndCat = pd.DataFrame()
            self.NewIndCat['Indiv'] = self.ped.Indiv
            self.NewIndCat['cat' + str(max(self.gens()))] = self.ped.cat
            self.IndCat = pd.merge(self.OldIndCat, self.NewIndCat, on='Indiv', how='outer')
        self.IndCat.to_csv(Dir + '/IndCat.csv', index=None)

    def select_age_0_1(self, categories, nrFn, nrMn, telFn, vhlevljenin, potomciNPn,
                       telMn):  # tukaj odbereš iz novorojenih živali tel, ptel in mlade bike, pripust1
        # FEMALES
        self.set_cat_sex_old("F", "potomciNP", "telF", categories)
        izlF = nrFn - telFn  # koliko jih izločiš
        self.izloci_poEBV("F", izlF, "nr", categories)  # tukaj jih izloči, funkcija v modulu

        self.izberi_poEBV_top("F", (nrFn - izlF), "nr", "telF", categories)  # izberi telice, ki jih osemeniš --> krave

        # MALES
        self.izberi_poEBV_top("M", vhlevljenin, "potomciNP", "vhlevljeni",
                              categories)  # odberi mlade TO SAMO NA ZAČETKU; POTEM POTOMCI BM IN ELITE!
        self.izloci_poEBV("M", int(potomciNPn - vhlevljenin), 'potomciNP', categories)
        self.izberi_random("M", telMn, "nr", "telM", categories)
        self.izloci_random("M", int(nrMn - telMn), "nr", categories)

    def select_age_1_2(self, categories, ptn, mladin, vhlevljenin,
                       bik12n):  # tukaj odbereš nič pri kravah - razen, če so že bikovske matere, pripust 2, bike24
        # FEMALES
        self.izberi_poEBV_top("F", ptn, 'telF', 'pt', categories)
        self.izloci_poEBV("F", (len(categories['telF']) - ptn), 'telF', categories)  # terlice postanejo

        # MALES
        self.izberi_poEBV_top("M", mladin, "vhlevljeni", "mladi", categories)  # odberi mlade
        self.izberi_poEBV_OdDo("M", mladin, vhlevljenin, "vhlevljeni", "pripust1", categories)  # odberi v pripustu
        self.izberi_random("M", bik12n, 'telM', 'bik12', categories)
        self.izloci_random("M", (len(categories['telM']) - bik12n), 'telM', categories)

    def select_age_2_3(self, categories, kraveUp, bmOdbira, bmn, bmUp, cak, pripust1n, pripust2n, mladin):
        # FEMALES
        # najprej dodaj nove krave
        self.set_cat_old('pt', 'k', categories)  # osemenjene telice postanejo krave - predpostavimo, da vse
        # potem izloči najstarejše krave - po 4. laktaciji
        if ('k' in categories.keys()) and ((kraveUp + 2) in self.age()):  # izloči koliko laktacij + 2 leti
            self.izloci_age_cat((kraveUp + 2), 'k', categories)
        # ostale krave prestavi naprej v krave - OZIROMA PODALJŠAJ STATUS!
        self.set_cat_age_old(3, 'k', 'k', categories)
        self.set_cat_age_old(4, 'k', 'k', categories)
        self.set_cat_age_old(5, 'k', 'k', categories)

        # če imaš že dovolj stare krave, potem odberi BM
        # BM se odbira po drugi laktaciji - to je starost 3 - 4 (starost v pedigreju = 3, ker imaš tudi 0)
        if ('k' in categories.keys()) and ((1 + bmOdbira) in self.age()):
            self.izberi_poEBV_top_age("F", 3, int(bmn / bmUp), 'k', 'pBM', categories)  # izberi bikovske matere
        # in izloči najastarejše BM, če jih imaš
        if ('bm' in categories.keys()):
            self.izloci_cat('bm', categories)
        # ostale BM prestavi naprej
        for i in range((1 + bmOdbira + 1), (
                        1 + bmOdbira + bmUp)):  # 1 leto prva osemenitev, bm odbrane po 2. laktaciji, +1 da začneš prestavljat
            self.set_cat_age_old(i, 'pBM', 'pBM', categories)
        self.set_cat_age_old((1 + bmOdbira + bmUp), 'pBM', 'bm',
                             categories)  # spremeni kategorijo iz plemenskih BM v bm v zadnji laktaciji

        # MALES
        # mladi biki postanejo čakajoči (~1 leto, da se osemeni krave s semenom oz. malo po 2. letu)
        self.set_cat_old('mladi', 'cak', categories)
        self.set_active_cat('mladi', 2, categories)

        # čakajočim bikov podaljšaj status (do starosti 5 let)
        # hkrati jim tudi nastavi status izl
        # ped.set_cat_age_old(2, 'cak', 'cak', categories)
        for i in range((2 + 1), (
                    2 + cak)):  # 1 leto, ko začnejo semenit in so mladi biki, 3 so čakajoči, +1 da začneš prestavlajt
            self.set_cat_age_old(i, 'cak', 'cak', categories)

        # povprečna doba v pripustu - glede na to odberi bike, ki preživijo še eno leto
        if 'pripust1' in categories.keys():
            self.izberi_random("M", pripust2n, 'pripust1', 'pripust2', categories)
            self.izloci_random("M", (pripust1n - pripust2n), 'pripust1', categories)

        # plemenske bike prestavljaj naprej
        self.set_cat_old('pb', 'pb', categories)
        self.izloci_cat('bik12', categories)
        self.izloci_cat('pripust2', categories)
        if ('cak' in categories.keys()) and ((cak + 2) in self.age()):  # +2 - eno leto so teleta, eno leto mladi biki
            self.izberi_poEBV_top_age("M", (cak + 2), int(mladin * 0.5), 'cak', 'pb', categories)
            self.set_active_cat('cak', 2,
                                categories)  # tukaj moraš to nastaviti, zato ker fja izberi avtomatsko nastavi na active=1
            self.izloci_poEBV_age("M", (cak + 2), int(mladin * 0.5), 'cak',
                                  categories)  # TUKAJ MORA BITI ŠE STAROST!!!!!!!!!!!

    def doloci_matere(self, stNB, potomciNPn, ptn, kraveUp):
        # MATERE
        sTbmMother = (potomciNPn * 2) if len(self.catCurrent_indiv('pBM')) >= (potomciNPn * 2) else len(
            self.catCurrent_indiv('pBM'))
        print sTbmMother
        print self.cat()
        print self.active()
        if sTbmMother != 0:
            bmMother = self.select_mother_random('pBM', sTbmMother)
            self.set_mother_catPotomca(bmMother, 'potomciNP')
        #

        if 'k' in self.cat():  # TUKAJ SO DOLOČENE SEDAJ VSE MATERE!!!
            mother = self.select_mother_EBV_top('k', int(
                (ptn * kraveUp * 0.9)))  # tukaj odberi brez tistih, ki so za gospodarsko križanje
            if len(mother) >= (
                        stNB - sTbmMother):  # če že imaš dovolj krav, določi matere vsem novorojenim oz. odbiraš matere, saj jih imaš preveč!
                motherOther = random.sample(mother, (stNB - sTbmMother))
                self.set_mother_catPotomca(motherOther, 'nr')  # TUKAJ SO DOLOČENE SEDAJ VSE MATERE!!!
            elif len(mother) < (
                        stNB - sTbmMother):  # če jih še ni dovolj, ne odbiraš mater, ampak uporabiš vse MINUS gosp. križanmja
                self.set_mother_catPotomca(mother, 'nr')

    def izberi_matere(self, stNB, potomciNPn, ptn, kraveUp):
        # MATERE
        sTbmMother = (potomciNPn * 2) if len(self.catCurrent_indiv('pBM')) >= (potomciNPn * 2) else len(
            self.catCurrent_indiv('pBM'))
        print sTbmMother
        print self.cat()
        print self.active()
        if sTbmMother != 0:
            bmMother = self.select_mother_random('pBM', sTbmMother)

        #

        if 'k' in self.cat():  # TUKAJ SO DOLOČENE SEDAJ VSE MATERE!!!
            mother = self.select_mother_EBV_top('k', int(
                (ptn * kraveUp * 0.9)))  # tukaj odberi brez tistih, ki so za gospodarsko križanje
            if len(mother) >= (
                        stNB - sTbmMother):  # če že imaš dovolj krav, določi matere vsem novorojenim oz. odbiraš matere, saj jih imaš preveč!
                motherOther = random.sample(mother, (stNB - sTbmMother))

            elif len(mother) < (
                        stNB - sTbmMother):  # če jih še ni dovolj, ne odbiraš mater, ampak uporabiš vse MINUS gosp. križanmja
                motherOther = mother
            return (list(bmMother), list(motherOther))

    def doloci_ocete(self, stNB, potomciNPn, cak, pbUp, pripustDoz, mladiDoz, pozitivnoTestDoz, NbGenTest,
                     EliteDamsPTBulls, EliteDamsGenBulls, EliteDamsPABulls, gen_mladi, gen_gpb):
        # OČETJE
        mladiOce = self.catCurrent_indiv('mladi')
        pripustOce = self.catCurrent_indiv('pripust1') + self.catCurrent_indiv('pripust2')
        testiraniOce = list(chain.from_iterable([self.catCurrent_indiv_age('pb', (2 + cak + 1 + x)) for x in range(1,
                                                                                                                   pbUp + 1)]))  # v času, ko določaš potomce, so že eno leto starjši!!!
        gentestiraniOce = list(chain.from_iterable([self.catCurrent_indiv_age('gpb', x + 1) for x in range(1,
                                                                                                           pbUp + 1)]))  # v času, ko določaš potomce, so že eno leto starjši!!!
        mladiOceBest = self.catCurrent_indiv_sex_CriteriaEBV('mladi', 'M', 4)
        cakOcetjeBest = list(
            chain.from_iterable([self.catCurrent_indiv_age_CriteriaEBV('cak', (2 + x), 4) for x in range(1, cak + 1)]))
        print 'GenTest{0}'.format(str(len(gentestiraniOce)))

        # set fathers for offspring of contracted mating
        bmMother = (potomciNPn * 2) if len(self.catCurrent_indiv('pBM')) >= (potomciNPn * 2) else len(
            self.catCurrent_indiv('pBM'))  # the number of elite dams - they are the limiting factor
        # for classical testing - if you already have elite bulls
        if EliteDamsPTBulls:  # whether the elite dams are inseminated with genomicaly tested or progeny tested bulls
            if testiraniOce:
                elita = np.random.choice(testiraniOce, bmMother, replace=True)
                #       pd.Series(elita).value_counts()#preveri zastopanost po bikih
                # naštimaj očete elite --> BM
        if EliteDamsGenBulls:  # if with genomically tested
            if gen_mladi:
                genMladiOcetje = mladiOceBest + cakOcetjeBest
                elita = np.random.choice(genMladiOcetje, bmMother, replace=True)
            if gen_gpb:
                if gentestiraniOce:  # if you already have genomically proven bulls
                    elita = np.random.choice(gentestiraniOce, bmMother, replace=True)
                else:  # otherwise inseminate with progeny tested until genomically tested become proven
                    elita = np.random.choice(testiraniOce, bmMother, replace=True)

        if EliteDamsPABulls:
            self.computePA_previousCat('potomciNP', 'M', categories)
            selRow = list(
                self.ped.loc[(self.ped.cat == cat) & (self.ped.sex == sex), 'PA'].sort_values(
                    ascending=False)[:number].index)  # katere izbereš
            return list(self.ped.Indiv[selRow])

        self.set_father_catPotomca(elita, 'potomciNP')

        # this if for the rest of the new born population - 'mix' semen of all potential fathers
        # first - according to the given % - how many offsprign will be produced by genomically tested bulls

	 # here choose the classically tested bulls --> get the doses
        ClassOcetje = list(
            testiraniOce * pozitivnoTestDoz + mladiOce * mladiDoz) if testiraniOce else []  # progeny teste fathers - keep the ratios between young / natural service / PT
        PripustOcetje = list(pripustOce * pripustDoz)
        if NbGenTest == stNB:  # če naj bo kompletna splošna populacija semenjena samo z genomskimi (mlade, čakajoče si tako dala v gpb)
            GenOcetje = list(
                gentestiraniOce * pozitivnoTestDoz) if 'gpb' in self.cat() else []  # dokler še imaš progene, uporabljaj mešano seme, potem ostanejo tako samo genomsko
            if len(GenOcetje) + len(PripustOcetje) < NbGenTest:
                classNB = NbGenTest - (len(GenOcetje) + len(PripustOcetje))
                Ocetje = random.sample(ClassOcetje, classNB - potomciNPn * 2)  + GenOcetje + PripustOcetje
            if len(GenOcetje) + len(PripustOcetje) >= NbGenTest:
                Ocetje = random.sample(GenOcetje + PripustOcetje, stNB - potomciNPn * 2)
        if NbGenTest == 0:  # če je % semenjenih z genomskimi 0, semeni samo z pb, mladimi in pripustom
            Ocetje = random.sample(ClassOcetje + PripustOcetje, stNB - potomciNPn * 2)
        if NbGenTest > 0 and NbGenTest < stNB:  # če je odstotek nekje med 0 in 100, semeni točno določen del z genomskimi, preostalo mix klasike in pripusta
            GenOcetje = list(
                np.random.choice(gentestiraniOce, (NbGenTest), replace=True)) if 'gpb' in self.cat() else []
            ClassOcetje = random.sample(ClassOcetje + PripustOcetje, stNB - NbGenTest - potomciNPn * 2)
            Ocetje = GenOcetje + ClassOcetje

        self.set_father_catPotomca(Ocetje, 'nr')

    def izberi_ocete_gen(self, stNB, potomciNPn, cak, pbUp, pripustDoz, mladiDoz, pozitivnoTestDoz, NbGenTest,
                     EliteDamsPTBulls, EliteDamsGenBulls, EliteDamsPABulls, gen_mladi, gen_gpb):
        # OČETJE
        gentestiraniOce = list(chain.from_iterable([self.catCurrent_indiv_age('gpb', x + 1) for x in range(1, pbUp + 1)]))  # v času, ko določaš potomce, so že eno leto starjši!!!
        return (gentestiraniOce)

    def save_cat_DF(self):
        categoriesDF = pd.DataFrame.from_dict(self.save_cat(), orient='index').transpose()
        categoriesDF.to_csv('Categories_gen' + str(max(self.gens())) + 'DF.csv', index=None)

    def save_sex_DF(self):
        sexDF = pd.DataFrame.from_dict(self.save_sex(), orient='index').transpose()
        sexDF.to_csv('Sex_gen' + str(max(self.gens())) + 'DF.csv', index=None)

    def save_active_DF(self):
        activeDF = pd.DataFrame.from_dict(self.save_active(), orient='index').transpose()
        activeDF.to_csv('Active_gen' + str(max(self.gens())) + 'DF.csv', index=None)

    def create_categoriesDict(self, catDFEx):
        categories = defaultdict(list)
        catDF = pd.read_csv(catDFEx)
        for cat in catDF.columns:
            values = [int(i) for i in catDF[cat] if not math.isnan(i)]
            categories[cat] = values
        return categories

    def create_sexDict(self, sexDFEx):
        sexDict = defaultdict(list)
        sexDF = pd.read_csv(sexDFEx)
        for sex in sexDF.columns:
            values = [int(i) for i in sexDF[sex] if not math.isnan(i)]
            sexDict[sex] = values
        return sexDict

    def create_activeDict(self, activeDFEx):
        activeDict = defaultdict(list)
        activeDF = pd.read_csv(activeDFEx)
        for active in activeDF.columns:
            values = [int(i) for i in activeDF[active] if not math.isnan(i)]
            activeDict[active] = values
        return activeDict

    def saveIndForGeno(self, genotypedCat):
        if os.path.isfile('IndForGeno.txt'):
            pd.DataFrame({0: sorted(list(set
                                         (chain.from_iterable([self.catCurrent_indiv_sex_criteria(x, sex, xC,
                                                                                                  int((len(
                                                                                                      self.catCurrent_indiv_sex(
                                                                                                          x, sex)) * (
                                                                                                       xP / 100))))
                                                               if xP != 100 else self.catCurrent_indiv_sex(x, sex) for
                                                               (x, xP, xC, sex) in genotypedCat]))))}).to_csv(
                # if self.sel == 'class': #ZDJ TEGA NI, KER JE POSEBEJ FILE!
                'IndForGeno_new.txt', index=None, header=None)
            os.system("grep -v -f IndForGeno.txt IndForGeno_new.txt > uniqNew && mv uniqNew IndForGeno_new.txt")
            os.system(
                'cat IndForGeno_new.txt IndForGeno.txt | sort -n| uniq > IndGenTmp && mv IndGenTmp IndForGeno.txt')
        else:
            pd.DataFrame({0: sorted(list(set
                                         (chain.from_iterable([self.catCurrent_indiv_sex_criteria(x, sex, xC,
                                                                                                  int((len(
                                                                                                      self.catCurrent_indiv_sex(
                                                                                                          x, sex)) * (
                                                                                                       xP / 100.0))))
                                                               for (x, xP, xC, sex) in genotypedCat]))))}).to_csv(
                'IndForGeno.txt', index=None, header=None)

    def updateAndSaveIndForGeno(self, genotypedCat, rmNbGen, removesex, AlphaSimDir):
        if os.path.isfile('IndForGeno.txt'):
            # first obtain new animals for genotypisation
            pd.DataFrame({0: sorted(list(set
                                         (chain.from_iterable([self.catCurrent_indiv_sex_criteria(x, sex, xC,
                                                                                                  int((len(
                                                                                                      self.catCurrent_indiv_sex(
                                                                                                          x, sex)) * (
                                                                                                       xP / 100))))
                                                               if xP != 100 else self.catCurrent_indiv_sex(x, sex) for
                                                               (x, xP, xC, sex) in genotypedCat]))))}).to_csv(
                # if self.sel == 'class': #ZDJ TEGA NI, KER JE POSEBEJ FILE!
                'IndForGeno_new.txt', index=None, header=None)
            # then remove old animals from the previous file
            inds = pd.read_table('IndForGeno.txt', header=None, names=['Indiv'])
            ped = pd.read_table(AlphaSimDir + '/SimulatedData/PedigreeAndGeneticValues_cat.txt', sep='\s+')
            inds = pd.merge(inds, ped[['Indiv', 'Generation']],
                            on='Indiv')  # obtain generation of the previously genotyped Individuals
            inds = pd.merge(inds, ped[['Indiv', 'sex']],
                            on='Indiv')  # obtain sex of the previously genotyped Individuals
            sexDF = inds.loc[(inds.sex == removesex)]
            pd.concat([inds.loc[inds.sex != removesex], sexDF.loc[
                sexDF.Generation.isin(sorted(list(set(sexDF.Generation)))[rmNbGen:])]]).sort_values(by='Indiv')[
                'Indiv'].to_csv(
                'IndForGeno.txt', index=None, header=None)
            os.system(
                "grep -v -f IndForGeno.txt IndForGeno_new.txt > uniqNew && mv uniqNew IndForGeno_new.txt")  # obtain only new ones - do you dont get duplicate cows - ni nepotrebno - to zato, da ti potegne le nove genotipe za GenoFiel
            os.system(
                'cat IndForGeno_new.txt IndForGeno.txt | sort -n| uniq > IndGenTmp && mv IndGenTmp IndForGeno.txt')
        else:
            pd.DataFrame({0: sorted(list(set
                                         (chain.from_iterable([self.catCurrent_indiv_sex_criteria(x, sex, xC,
                                                                                                  int((len(
                                                                                                      self.catCurrent_indiv_sex(
                                                                                                          x, sex)) * (
                                                                                                       xP / 100.0))))
                                                               for (x, xP, xC, sex) in genotypedCat]))))}).to_csv(
                'IndForGeno.txt', index=None, header=None)


class OrigPed(object):
    def __init__(self, AlphaSimDir, codeDir):
        self.name = AlphaSimDir + '/SimulatedData/PedigreeAndGeneticValues.txt'
        self.pdPed = pd.read_table(self.name, sep='\s+')
        self.AlphaSimDir = AlphaSimDir
        self.codeDir = codeDir

    def addInfo(self):
        pedTotal = pd.read_csv(self.AlphaSimDir + 'ExternalPedigreeTotal.txt')
        self.pdPed.loc[:, 'cat'] = pedTotal.cat
        self.pdPed.loc[:, 'age'] = pedTotal.age
        self.pdPed.loc[:, 'sex'] = pedTotal.sex
        self.pdPed.loc[:, 'active'] = pedTotal.active
        self.pdPed.to_csv(self.AlphaSimDir + '/SimulatedData/PedigreeAndGeneticValues_cat.txt', index=None, sep=" ")

    def computeEBV(self, cor):
        # precacunas EBV v Ru in zapises PEDIGRE
        shutil.copy(self.codeDir + "/Rcorr_PedEBV.R", "Rcorr_PedEBV_ThisGen.R")
        os.system('sed -i "s|AlphaSimPed|' + self.name + '|g" Rcorr_PedEBV_ThisGen.R')
        os.system('sed -i "s|setCor|' + str(cor) + '|g" Rcorr_PedEBV_ThisGen.R')
        call('Rscript Rcorr_PedEBV_ThisGen.R', shell=True)


class blupf90:
    def __init__(self, AlphaSimDir, codeDir, way=None):
        self.blupgenParamFile = codeDir + '/renumf90.par'
        self.blupgenParamFile_Clas = codeDir + '/renumf90_Clas.par'
        # self.blupgenParamFile = '/home/jana/Genotipi/Genotipi_CODES/blupf90_Selection'
        # self.blupParamFile = AlphaSimDir + 'blupf90_Selection'
        if way == 'milk':
            self.AlphaPed = pd.read_table(AlphaSimDir + '/SimulatedData/PedigreeAndGeneticValues_cat.txt', sep=' ')
            # self.AlphaGender = pd.read_table(AlphaSimDir + '/SimulatedData/Gender.txt', sep='\s+')
            self.AlphaSimDir = AlphaSimDir
            self.gen = max(self.AlphaPed['Generation'])
            self.animals = len(self.AlphaPed)
            self.blupPed = self.AlphaPed.loc[:, ['Indiv', 'Father', 'Mother']]
            self.blupDatT = self.AlphaPed.loc[:, ['Indiv', 'phenoNormUnres1', 'cat', 'sex', 'age', 'active']]
        if way == 'burnin_milk':
            self.AlphaSimDir = AlphaSimDir
            self.AlphaPed = pd.read_table(AlphaSimDir + '/SimulatedData/PedigreeAndGeneticValues.txt', sep='\s+')
            self.AlphaGender = pd.read_table(AlphaSimDir + '/SimulatedData/Gender.txt', sep='\s+')
            self.AlphaPed.loc[:, 'sex'] = self.AlphaGender.Gender
            self.gen = max(self.AlphaPed['Generation'])
            self.animals = len(self.AlphaPed)
            self.blupDatT = self.AlphaPed.loc[:, ['Indiv', 'phenoNormUnres1', 'sex']]
            self.blupPed = self.AlphaPed.loc[:, ['Indiv', 'Father', 'Mother']]

    def makePed_class(self):
        self.blupPed.loc[((self.blupPed['Mother'] != 0) & (self.blupPed['Father'] != 0)), 'Code'] = 1
        self.blupPed.loc[((self.blupPed['Mother'] == 0) & (self.blupPed['Father'] != 0)), 'Code'] = 2
        self.blupPed.loc[((self.blupPed['Mother'] != 0) & (self.blupPed['Father'] == 0)), 'Code'] = 2
        self.blupPed.loc[((self.blupPed['Mother'] == 0) & (self.blupPed['Father'] == 0)), 'Code'] = 3
        # df['Code'] = df.Code.astype(int)
        self.blupPed.to_csv(self.AlphaSimDir + 'Blupf90.ped', float_format="%.0f", header=None, index=False, sep=" ")

    def makePed_gen(self):
        self.blupPed.to_csv(self.AlphaSimDir + 'Blupf90.ped', float_format="%.0f", header=None, index=False, sep=" ")

    def deletePhenotype_cat(self, listUnphenotyped):
        for i in listUnphenotyped:
            self.blupDatT.loc[self.blupDatT.cat == i, 'phenoNormUnres1'] = 0

    def deletePhenotype_sexage(self, sex, age):
        self.blupDatT.loc[(self.blupDatT.sex == sex) & (self.blupDatT.age == age), 'phenoNormUnres1'] = 0

    def deletePhenotype_sex(self, sex):
        self.blupDatT.loc[self.blupDatT.sex == sex, 'phenoNormUnres1'] = 0

    def deleteInd_sex(self, sex):
        self.blupDatT = self.blupDatT.loc[self.blupDatT.sex != sex]

    def makeDat_sex(self, sex):  # THIS IS not for the repeatability model
        # first remove animals that do not have phenotype (chosen sex)
        self.blupDatT = self.blupDatT.loc[self.blupDatT.sex == sex]
        self.blupDatT[['Indiv', 'phenoNormUnres1', 'sex']].to_csv(self.AlphaSimDir + 'Blupf90.dat', header=None,
                                                                  index=False, sep=" ")

    # this is for the repeatability model since it merges previous dat file with the new one
    # this one deletes phenotypes of the given sex
    def makeDat_removePhen_milk(self):
        # first remove phenotype from animals that do not have phenotype
        # self.deletePhenotype_sex('M') # odstrani fenotip moškim živalim
        # self.deletePhenotype_cat(['potomciNP', 'nr', 'telF', 'pt']) #odstrani fenotip ženskim telicam
        # merge with previous dat file - if there is one - add only the phenotypes of ACTIVE individuals!!!
        if os.path.isfile(self.AlphaSimDir + 'Blupf90.dat'):
            blupDatOld = pd.read_csv('Blupf90.dat', sep=" ", names=['Indiv', 'phenoNormUnres1',
                                                                    'sex'])  # to je dat iz prejšnjega kroga selekcija
            # print 'datOld{0}'.format(str(len(blupDatOld)))
            # pusti samo živali, ki imajo fenotipe (za mleko so to krave in bikovske matere)
            blupDatNew = self.blupDatT.loc[
                (self.blupDatT.active == 1) & ((self.blupDatT.cat == 'k') | (self.blupDatT.cat == 'bm')
                                               | (self.blupDatT.cat == 'pBM')), ['Indiv', 'phenoNormUnres1', 'sex']]
            # blupDatNew = self.blupDatT.loc[(self.blupDatT.active == 1), ['Indiv', 'phenoNormUnres1', 'sex']] #tukaj izberi samo krave in bikovske matere - aktivne!
            # print 'datNew{0}'.format(str(len(blupDatNew)))
            # blupDatNew.loc[blupDatNew.sex == 'M', 'sex'] = 1
            blupDatNew.loc[blupDatNew.sex == 'F', 'sex'] = 2
            pd.concat([blupDatOld, blupDatNew]).to_csv(self.AlphaSimDir + 'Blupf90.dat', header=None, index=False,
                                                       sep=" ")  # dodaj fenotip
            # os.system('sed -i "s/ 0.0 / 0 /g" Blupf90.dat')
            # else: #if there is none create the file with all individuals! --> This is after burn-in
            # self.blupDatT = self.blupDatT.loc[self.blupDatT.sex == 'F']
            # self.blupDatT.loc[self.blupDatT.sex == 'F', 'sex'] = 2
            # self.blupDatT[['Indiv', 'phenoNormUnres1', 'sex']].to_csv(self.AlphaSimDir + 'Blupf90.dat', header=None, index=False, sep=" ")    #this is for the repeatability model since it merges previous dat file with the new one
            # os.system('sed -i "s/ 0.0 / 0 /g" Blupf90.dat')

    # this one deletes phenotypes according to given categories
    def makeDat_removePhen_cat(self, listUnphenotyped):
        # first remove phenotype from animals that do not have phenotype
        self.deletePhenotype_cat(listUnphenotyped)
        # merge with previous dat file - if there is one
        if os.path.isfile(self.AlphaSimDir + 'Blupf90.dat'):
            blupDatOld = pd.read_csv('Blupf90.dat', sep=" ", names=['Indiv', 'phenoNormUnres1', 'sex'])
            blupDatNew = blupDatT[['Indiv', 'phenoNormUnres1', 'sex']]
            pd.concat([blupDatOld, blupDatNew]).to_csv(self.AlphaSimDir + 'Blupf90.dat', header=None, index=False,
                                                       sep=" ")
        else:
            self.blupDatT[['Indiv', 'phenoNormUnres1', 'sex']].to_csv(self.AlphaSimDir + 'Blupf90.dat', header=None,
                                                                      index=False, sep=" ")

    def setNumberAnimals(self, blupParamFile):
        os.system('sed -i "s|NumberOfAnimals|' + str(self.animals) + '|g" ' + blupParamFile)

    def setResidualVariance(self, resvar, blupParamFile):
        os.system('sed -i "s|ResidualVariance|' + str(resvar) + '|g" ' + blupParamFile)

    def setGeneticVariance(self, genvar, blupParamFile):
        os.system('sed -i "s|GeneticVariance|' + str(genvar) + '|g" ' + blupParamFile)

    def prepareSelPed(self):
        blupSol = pd.read_csv(self.AlphaSimDir + '/renumbered_Solutions_' + str(self.gen), header=None,
                              sep='\s+', names=['renID', 'ID', 'Solution'])
        AlphaSelPed = self.AlphaPed.loc[:, ['Generation', 'Indiv', 'Father', 'Mother', 'gvNormUnres1']]
        # blupSolRandom = blupSol.loc[blupSol.Effect == 1] Če imaš še fixed effect
        AlphaSelPed.loc[:, 'EBV'] = blupSol.Solution
        AlphaSelPed.to_csv(self.AlphaSimDir + 'GenPed_EBV.txt', index=None)

    def preparePedDat_cat(self, listUnphenotyped):
        self.makePed()
        assert isinstance(listUnphenotyped, object)
        self.deletePhenotype_cat(listUnphenotyped)
        self.makeDat()
        shutil.copy(self.blupgenParamFile, self.AlphaSimDir)

    def preparePedDat_sex(self):
        self.makePed()
        self.deletePhenotype_sex('M')
        self.makeDat()
        shutil.copy(self.blupgenParamFile, self.AlphaSimDir)

    def prepareParamFiles(self, genvar, resvar, blupParamFile):
        self.setNumberAnimals(blupParamFile)
        self.setGeneticVariance(genvar, blupParamFile)
        self.setResidualVariance(resvar, blupParamFile)


class accuracies:
    def __init__(self, AlphaSimDir):
        self.AlphaSimDir = AlphaSimDir
        self.accuracies = defaultdict()
        if os.path.isfile(self.AlphaSimDir + 'Accuracies_Gen.csv'):
            self.accuraciesGen = pd.read_table(self.AlphaSimDir + 'Accuracies_Gen.csv', sep=",")
        else:
            self.accuraciesGen = pd.DataFrame()
        if os.path.isfile(self.AlphaSimDir + 'Accuracies_Cat.csv'):
            self.accuraciesCat = pd.read_table(self.AlphaSimDir + 'Accuracies_Cat.csv', sep=",")
        else:
            self.accuraciesCat = pd.DataFrame()

    def saveAcc(self):
        name = self.AlphaSimDir + '/SimulatedData/PedigreeAndGeneticValues_cat.txt'
        pdPed = pd.read_table(name, sep=' ')
        gen = max(pdPed['Generation'])
        EBV = pd.read_table(self.AlphaSimDir + 'renumbered_Solutions_' + str(gen), header=None,
                            sep='\s+', names=['renID', 'ID', 'solution'])
        pdPed.loc[:, 'EBV'] = list(EBV.solution)
        cor = stats.pearsonr(pdPed.EBV, pdPed.gvNormUnres1)[0]
        self.accuracies[gen] = cor
        corGenX = pd.DataFrame({gen: list(pdPed.groupby('Generation')[['gvNormUnres1', 'EBV']].corr().ix[0::2, 'EBV'])})
        self.accuraciesGen = pd.concat([self.accuraciesGen, corGenX], ignore_index=True, axis=1)
        corCatX = pd.DataFrame(
            {'cat': list((pdPed.groupby('cat')[['gvNormUnres1', 'EBV']].corr().ix[0::2, 'EBV'].index.levels)[0]),
             gen: list(pdPed.groupby('cat')[['gvNormUnres1', 'EBV']].corr().ix[0::2, 'EBV'])})
        try:
            self.accuraciesCat = pd.merge(self.accuraciesCat, corCatX, on='cat', how='outer')
        except:
            self.accuraciesCat = pd.concat([self.accuraciesCat, corCatX], ignore_index=True)

    def writeAcc(self):
        pd.DataFrame({'Gen': self.accuracies.keys(), 'r': self.accuracies.values()}) \
            .to_csv(self.AlphaSimDir + 'AccuraciesBV.csv', index=None)
        self.accuraciesGen.to_csv(self.AlphaSimDir + 'Accuracies_Gen.csv', sep=",", index=None)
        self.accuraciesCat.to_csv(self.AlphaSimDir + 'Accuracies_Cat.csv', sep=",", index=None)

    def removeFiles(self):
        if os.path.isfile(self.AlphaSimDir + 'Accuracies_Cat.csv'):
            os.remove(self.AlphaSimDir + 'Accuracies_Cat.csv')
        if os.path.isfile(self.AlphaSimDir + 'Accuracies_Gen.csv'):
            os.remove(self.AlphaSimDir + 'Accuracies_Gen.csv')
            # if os.path.isfile(self.AlphaSimDir + 'AccuraciesBV.csv'):
            # os.remove(self.AlphaSimDir + 'AccuraciesBV.csv')


class AlphaSim_OutputFile:
    def __init__(self, AlphaSimDir):
        self.AlphaSimDir = AlphaSimDir
        self.TraitVariance = pd.read_csv(self.AlphaSimDir + '/SimulatedData/TraitVarianceComponents.txt', skiprows=1,
                                         sep='\s+', skipfooter=2, engine='python')

    def getAddVar(self):
        return self.TraitVariance.NormalModelUnres['VarA']

    def getResVar(self):
        return self.TraitVariance.NormalModelUnres['VarE']

    def getGenVar(self):
        return self.TraitVariance.NormalModelUnres['VarG']

    def getDomVar(self):
        return self.TraitVariance.NormalModelUnres['VarD']


################################################################
# FUNKCIJE
###################################################################
def selekcija_total(pedFile, **kwargs):
    print kwargs
    ped = pedigree(pedFile)

    # določi spol
    # ped.set_sex_AlphaSim(kwargs.get('AlphaSimDir'))

    # tukaj potem pridobi kategorije - če imaš samo eno burn-in, štartaš iz nule
    if max(ped.gen) == 1:  # ČE SAMO ENA GENERACIJA V PEDIGREJU! - to zdj ne delaš več tako
        ped.set_cat_gen(max(ped.gen), "nr")  # to je samo na prvem loopu
        ped.set_sex_list([x for x in range(0, ped.rows()) if x % 2 == 0], "F")
        ped.set_sex_list([x for x in range(0, ped.rows()) if x % 2 != 0], "M")
        # določi spol
        # ped.set_sex_AlphaSim(kwargs.get('AlphaSimDir'))
        ped.izberi_poEBV_top_catCurrent("F", int(kwargs.get('potomciNPn')), 'nr', 'potomciNP')
        ped.izberi_poEBV_top_catCurrent("M", int(kwargs.get('potomciNPn')), 'nr', 'potomciNP')

        # global categories #to moraš dat global samo v prvenm loopu, drugje dobiš return
        categories = ped.save_cat()
        # global sex
        sex = ped.save_sex()
        active = ped.save_active()
        ped = pedigree(pedFile)

    elif max(ped.gens()) > 1:
        categories = ped.create_categoriesDict('Categories_gen' + str(max(ped.gens())) + 'DF.csv')
        sex = ped.create_sexDict('Sex_gen' + str(max(ped.gens())) + 'DF.csv')
        active = ped.create_activeDict('Active_gen' + str(max(ped.gens())) + 'DF.csv')

        # določi spol po AlphaSIm - ampak to ni potrebno, ker je vseskozi external pedigre
    # ped.set_sex_AlphaSim(kwargs.get('AlphaSimDir'))

    ped.set_sex_prevGen(sex)  # add sex information for individuals from prevGen
    ped.set_active_prevGen(active)  # add active information for individuals from prevGen

    # remove category information from the ped itself
    for i in ped.gens():
        ped.set_cat_gen(i, "")

    # transfer culled (izlocene) category from prevGen
    ped.set_cat_old('izl', 'izl', categories)

    # compute age of the animals in the current selection year
    ped.compute_age()

    #################################################
    # FEMALES
    #################################################
    # age 0 - here you have newborn females (NB & potomkeNP) --> nekaj jih izloči, druge gredo naprej do ženskih telet
    ped.set_cat_sex_old("F", "potomciNP", "telF", categories)  # potomke načrtnih parjenj gredo v telice
    ped.set_active_catCurrent('telF', 1)
    izlF = int(kwargs.get('nrFn')) - int(kwargs.get('telFn'))  # number of culles NB females
    ped.izberi_poEBV_top("F", kwargs.get('telFn'), "nr", "telF",
                         categories)  # izberi NB ženske, ki preživijo in postanejo telice
    ped.izloci_poEBV("F", izlF, "nr", categories)  # cull females (lowest on EBV) tukaj jih izloči, funkcija v modulu

    # age 1 - pri enem letu osemeni določeno število telic (% določen zgoraj), druge izloči
    if 'telF' in categories.keys():
        ped.izberi_poEBV_top("F", kwargs.get('ptn'), 'telF', 'pt', categories)  # plemenske telice
        ped.izloci_poEBV("F", (len(categories['telF']) - kwargs.get('ptn')), 'telF', categories)  # preostale izloči

    # age > 2 - tukaj odbiraš in izločaš krave, odbiraš in izločaš BM
    # najprej dodaj nove krave, če jih že imaš v populaciji
    if ('pt' in categories.keys()):  # če imaš v pedigreju plemenske telice
        ped.set_cat_old('pt', 'k', categories)  # osemenjene telice postanejo krave - predpostavimo, da vse
        ped.set_active_cat('pt', 1, categories)
        ped.set_active_cat('k', 1, categories)

    # če imaš že dovolj stare krave, potem odberi BM - NAJPREJ - kasneje funkcija random polni prazne prostore (cat)!!!
    # BM se odbira po drugi laktaciji - to je starost 3 - 4 (starost v pedigreju = 3, ker imaš tudi 0)
    if ('k' in categories.keys()) and ((1 + kwargs.get('bmOdbira')) in ped.age()):
        ped.izberi_poEBV_top_age("F", kwargs.get('bmOdbira') + 1, int(kwargs.get('bmn') / kwargs.get('bmUp')), 'k',
                                 'pBM',
                                 categories)  # izberi BM, ki jih osemeniš (plemenske BM = pBM) iz krav po 2. laktaciji
        ped.set_active_cat('pBM', 1, categories)
    # ostale BM prestavi naprej - BM po 1. do izločitvene laktacije
    if 'pBM' in categories.keys():
        for i in range((1 + kwargs.get('bmOdbira') + 1), (
                        1 + kwargs.get('bmOdbira') + kwargs.get(
                    'bmUp'))):  # 1 leto prva osemenitev, bm odbrane po 2. laktaciji, +1 da začneš prestavljat
            ped.set_cat_age_old(i, 'pBM', 'pBM', categories)
        # spremeni kategorijo iz plemenskih BM v bm v zadnji laktaciji
        ped.set_cat_age_old((1 + kwargs.get('bmOdbira') + kwargs.get('bmUp')), 'pBM', 'bm',
                            categories)



        # potem odberi krave - na random, vsako leto za določeno število manj krav
    for i in range(2 + 1, (
                2 + kwargs.get('kraveUp'))):  # 2 + 1 - pri dveh letih prva laktacija, prestavljati začneš leto po tem
        ped.izberi_random_age('F', i, (
            kwargs.get('ptn') - (kwargs.get('MinusDamLact') * (i - 2)) - int(kwargs.get('bmn') / kwargs.get('bmUp'))),
                              'k',
                              'k', categories)
        ped.izloci_random_age('F', i, (kwargs.get('MinusDamLact')), 'k', categories)

    # potem izloči najstarejše krave - po 4. laktaciji
    if ('k' in categories.keys()) and ((kwargs.get('kraveUp') + 2) in ped.age()):  # izloči koliko laktacij + 2 leti
        ped.izloci_age_cat((kwargs.get('kraveUp') + 2), 'k', categories)

    # in izloči najastarejše BM, če jih imaš
    if ('bm' in categories.keys()):
        ped.izloci_cat('bm', categories)

    #################################################################
    # MALES
    #################################################################
    # age 0: štartaš z NB in potomci NP --> odbereš vhlevljene iz potomcev NP in moška teleta in NB
    # NAJPREJ DELI, KI SO SKUPNI PROGENEMU IN GENOMSKEMU TESTIRANJU
    ped.izberi_random("M", kwargs.get('telMn'), "nr", "telM", categories)  # izberi moška teleta, ki preživijo (random)
    ped.izloci_random("M", int(kwargs.get('nrMn') - kwargs.get('potomciNPn') - kwargs.get('telMn')), "nr",
                      categories)  # druga teleta izloči

    if 'telM' in categories.keys():
        ped.izberi_random("M", kwargs.get('bik12n'), 'telM', 'bik12',
                          categories)  # random odberi bike, ki preživijo do 2. leta
        ped.izloci_random("M", (len(categories['telM']) - kwargs.get('bik12n')), 'telM', categories)  # izloči preostale

    if 'bik12' in categories.keys():  # izloci bike nad 2. leti
        ped.izloci_cat('bik12', categories)

    # Tukaj deli, kjer se progeni in genomski testi razlikujeta
    # progeni: potomciNP -> vhlevljeni -> mladi -> čakajoči -> pozitivno testirani
    # če je PROGENI TESTIRANJE
    if kwargs.get('EBV'):
        ped.izberi_random("M", kwargs.get('vhlevljenin'), "potomciNP", "vhlevljeni",
                          categories)  # vhlevi najboljše potomceNP - RANDOM
        ped.izloci_random('M', int(kwargs.get('potomciNPn') - kwargs.get('vhlevljenin')), 'potomciNP',
                          categories)  # druge potomceNP izloči


        # age1: tukaj odbereš mlade iz vhlevljenih bikov in bike, ki preživijo do drugega leta
    if 'vhlevljeni' in categories.keys():  # če imaš vhlevljene bike (samo v progenem testu)
        #############################
        # to je za potrebe prehoda
        if kwargs.get("genTest_gpb") and 'vhlevljeni' in [i[0] for i in kwargs[
            'genotyped']]:  # če se greš genomsko shemo, jih takoj pogenotipiziraj in prestavi v gpb (5)
            ped.izberi_poEBV_top("M", kwargs.get('genpbn'), "vhlevljeni", "gpb", categories)  # odberi gpb
            ped.izberi_poEBV_OdDo("M", kwargs.get('genpbn'), (kwargs.get('genpbn') + kwargs.get('pripust1n')),
                                  'vhlevljeni', 'pripust1',
                                  categories)  # preostali genomsko testirani gredo v pripust
            ped.izloci_poEBV("M", (kwargs.get("vhlevljenin") - (kwargs.get('genpbn') + kwargs.get('pripust1n'))),
                             "vhlevljeni", categories)
        #############################
        else:  # drugače 8 v mlade
            ped.izberi_poEBV_top("M", kwargs.get('mladin'), "vhlevljeni", "mladi", categories)  # odberi mlade
            ped.izberi_poEBV_OdDo("M", kwargs.get('mladin'), kwargs.get('vhlevljenin'), "vhlevljeni", "pripust1",
                                  categories)  # preostali vhlevljeni gredo v pripust

    # age > 2: tukaj mladi biki postanejo cakajoci in cakajo v testu
    # po koncanem testu odberes pozitivno testirane PB
    # mladi biki postanejo čakajoči (~1 leto, da se osemeni krave s semenom oz. malo po 2. letu)
    if 'mladi' in categories.keys():
        #############################
        # to je za potrebe prehoda
        if kwargs.get("genTest_gpb") and 'mladi' in [i[0] for i in kwargs[
            'genotyped']]:  # če se greš genomsko shemo, jih takoj pogenotipiziraj in prestavi v gpb (5)
            ped.izberi_poEBV_top("M", kwargs.get('genpbn'), "mladi", "gpb", categories)
            ped.izloci_poEBV("M", (kwargs.get("mladin") - kwargs.get("genpbn")), "mladi", categories)
        ###########################
        else:
            ped.set_cat_old('mladi', 'cak', categories)  # mlade prestavi v cakajoce in jih izloci iz populacije
            ped.set_active_cat('mladi', 2, categories)

    # čakajočim bikov podaljšaj status (do starosti 5 let oz. kolikor let v testu)
    # hkrati jim tudi nastavi status izl
    # ped.set_cat_age_old(2, 'cak', 'cak', categories)
    if 'cak' in categories.keys():
        #############################
        # to je za potrebe prehoda
        if kwargs.get("genTest_gpb") and 'cak' in [i[0] for i in kwargs[
            'genotyped']]:  # če se greš genomsko shemo, jih takoj pogenotipiziraj in prestavi v gpb (5)
            ped.izberi_poEBV_top("M", kwargs.get('genpbn') * kwargs.get("cak"), "cak", "gpb", categories)
            ped.izloci_poEBV("M", ((kwargs.get("mladin") - kwargs.get('genpbn')) * kwargs.get("cak")), "cak",
                             categories)
        ##########################
        else:
            for i in range((2 + 1), (2 + kwargs.get(
                    'cak'))):  # 1 leto, ko začnejo semenit in so mladi biki, 3 so čakajoči, +1 da začneš prestavlajt
                ped.set_cat_age_old(i, 'cak', 'cak', categories)

    # če že imaš bike dovolj dolgo v testu, odberi pozitivno testirane bike
    if ('cak' in categories.keys()) and (
                (kwargs.get('cak') + 2) in ped.age()) and (
            kwargs.get("EBV") or kwargs.get("gpb_pb") or kwargs.get("genTest_mladi")) and 'cak' not in [i[0] for i in
                                                                                                        kwargs[
                                                                                                            'genotyped']]:  # +2 - eno leto so teleta, eno leto mladi biki, če imaš genomsko, gredo že pred do gpb
        ped.izberi_poEBV_top_age("M", (kwargs.get('cak') + 2), kwargs.get('pbn'), 'cak', 'pb', categories)
        ped.set_active_cat('cak', 2,
                           categories)  # tukaj moraš to nastaviti, zato ker fja izberi avtomatsko nastavi na active=1, vsi cakajoci so izloceni
        ped.izloci_poEBV_age("M", (kwargs.get('cak') + 2), (kwargs.get('mladin') - kwargs.get('pbn')), 'cak',
                             categories)  # TUKAJ MORA BITI ŠE STAROST!!!!!!!!!!!

    # genomski test: potomciNP = genomsko testiranje -> pozitivno testirani
    if kwargs.get('gEBV'):  # v prvem letu so vsi potomciNP v genomskem testiranju oz. pridobivanju gEBV
        ped.set_cat_sex_old('M', "potomciNP", "genTest", categories)
        ped.set_active_cat('potomciNP', 1, categories)

    if 'genTest' in categories.keys():  # če imaš genomsko testirane bike
        if kwargs.get('genTest_mladi'):
            ped.set_active_cat('genTest', 2, categories)
            ped.izberi_poEBV_top("M", kwargs.get('mladin'), 'genTest', 'mladi',
                                 categories)  # naslednji najboljši gredo v test
            ped.izberi_poEBV_OdDo("M", kwargs.get('mladin'), kwargs.get('vhlevljenin'), 'genTest', 'pripust1',
                                  categories)  # preostali genomsko testirani gredo v pripust
            ped.izloci_poEBV("M", (kwargs.get('potomciNPn') - kwargs.get('vhlevljenin')), "genTest", categories)

        if kwargs.get('genTest_gpb'):
            ped.izberi_poEBV_top("M", kwargs.get('genpbn'), "genTest", "gpb",
                                 categories)  # odberi genomsko testirane bike za AI
            ped.izberi_poEBV_OdDo("M", kwargs.get('genpbn'), kwargs.get('vhlevljenin'), 'genTest', 'pripust1',
                                  categories)  # preostali genomsko testirani gredo v pripust
            ped.izloci_poEBV("M", (kwargs.get('potomciNPn') - kwargs.get('vhlevljenin')), "genTest", categories)

    # prestavi jih naprej predno odbereš progeno testirane
    if 'gpb' in categories.keys():
        ped.set_cat_old('gpb', 'gpb', categories)
        ped.set_active_cat('gpb', 2, categories)

    # to je OK, ker ta vleče kategorije iz slovarja prejšnje generacije, ne direkt iz pedigreja
    # Tukaj glede na to, ali je označeno, da gredo gpb v pb, izvedi ta korak
    if kwargs.get('gpb_pb'):  # ali naj gredo gpb v progeno testirane
        if 'gpb' in categories.keys() and ((kwargs.get('cak') + 2) in ped.age()):  # če da, prestavi, ko so dovolj stari
            ped.set_cat_age_old((kwargs.get('cak') + 2), 'gpb', 'pb', categories)

    # pripust in pb so spet enaki pri obeh testiranjih
    # povprečna doba v pripustu - glede na to odberi bike, ki preživijo še eno leto
    if 'pripust1' in categories.keys():
        pripust2 = int(round(len(categories['pripust1']) * kwargs.get('pripust2n')))
        ped.izberi_random("M", pripust2, 'pripust1', 'pripust2',
                          categories)  # prestavi v 2. leto pripusta (ne vse - % glede na leta v UP)
        ped.izloci_random("M", (len(categories['pripust1']) - pripust2), 'pripust1',
                          categories)  # preostale iz pripusta izloci

    if 'pripust2' in categories.keys():  # izloci po 2. letu v pripustu
        ped.izloci_cat('pripust2', categories)

    # plemenske bike prestavljaj naprej - zvseskozi, za očete pa nato uporabi le tiste iz željenih let
    if 'pb' in categories.keys():
        ped.set_cat_old('pb', 'pb', categories)
        ped.set_active_cat('pb', 2, categories)

    print ped.cat()
    #########################################################
    # add new generation
    #########################################################
    # tukaj potem dodaj eno generacijo novorojenih
    ped.add_new_gen_naive(kwargs.get('stNBn'), kwargs.get('potomciNPn') * 2)
    # določi starost glede na te novorojene
    ped.compute_age()
    # dodaj matere
    ped.doloci_matere(kwargs.get('stNBn'), kwargs.get('potomciNPn'), kwargs.get('ptn'), kwargs.get('kraveUp'))
    # preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
    ped.mother_nr_blank()
    # dodaj očete
    ped.doloci_ocete(kwargs.get('stNBn'), kwargs.get('potomciNPn'), kwargs.get('cak'), kwargs.get('pbUp'),
                     kwargs.get('pripustDoz'), kwargs.get('mladiDoz'), kwargs.get('pozitivnoTestDoz'),
                     kwargs.get('CowsGenBulls_Per'), kwargs.get('EliteDamsPTBulls'), kwargs.get('EliteDamsGenBulls'),
                     kwargs.get('EliteDamsPABulls'),
                     kwargs.get('genTest_mladi'), kwargs.get('genTest_gpb'))

    # preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
    ped.mother_nr_blank()

    # ped.UpdateIndCat('/home/jana/')
    categories.clear()
    ped.save_cat_DF()
    ped.save_sex_DF()
    ped.save_active_DF()
    if kwargs.get('gEBV'):
        if kwargs.get('UpdateGenRef'):
            ped.updateAndSaveIndForGeno(kwargs.get('genotyped'), kwargs.get('NbUpdatedGen'), kwargs.get('sexToUpdate'),
                                        kwargs.get('AlphaSimDir'))
        if not kwargs.get('UpdateGenRef'):
            ped.saveIndForGeno(kwargs.get('genotyped'))
        os.system(
            'less IndForGeno.txt | wc -l > ReferenceSize_new.txt && cat ReferenceSize_new.txt ReferenceSize.txt > Reftmp && mv Reftmp ReferenceSize.txt')
    ped.write_ped(kwargs.get('AlphaSimDir') + "/ExternalPedigree.txt")
    ped.write_pedTotal(kwargs.get('AlphaSimDir') + "/ExternalPedigreeTotal.txt")
#    ped.write_pedTotal("/home/jana/PedTotal.txt")

    return ped, ped.save_cat(), ped.save_sex(), ped.save_active()

    #################


def odbira_mater(pedFile, **kwargs):
    print kwargs
    ped = pedigree(pedFile)

    # določi spol
    # ped.set_sex_AlphaSim(kwargs.get('AlphaSimDir'))

    # tukaj potem pridobi kategorije - če imaš samo eno burn-in, štartaš iz nule
    if max(ped.gen) == 1:  # ČE SAMO ENA GENERACIJA V PEDIGREJU! - to zdj ne delaš več tako
        ped.set_cat_gen(max(ped.gen), "nr")  # to je samo na prvem loopu
        ped.set_sex_list([x for x in range(0, ped.rows()) if x % 2 == 0], "F")
        ped.set_sex_list([x for x in range(0, ped.rows()) if x % 2 != 0], "M")
        # določi spol
        # ped.set_sex_AlphaSim(kwargs.get('AlphaSimDir'))
        ped.izberi_poEBV_top_catCurrent("F", int(kwargs.get('potomciNPn')), 'nr', 'potomciNP')
        ped.izberi_poEBV_top_catCurrent("M", int(kwargs.get('potomciNPn')), 'nr', 'potomciNP')

        # global categories #to moraš dat global samo v prvenm loopu, drugje dobiš return
        categories = ped.save_cat()
        # global sex
        sex = ped.save_sex()
        active = ped.save_active()
        ped = pedigree(pedFile)

    elif max(ped.gens()) > 1:
        categories = ped.create_categoriesDict('Categories_gen' + str(max(ped.gens())) + 'DF.csv')
        sex = ped.create_sexDict('Sex_gen' + str(max(ped.gens())) + 'DF.csv')
        active = ped.create_activeDict('Active_gen' + str(max(ped.gens())) + 'DF.csv')

        # določi spol po AlphaSIm - ampak to ni potrebno, ker je vseskozi external pedigre
    # ped.set_sex_AlphaSim(kwargs.get('AlphaSimDir'))

    ped.set_sex_prevGen(sex)  # add sex information for individuals from prevGen
    ped.set_active_prevGen(active)  # add active information for individuals from prevGen

    # remove category information from the ped itself
    for i in ped.gens():
        ped.set_cat_gen(i, "")

    # transfer culled (izlocene) category from prevGen
    ped.set_cat_old('izl', 'izl', categories)

    # compute age of the animals in the current selection year
    ped.compute_age()

    #################################################
    # FEMALES
    #################################################
    # age 0 - here you have newborn females (NB & potomkeNP) --> nekaj jih izloči, druge gredo naprej do ženskih telet
    ped.set_cat_sex_old("F", "potomciNP", "telF", categories)  # potomke načrtnih parjenj gredo v telice
    ped.set_active_catCurrent('telF', 1)
    izlF = int(kwargs.get('nrFn')) - int(kwargs.get('telFn'))  # number of culles NB females
    ped.izberi_poEBV_top("F", kwargs.get('telFn'), "nr", "telF",
                         categories)  # izberi NB ženske, ki preživijo in postanejo telice
    ped.izloci_poEBV("F", izlF, "nr", categories)  # cull females (lowest on EBV) tukaj jih izloči, funkcija v modulu

    # age 1 - pri enem letu osemeni določeno število telic (% določen zgoraj), druge izloči
    if 'telF' in categories.keys():
        ped.izberi_poEBV_top("F", kwargs.get('ptn'), 'telF', 'pt', categories)  # plemenske telice
        ped.izloci_poEBV("F", (len(categories['telF']) - kwargs.get('ptn')), 'telF', categories)  # preostale izloči

    # age > 2 - tukaj odbiraš in izločaš krave, odbiraš in izločaš BM
    # najprej dodaj nove krave, če jih že imaš v populaciji
    if ('pt' in categories.keys()):  # če imaš v pedigreju plemenske telice
        ped.set_cat_old('pt', 'k', categories)  # osemenjene telice postanejo krave - predpostavimo, da vse
        ped.set_active_cat('pt', 1, categories)
        ped.set_active_cat('k', 1, categories)

    # če imaš že dovolj stare krave, potem odberi BM - NAJPREJ - kasneje funkcija random polni prazne prostore (cat)!!!
    # BM se odbira po drugi laktaciji - to je starost 3 - 4 (starost v pedigreju = 3, ker imaš tudi 0)
    if ('k' in categories.keys()) and ((1 + kwargs.get('bmOdbira')) in ped.age()):
        ped.izberi_poEBV_top_age("F", kwargs.get('bmOdbira') + 1, int(kwargs.get('bmn') / kwargs.get('bmUp')), 'k',
                                 'pBM',
                                 categories)  # izberi BM, ki jih osemeniš (plemenske BM = pBM) iz krav po 2. laktaciji
        ped.set_active_cat('pBM', 1, categories)
    # ostale BM prestavi naprej - BM po 1. do izločitvene laktacije
    if 'pBM' in categories.keys():
        for i in range((1 + kwargs.get('bmOdbira') + 1), (
                        1 + kwargs.get('bmOdbira') + kwargs.get(
                    'bmUp'))):  # 1 leto prva osemenitev, bm odbrane po 2. laktaciji, +1 da začneš prestavljat
            ped.set_cat_age_old(i, 'pBM', 'pBM', categories)
        # spremeni kategorijo iz plemenskih BM v bm v zadnji laktaciji
        ped.set_cat_age_old((1 + kwargs.get('bmOdbira') + kwargs.get('bmUp')), 'pBM', 'bm',
                            categories)



        # potem odberi krave - na random, vsako leto za določeno število manj krav
    for i in range(2 + 1, (
                2 + kwargs.get('kraveUp'))):  # 2 + 1 - pri dveh letih prva laktacija, prestavljati začneš leto po tem
        ped.izberi_random_age('F', i, (
            kwargs.get('ptn') - (kwargs.get('MinusDamLact') * (i - 2)) - int(kwargs.get('bmn') / kwargs.get('bmUp'))),
                              'k',
                              'k', categories)
        ped.izloci_random_age('F', i, (kwargs.get('MinusDamLact')), 'k', categories)

    # potem izloči najstarejše krave - po 4. laktaciji
    if ('k' in categories.keys()) and ((kwargs.get('kraveUp') + 2) in ped.age()):  # izloči koliko laktacij + 2 leti
        ped.izloci_age_cat((kwargs.get('kraveUp') + 2), 'k', categories)

    # in izloči najastarejše BM, če jih imaš
    if ('bm' in categories.keys()):
        ped.izloci_cat('bm', categories)

    # tukaj določi, katere bodo matere
    return (
    ped, (ped.izberi_matere(kwargs.get('stNBn'), kwargs.get('potomciNPn'), kwargs.get('ptn'), kwargs.get('kraveUp'))))


# odberi testirane ocete
def odberi_testOce_gen(ped, **kwargs):
    print kwargs
    # ped = pedigree(pedFile)

    # določi spol
    # ped.set_sex_AlphaSim(kwargs.get('AlphaSimDir'))

    # tukaj potem pridobi kategorije - če imaš samo eno burn-in, štartaš iz nule
    if max(ped.gen) == 1:  # ČE SAMO ENA GENERACIJA V PEDIGREJU! - to zdj ne delaš več tako
        ped.set_cat_gen(max(ped.gen), "nr")  # to je samo na prvem loopu
        ped.set_sex_list([x for x in range(0, ped.rows()) if x % 2 == 0], "F")
        ped.set_sex_list([x for x in range(0, ped.rows()) if x % 2 != 0], "M")
        # določi spol
        # ped.set_sex_AlphaSim(kwargs.get('AlphaSimDir'))
        ped.izberi_poEBV_top_catCurrent("F", int(kwargs.get('potomciNPn')), 'nr', 'potomciNP')
        ped.izberi_poEBV_top_catCurrent("M", int(kwargs.get('potomciNPn')), 'nr', 'potomciNP')

        # global categories #to moraš dat global samo v prvenm loopu, drugje dobiš return
        categories = ped.save_cat()
        # global sex
        sex = ped.save_sex()
        active = ped.save_active()
        ped = pedigree(pedFile)

    elif max(ped.gens()) > 1:
        categories = ped.create_categoriesDict('Categories_gen' + str(max(ped.gens())) + 'DF.csv')
        sex = ped.create_sexDict('Sex_gen' + str(max(ped.gens())) + 'DF.csv')
        active = ped.create_activeDict('Active_gen' + str(max(ped.gens())) + 'DF.csv')

    #################################################################
    # MALES
    #################################################################
    # age 0: štartaš z NB in potomci NP --> odbereš vhlevljene iz potomcev NP in moška teleta in NB
    # NAJPREJ DELI, KI SO SKUPNI PROGENEMU IN GENOMSKEMU TESTIRANJU
    ped.izberi_random("M", kwargs.get('telMn'), "nr", "telM", categories)  # izberi moška teleta, ki preživijo (random)
    ped.izloci_random("M", int(kwargs.get('nrMn') - kwargs.get('potomciNPn') - kwargs.get('telMn')), "nr",
                      categories)  # druga teleta izloči

    if 'telM' in categories.keys():
        ped.izberi_random("M", kwargs.get('bik12n'), 'telM', 'bik12',
                          categories)  # random odberi bike, ki preživijo do 2. leta
        ped.izloci_random("M", (len(categories['telM']) - kwargs.get('bik12n')), 'telM', categories)  # izloči preostale

    if 'bik12' in categories.keys():  # izloci bike nad 2. leti
        ped.izloci_cat('bik12', categories)

    # Tukaj deli, kjer se progeni in genomski testi razlikujeta
    # progeni: potomciNP -> vhlevljeni -> mladi -> čakajoči -> pozitivno testirani
    # če je PROGENI TESTIRANJE
    if kwargs.get('EBV'):
        ped.izberi_random("M", kwargs.get('vhlevljenin'), "potomciNP", "vhlevljeni",
                          categories)  # vhlevi najboljše potomceNP - RANDOM
        ped.izloci_random('M', int(kwargs.get('potomciNPn') - kwargs.get('vhlevljenin')), 'potomciNP',
                          categories)  # druge potomceNP izloči


        # age1: tukaj odbereš mlade iz vhlevljenih bikov in bike, ki preživijo do drugega leta
    if 'vhlevljeni' in categories.keys():  # če imaš vhlevljene bike (samo v progenem testu)
        #############################
        # to je za potrebe prehoda
        ped.set_cat_old('vhlevljeni', 'kandidati', categories)

    # age > 2: tukaj mladi biki postanejo cakajoci in cakajo v testu
    # po koncanem testu odberes pozitivno testirane PB
    # mladi biki postanejo čakajoči (~1 leto, da se osemeni krave s semenom oz. malo po 2. letu)
    if 'mladi' in categories.keys():
        #############################
        ped.set_cat_old('mladi', 'kandidati', categories)


    # čakajočim bikov podaljšaj status (do starosti 5 let oz. kolikor let v testu)
    # hkrati jim tudi nastavi status izl
    # ped.set_cat_age_old(2, 'cak', 'cak', categories)
    if 'cak' in categories.keys():
        #############################
        # to je za potrebe prehoda
        ped.set_cat_old('cak', 'kandidati', categories)

    # genomski test: potomciNP = genomsko testiranje -> pozitivno testirani
    if kwargs.get('gEBV'):  # v prvem letu so vsi potomciNP v genomskem testiranju oz. pridobivanju gEBV
        ped.set_cat_sex_old('M', "potomciNP", "genTest", categories)
        ped.set_active_cat('potomciNP', 1, categories)

    if 'genTest' in categories.keys():  # če imaš genomsko testirane bike
        ped.set_cat_old("genTest", "kandidati", categories)

    # prestavi jih naprej predno odbereš progeno testirane
    if 'gpb' in categories.keys():
        ped.set_cat_old('gpb', 'gpb', categories)
        ped.set_active_cat('gpb', 2, categories)

    # to je OK, ker ta vleče kategorije iz slovarja prejšnje generacije, ne direkt iz pedigreja
    # Tukaj glede na to, ali je označeno, da gredo gpb v pb, izvedi ta korak
    if kwargs.get('gpb_pb'):  # ali naj gredo gpb v progeno testirane
        if 'gpb' in categories.keys() and ((kwargs.get('cak') + 2) in ped.age()):  # če da, prestavi, ko so dovolj stari
            ped.set_cat_age_old((kwargs.get('cak') + 2), 'gpb', 'pb', categories)

    # pripust in pb so spet enaki pri obeh testiranjih
    # povprečna doba v pripustu - glede na to odberi bike, ki preživijo še eno leto
    if 'pripust1' in categories.keys():
        pripust2 = int((len(categories['pripust1']) * kwargs.get('pripust2n')))
        ped.izberi_random("M", pripust2, 'pripust1', 'pripust2',
                          categories)  # prestavi v 2. leto pripusta (ne vse - % glede na leta v UP)
        ped.izloci_random("M", (len(categories['pripust1']) - pripust2), 'pripust1',
                          categories)  # preostale iz pripusta izloci

    if 'pripust2' in categories.keys():  # izloci po 2. letu v pripustu
        ped.izloci_cat('pripust2', categories)

    # plemenske bike prestavljaj naprej - zvseskozi, za očete pa nato uporabi le tiste iz željenih let
    if 'pb' in categories.keys():
        ped.set_cat_old('pb', 'pb', categories)
        ped.set_active_cat('pb', 2, categories)

    print ped.cat()
    return (ped, (ped.catCurrent_indiv('kandidati'),
                  ped.izberi_ocete_gen(kwargs.get('stNBn'), kwargs.get('potomciNPn'), kwargs.get('cak'), kwargs.get('pbUp'),
                                   kwargs.get('pripustDoz'), kwargs.get('mladiDoz'), kwargs.get('pozitivnoTestDoz'),
                                   kwargs.get('CowsGenBulls_Per'), kwargs.get('EliteDamsPTBulls'),
                                   kwargs.get('EliteDamsGenBulls'), kwargs.get('EliteDamsPABulls'),
                                   kwargs.get('genTest_mladi'), kwargs.get('genTest_gpb'))))

def odberiStarse_OCSgen(pedigree_genEBV, AlphaRelateDir, **selPar):
    # tukaj odberi metere
    ped, (bm, motherOther) = odbira_mater(pedigree_genEBV, **selPar)
    # tukaj odberi random sample krav
    motherSample = random.sample(motherOther, int(len(motherOther) * 0.2)) + random.sample(bm, int(len(bm) * 0.2))

    #tukaj odberi očete - kandidati so drugačni v prvem letu (mladi + vhlevljeni + čakajoči), kasneje pa so to genTest
    ped, (kandidati, (genOce)) = odberi_testOce_gen(ped, **selPar)

    pd.DataFrame({"ID": list(motherSample) + list(kandidati) + list(genOce)}).to_csv(AlphaRelateDir + "/IndOpt.txt", index=None, header=None)

    return ped

class AlphaRelate(object):
        def __init__(self, AlphaRelateDir, AlphaSimDir):
            self.AlphaRelateDir = AlphaRelateDir
            self.AlphaSimDir = AlphaSimDir
            self.AlphaRelateSpec = self.AlphaRelateDir + "/AlphaRelateSpec.txt"
            try:
		shutil.copy(AlphaSimDir + "/IndOpt.txt", AlphaRelateDir)
	    except:
		pass

        def preparePedigree(self):
            ped = pd.read_csv(self.AlphaSimDir + "/SimulatedData/PedigreeAndGeneticValues_cat.txt", sep="\s+")
            ped[["Indiv","Father", "Mother"]].to_csv(self.AlphaRelateDir + "/PEDIGREE.txt", sep=",", index=None, header=None)
            ped.loc[:, "sex1"] = [1 if x == "M" else 2 for x in ped.sex ]
            ped[["Indiv", "sex1"]].to_csv(self.AlphaRelateDir + "/GENDER.txt", sep=" ", index=None, header=None)

        def runAlphaRelate(self):
            os.chdir(self.AlphaRelateDir)
            os.system("./AlphaRelate_Linux")



class AlphaMate(object):
    def __init__(self, AlphaMateDir, AlphaSimDir, round):
        self.AlphaMateDir = AlphaMateDir
        self.AlphaSimDir = AlphaSimDir
        self.AlphaMateSpec_gen = AlphaMateDir + "/AlphaMateSpec_gen.txt"
        shutil.copy(self.AlphaMateSpec_gen, AlphaMateDir + "/AlphaMateSpec.txt")
        self.AlphaMateSpec = AlphaMateDir + "/AlphaMateSpec.txt"
	try:      
		shutil.copy(AlphaSimDir + "/IndOpt.txt", AlphaMateDir)
	except:
		pass
        self.indopt = sorted(list(pd.read_table("IndOpt.txt", header=None).loc[:, 0]))
        self.ped = pd.read_csv(self.AlphaSimDir + "/SimulatedData/PedigreeAndGeneticValues_cat.txt", sep="\s+")
#        self.round = max(self.ped.Generation)
	self.round = round




    def prepareGender(self):
        ped = self.ped
	ped.loc[:, "sex1"] = [1 if x == "M" else 2 for x in ped.sex]
        ped[ped.Indiv.isin(self.indopt)][["Indiv", "sex1"]].to_csv(self.AlphaMateDir + "/GENDER.txt", sep=" ", index=None, header=None)

    def countFemaleSel(self):
        gender = pd.read_table(self.AlphaMateDir + "/GENDER.txt", header=None, sep=" ")
        return (int(sum(gender[[1]] == 2)))

    def prepareCriterionFile(self):
        sol = pd.read_csv(self.AlphaSimDir + "/renumbered_Solutions_" + str(self.round), header=None, sep=" ")
        sol.columns = ["renID", "ID", "EBV"]
        sol.loc[sol.ID.isin(self.indopt)][["ID", "EBV"]].to_csv(self.AlphaMateDir + "/CRITERION.txt", header=None, index=None)

    def prepareSpecFile(self, NoMatings, NoMaleParents, NoFemaleParents, Degree):
        os.system('sed -i "s|NoMatings|' + str(NoMatings) + '|g" ' + self.AlphaMateSpec)
        os.system('sed -i "s|NoMaleParents|' + str(NoMaleParents) + '|g" ' + self.AlphaMateSpec)
        os.system('sed -i "s|NoFemaleParents|' + str(NoFemaleParents) + '|g" ' + self.AlphaMateSpec)
        os.system('sed -i "s|SetDegree|' + str(Degree) + '|g" ' + self.AlphaMateSpec)

    def runAlphaMate(self):
        os.chdir(self.AlphaMateDir)
        os.system("./AlphaMate")

    def obtainSelMales(self):
        Cont = pd.read_table(self.AlphaMateDir + "/ContributionsModeOptTarget1.txt", sep="\s+")
        selM = Cont.loc[Cont.Gender == 1]
        selM.loc[:, "nMatingPop"] = (selM.loc[:, "nMating"] / 0.2).astype(int)
        return list(chain.from_iterable([[int(father)] * dose for (father, dose) in zip(selM.Id, selM.nMatingPop)])) #to so očeti



def finishPedigree_OCS():
    # dodaj novo generacijo
    ped.add_new_gen_naive(selPar['stNBn'], selPar['potomciNPn'] * 2)
    ped.compute_age()
    # dodaj matere
    ped.doloci_matere(kwargs.get('stNBn'), kwargs.get('potomciNPn'), kwargs.get('ptn'), kwargs.get('kraveUp'))
    # preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
    ped.mother_nr_blank()
    # dodaj očete
    ped.ped.loc[ped.ped.cat.isin(['nr', 'potomciNP']), 'Father'] = Ocetje
    ped.ped.Father = ped.ped.Father.astype(int)
#    ped.write_ped("/home/jana/bin/AlphaMateLinux/OCSSloPop/CowSample//ExternalPedigree.txt")
#    ped.write_pedTotal("/home/jana/bin/AlphaMateLinux/OCSSloPop/CowSample//ExternalPedigreeTotal.txt")




# te funkcije so zdj redundant --> vse zdruzeno v selection_total
def selekcija_ena_gen(pedFile, categories=None, sex=None, active=None, stNB=None, nrFn=None, \
                      nrMn=None, telFn=None, telMn=None, potomciNPn=None, vhlevljenin=None, ptn=None, mladin=None,
                      bik12n=None, \
                      pripust1n=None, pripust2n=None, cak=None, kraveUp=None, bmOdbira=None, bmn=None, bmUp=None,
                      pripustDoz=None, mladiDoz=None, \
                      pozitivnoTestDoz=None, pbUp=None):
    ped = pedigree(pedFile)

    if max(ped.gen) == 1:
        ped.set_cat_gen(max(ped.gen), "nr")  # to je samo na prvem loopu
        ped.set_sex_list([x for x in range(0, ped.rows()) if x % 2 == 0], "F")
        ped.set_sex_list([x for x in range(0, ped.rows()) if x % 2 != 0], "M")
        ped.izberi_poEBV_top_catCurrent("F", int(potomciNPn), 'nr', 'potomciNP')
        ped.izberi_poEBV_top_catCurrent("M", int(potomciNPn), 'nr', 'potomciNP')

        # global categories #to moraš dat global samo v prvenm loopu, drugje dobiš return
        categories = ped.save_cat()
        # global sex
        sex = ped.save_sex()

        ped = pedigree(pedFile)

        ped.set_sex_prevGen(sex)  # prva odbira
        ped.compute_age()
        ped.select_age_0_1(categories, nrFn, nrMn, telFn, vhlevljenin, potomciNPn, telMn)
        ped.add_new_gen_naive(stNB, potomciNPn * 2)

        ped.compute_age()
        # dodaj matere
        ped.doloci_matere(stNB, ptn, kraveUp)
        # preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
        ped.mother_nr_blank()
        # dodaj očete

        ped.doloci_ocete(stNB, potomciNPn, cak, pripustDoz, mladiDoz, pbUp, pozitivnoTestDoz)
        # preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
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
        ped.select_age_0_1(categories, nrFn, nrMn, telFn, vhlevljenin, potomciNPn, telMn)
        ped.select_age_1_2(categories, ptn, mladin, vhlevljenin, bik12n)
        ped.add_new_gen_naive(stNB, potomciNPn * 2)
        ped.compute_age()

        # dodaj matere
        ped.doloci_matere(stNB, ptn, kraveUp)
        # preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
        ped.mother_nr_blank()
        # dodaj očete
        ped.doloci_ocete(stNB, potomciNPn, cak, pripustDoz, mladiDoz, pbUp, pozitivnoTestDoz)
        # preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
        ped.mother_nr_blank()

        categories.clear()  # sprazni slovar od prejšnjega leta

    # dodaj starše novorojenim - VEDNO PRVO MAME, KER JE FUNKCIJA ZA OČETE NAŠTIMANA, DA SE RAVNA PO MAMAH!


    if max(ped.gen) >= 3:
        ped.set_sex_prevGen(sex)
        ped.set_active_prevGen(active)

        for i in ped.gens():
            ped.set_cat_gen(i, "")

        ped.set_cat_old('izl', 'izl', categories)

        ped.compute_age()
        ped.select_age_0_1(categories, nrFn, nrMn, telFn, vhlevljenin, potomciNPn, telMn)
        ped.select_age_1_2(categories, ptn, mladin, vhlevljenin, bik12n)
        ped.select_age_2_3(categories, kraveUp, bmOdbira, bmUp, cak, pripust1n, pripust2n, mladin)
        ped.add_new_gen_naive(stNB, potomciNPn * 2)
        ped.compute_age()

        # dodaj matere
        ped.doloci_matere(stNB, ptn, kraveUp)
        # preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
        ped.mother_nr_blank()
        # dodaj očete
        ped.doloci_ocete(stNB, potomciNPn, cak, pripustDoz, mladiDoz, pbUp, pozitivnoTestDoz)
        # preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
        ped.mother_nr_blank()

        categories.clear()  # sprazni slovar od prejšnjega leta

    ped.write_ped(kwargs.get('AlphaSimDir') + "/ExternalPedigree.txt")

    return ped, ped.save_cat(), ped.save_sex(), ped.save_active()


def nastavi_cat(PedFile, **kwargs):
    ped = pedigree(PedFile)
    ped.compute_age()

    # določi spol
    ped.set_sex_AlphaSim(kwargs.get('AlphaSimDir'))

    # ped.set_sex_list([i for i in range(ped.rows()) if i%2==0], "F")
    # ped.set_sex_list([i for i in range(ped.rows()) if i%2!=0], "M")

    # MALES FIRST
    # najprej deli, ki so skupni progenemu testi in genomskemu testiranju - to je neselekcionirana populacija
    # določi moška teleta pod 12
    ped.izberi_random_age_naive('M', 0, kwargs.get('telMn'), 'telM')
    # določi bike nad 12 m
    ped.izberi_random_age_naive('M', 1, kwargs.get('bik12n'), 'bik12')

    # PROGENI TEST
    if kwargs.get('EBV'):
        # age 0
        ped.izberi_poEBV_top_age_naive('M', 0, kwargs.get('vhlevljenin'), 'vhlevljeni')  # določi vhlevljene
        # age1
        # določi mlade
        ped.izberi_poEBV_top_age_naive('M', 1, kwargs.get('mladin'), 'mladi')
        # določi pripust - 1. leto
        ped.izberi_poEBV_OdDo_age_naive('M', 1, kwargs.get('mladin'), kwargs.get('vhlevljenin'), 'pripust1')
        # age2,3,4
        for i in range(2, 2 + kwargs.get('cak')):  # leta, ko so cakajoci
            ped.izberi_poEBV_top_age_naive('M', i, kwargs.get('mladin'), 'cak')

        # od 1-2 leta v pripustu
        pripust2 = int(round(kwargs.get('pripust2n') * len(ped.catCurrent_indiv('pripust1'))))
        ped.izberi_poEBV_OdDo_age_naive('M', 2, kwargs.get('mladin'), (kwargs.get('mladin') + pripust2), 'pripust2')

        # age 5 - 10: pb
        pbAge = range((2 + kwargs.get('cak')), (2 + kwargs.get('cak') + kwargs.get('pbUp'))) if (2 + kwargs.get(
            'cak') + kwargs.get('pbUp')) <= max(ped.gens()) else range((2 + kwargs.get('cak')), max(ped.gens()))
        for i in pbAge:
            ped.izberi_poEBV_top_age_naive('M', i, kwargs.get('pbn'), 'pb')

    if kwargs.get('gEBV'):
        ped.izberi_poEBV_top_age_naive('M', 0, kwargs.get('potomciNPn'), 'genTest')

        # določi pripust - 1. leto
        ped.izberi_poEBV_OdDo_age_naive('M', 1, kwargs.get('genpbn'), kwargs.get('vhlevljenin'),
                                        'pripust1')  # kateri niso odbrani po genomskih, gredo za pripust

        # od 1-2 leta v pripustu
        pripust2 = int(round(kwargs.get('pripust2n') * len(ped.catCurrent_indiv('pripust1'))))
        ped.izberi_poEBV_OdDo_age_naive('M', 2, kwargs.get('genpbn'), (kwargs.get('genpbn') + pripust2), 'pripust2')

        for i in range(1, (
                    2 + kwargs.get(
                    'cak'))):  # od enega leta pa do konca progenega testa so genomsko testirani, kasneje progeno
            ped.izberi_poEBV_top_age_naive('M', i, kwargs.get('genpbn'), 'gpb')  # odberi genomsko testirane bike za AI

        # odberi tudi tiste, ki so že tudi progeno testirani
        pbAge = range((2 + kwargs.get('cak')), (2 + kwargs.get('cak') + kwargs.get('pbUp'))) \
            if (2 + kwargs.get('cak') + kwargs.get('pbUp')) <= max(ped.gens()) else range((2 + kwargs.get('cak')),
                                                                                          max(ped.gens()))
        for i in pbAge:
            ped.izberi_poEBV_top_age_naive('M', i, kwargs.get('genpbn'), 'pb')

    # FEMALES
    # age 0
    # določi ženska teleta pod 12
    ped.izberi_poEBV_top_age_naive('F', 0, kwargs.get('telFnTotal'), 'telF')

    # age1
    # določi plemenske telice
    ped.izberi_poEBV_top_age_naive('F', 1, kwargs.get('ptn'), 'pt')

    # age2 #pri starosti 2 let so samo plemenske telice, nima tudi BM
    for i in range(2, (1 + kwargs.get('bmOdbira'))):
        ped.izberi_poEBV_top_age_naive('F', i, kwargs.get('ptn'), 'k')

    # age3,4,5 #tukaj imaš enako število krav - minus tisti del, ki postane BM
    # odberi plemenske bm najprej
    for i in range((1 + kwargs.get('bmOdbira')), (1 + kwargs.get('bmOdbira') + kwargs.get('bmUp'))):
        ped.izberi_poEBV_top_age_naive('F', i, int(kwargs.get('bmn') / kwargs.get('bmUp')), 'pBM')
        ped.izberi_poEBV_top_age_naive('F', i, (
            (kwargs.get('ptn') - kwargs.get('MinusDamLact') * (i - 2)) - int(kwargs.get('bmn') / kwargs.get('bmUp'))),
                                       'k')

    # age 6
    # izberi odslužene bm
    ped.izberi_poEBV_top_age_naive('F', (1 + kwargs.get('bmOdbira') + kwargs.get('bmUp')),
                                   int(kwargs.get('bmn') / kwargs.get('bmUp')), 'bm')

    # ostali so izločeni

    ped.set_active_list(ped.row_cat('gpb'), 2)

    # določi še izločene
    ped.set_active_list(ped.row_cat(""), 2)
    ped.set_cat_list(ped.row_cat(""), 'izl')

    ped.add_new_gen_naive(kwargs.get('stNBn'), kwargs.get('potomciNPn') * 2)

    ped.compute_age()
    # dodaj matere
    ped.doloci_matere(kwargs.get('stNBn'), kwargs.get('potomciNPn'), kwargs.get('ptn'), kwargs.get('kraveUp'))
    # preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
    ped.mother_nr_blank()
    # dodaj očete
    ped.doloci_ocete(kwargs.get('stNBn'), kwargs.get('potomciNPn'), kwargs.get('cak'), kwargs.get('pbUp'),
                     kwargs.get('pripustDoz'), kwargs.get('mladiDoz'), kwargs.get('pozitivnoTestDoz'),
                     kwargs.get('CowsGenBulls_Per'), kwargs.get('EliteDamsPTBulls'), kwargs.get('EliteDamsGenBulls'),
                     kwargs.get('EliteDamsPABulls'),
                     kwargs.get('genTest_mladi'), kwargs.get('genTest_gpb'))

    # ped.UpdateIndCat('/home/jana/')
    ped.save_cat_DF()
    ped.save_sex_DF()
    ped.save_active_DF()
    if kwargs.get('gEBV'):
        ped.saveIndForGeno(kwargs.get('genotyped'))
    ped.write_ped(kwargs.get('AlphaSimDir') + "/ExternalPedigree.txt")
    ped.write_pedTotal(kwargs.get('AlphaSimDir') + "/ExternalPedigreeTotal.txt")
#    ped.write_pedTotal("/home/jana/PedTotal.txt")

    return ped, ped.save_cat(), ped.save_sex(), ped.save_active()


class TBVGenTable:  # to je tabela za grafiranje genetskih trendov čez populacije
    def __init__(self, TBVTable):
        self.TBVtable = pd.read_table(TBVTable, header=None, sep='\s+', names=['Indiv', 'TBV'])
        self.TBVmean = np.mean(self.TBVtable.TBV)
        self.TBVsd = np.std(self.TBVtable.TBV)
        self.TBVvar = np.var(self.TBVtable.TBV)
        self.TBVse = stats.sem(self.TBVtable.TBV)


class TBVPed(object):  # to je tabela za grafiranje genetskih trendov čez populacije
    def __init__(self, AlphaSimDir):
        self.AlphaSimDir = AlphaSimDir

    def genTrend(self, table, startgen, stopgen):
        TBVtable = pd.read_table(table, sep='\s+')
        TBVtable = TBVtable.loc[TBVtable.Generation.isin(range(startgen, stopgen))]
        gens = list(set(TBVtable.Generation))
        means = TBVtable.gvNormUnres1.groupby(TBVtable.Generation).aggregate(np.mean)
        vars = TBVtable.gvNormUnres1.groupby(TBVtable.Generation).aggregate(np.var)
        return gens, means, vars

    def catTrend(self):
        TBVtable = pd.read_table(self.AlphaSimDir + 'SimulatedData/PedigreeAndGeneticValues_cat.txt', sep='\s+')
        gen = max(set(TBVtable.Generation))
        means = pd.DataFrame({'Mean' + str(gen): TBVtable.gvNormUnres1.groupby(TBVtable.cat).aggregate(np.mean)})
        vars = pd.DataFrame({'Var' + str(gen): TBVtable.gvNormUnres1.groupby(TBVtable.cat).aggregate(np.var)})
        return pd.merge(means, vars, left_index=True, right_index=True)

    def GenTrend(self):
        TBVtable = pd.read_table(self.AlphaSimDir + 'SimulatedData/PedigreeAndGeneticValues_cat.txt', sep='\s+')
        gen = max(set(TBVtable.Generation))
        means = pd.DataFrame({'Mean' + str(gen): TBVtable.gvNormUnres1.groupby(TBVtable.Generation).aggregate(np.mean)})
        vars = pd.DataFrame({'Var' + str(gen): TBVtable.gvNormUnres1.groupby(TBVtable.Generation).aggregate(np.var)})
        return pd.merge(means, vars, left_index=True, right_index=True)


class TBVCat(TBVPed):
    def __init__(self, AlphaSimDir):
        super(TBVCat, self).__init__(AlphaSimDir)
        self.AlphaSimDir = AlphaSimDir
        # is the files exist - read them in - you need this in case you stop selection and restart it (from GUI)
        if os.path.isfile(self.AlphaSimDir + 'GenTrends_cat.csv'):
            self.genTrends_cat = pd.read_table(self.AlphaSimDir + 'GenTrends_cat.csv', sep=",")
        else:
            self.genTrends_cat = pd.DataFrame()

        if os.path.isfile(self.AlphaSimDir + 'GenTrends_gen.csv'):
            self.genTrends_gen = pd.read_table(self.AlphaSimDir + 'GenTrends_gen.csv', sep=",")
        else:
            self.genTrends_gen = pd.DataFrame()

    def saveTrends(self):
        currentTrend_C = self.catTrend()
        if not self.genTrends_cat.empty:
            self.genTrends_cat = pd.merge(self.genTrends_cat, currentTrend_C, left_index=True, right_index=True)
        if self.genTrends_cat.empty:
            self.genTrends_cat = pd.concat([self.genTrends_cat, currentTrend_C])

        currentTrend_G = self.GenTrend()
        if not self.genTrends_gen.empty:
            self.genTrends_gen = pd.merge(self.genTrends_gen, currentTrend_G, left_index=True, right_index=True)
        if self.genTrends_gen.empty:
            self.genTrends_gen = pd.concat([self.genTrends_gen, currentTrend_G])

    def writeTrends(self):
        self.genTrends_cat.to_csv(self.AlphaSimDir + 'GenTrends_cat.csv', sep=",", index=None)
        self.genTrends_gen.to_csv(self.AlphaSimDir + 'GenTrends_gen.csv', sep=",", index=None)

    def removeFiles(self):
        if os.path.isfile(self.AlphaSimDir + 'GenTrends_gen.csv'):
            os.remove(self.AlphaSimDir + 'GenTrends_gen.csv')
        if os.path.isfile(self.AlphaSimDir + 'GenTrends_cat.csv'):
            os.remove(self.AlphaSimDir + 'GenTrends_cat.csv')


class AlphaSimSpec:
    def __init__(self, AlphaSimDir, CodeDir):
        self.genSpecFile = CodeDir + "/AlphaSimSpec.txt"
        shutil.copy(self.genSpecFile, AlphaSimDir)
        self.SpecFile = AlphaSimDir + "/AlphaSimSpec.txt"

    def setPedType(self, pedType):
        os.system('sed -i "s|PedigreeType|' + pedType + '|g" ' + self.SpecFile)

    def setBurnInGen(self, StBurnInGen):
        os.system('sed -i "s|EnterBurnInGenerationNumber|' + str(StBurnInGen) + '|g" ' + self.SpecFile)

    def setSelGen(self, StSelGen):
        os.system('sed -i "s|EnterSelectionGenerationNumber|' + str(StSelGen) + '|g" ' + self.SpecFile)

    def setNoSires(self, NumberOfSires):
        os.system('sed -i "s|EnterNumberOfSires|' + str(NumberOfSires) + '|g" ' + self.SpecFile)

    def setNoDams(self, NumberOfDams):
        os.system('sed -i "s|EnterNumberOfDams|' + str(NumberOfDams) + '|g" ' + self.SpecFile)

    def turnOnGenFlex(self):
        os.system('sed -i "s|TurnOnGenFlex|On|g" ' + self.SpecFile)

    def turnOffGenFlex(self):
        os.system('sed -i "s|TurnOnGenFlex|Off|g" ' + self.SpecFile)

    def turnOnSelFlex(self):
        os.system('sed -i "s|TurnOnSelFlex|On|g" ' + self.SpecFile)

    def turnOffSelFlex(self):
        os.system('sed -i "s|TurnOnSelFlex|Off|g" ' + self.SpecFile)

    def setExtPedForGen(self, gen):
        os.system('sed -i "s|TheImportedGenerationPed|' + str(gen) + '|g" ' + self.SpecFile)

    def setTBVComp(self, option):
        os.system('sed -i "s|TBVComputation|' + str(option) + '|g" ' + self.SpecFile)

    def setFlexGenToFrom(self, To, From):
        os.system('sed -i "s|StartFlexGen,StopFlexGen|' + str(To) + ',' + str(From) + '|g" ' + self.SpecFile)

    def setNB(self, stNB):
        os.system('sed -i "s|EnterIndividualInPopulation|' + str(stNB) + '|g" ' + self.SpecFile)


class test:
    def __init__(self):
        print 123


class genInterval():
    def __init__(self, AlphaSimDir):
        self.name = AlphaSimDir + '/SimulatedData/PedigreeAndGeneticValues_cat.txt'
        self.pdPed = pd.read_table(self.name, sep='\s+')
        self.AlphaSimDir = AlphaSimDir

    def obtainSelInd_Parents(self, listCatOffspring):
        i = self.pdPed.loc[self.pdPed.cat.isin(listCatOffspring)][
            ['Generation', 'Indiv', 'cat', 'sex', 'Mother', 'Father']]
        i.columns = ['GenerationInd', 'IID', 'cat', 'sex', 'MID', 'FID']
        return i

    def obtainFathers(self, listCatOffspring):  # argument catOffspring is a list
        f = self.pdPed.loc[self.pdPed.Indiv.isin(self.pdPed.loc[self.pdPed.cat.isin(listCatOffspring)]['Father'])][
            ['Generation', 'Indiv']]
        f.columns = ['GenerationFather', 'FID']
        return f

    def obtainMothers(self, listCatOffspring):  # argument catOffspring is a list
        m = self.pdPed.loc[self.pdPed.Indiv.isin(self.pdPed.loc[self.pdPed.cat.isin(listCatOffspring)]['Mother'])][
            ['Generation', 'Indiv']]
        m.columns = ['GenerationMother', 'MID']
        return m

    def computeFatherAge(self, i, f):
        FA = pd.merge(i, f, on='FID', how='outer')
        FA.loc[:, 'FAge'] = FA.GenerationInd - FA.GenerationFather
        return FA

    def computeMotherAge(self, i, m):
        MA = pd.merge(i, m, on='MID', how='outer')
        MA.loc[:, 'MAge'] = MA.GenerationInd - MA.GenerationMother
        return MA

    def write_GenInt(self, FA, MA):
        FADF = FA.FAge.groupby([FA.GenerationInd, FA.sex]).mean()
        MADF = MA.MAge.groupby([MA.GenerationInd, MA.sex]).mean()
        FADF = FADF.to_frame().reset_index(level=['GenerationInd', "sex"])
        MADF = MADF.to_frame().reset_index(level=['GenerationInd', "sex"])
        genIntsF = pd.DataFrame({'Gen': FADF.GenerationInd, 'sex': FADF.sex, 'genInt': FADF.FAge, 'line': 'sire'})
        genIntsM = pd.DataFrame({'Gen': MADF.GenerationInd, 'sex': MADF.sex, 'genInt': MADF.MAge, 'line': 'dam'})
        genInts = pd.concat([genIntsM, genIntsF])
        if os.path.isfile(self.AlphaSimDir + 'GenInts.txt'):
            GenIntsOld = pd.read_csv(self.AlphaSimDir + 'GenInts.txt', sep=" ")
            pd.concat([GenIntsOld, genInts]).to_csv(self.AlphaSimDir + 'GenInts.txt', index=False, sep=" ")
        else:
            genInts.to_csv(self.AlphaSimDir + 'GenInts.txt', index=False, sep=" ")

    def prepareGenInts(self, listCatOffspring):
        i = self.obtainSelInd_Parents(listCatOffspring)
        f = self.obtainFathers(listCatOffspring)
        m = self.obtainMothers(listCatOffspring)
        FA = self.computeFatherAge(i, f)
        MA = self.computeMotherAge(i, m)
        self.write_GenInt(FA, MA)

    def plotGenInt(self):
        GenInts = pd.read_csv(self.AlphaSimDir + 'GenInts.txt', sep=" ")
        GenInts.loc[:, 'label'] = GenInts['line'] + GenInts['sex']
        Gens = GenInts[['Gen', 'genInt', 'label']]
        Gens.groupby(['Gen', 'label']).sum().unstack().plot()
        legend(['dam>dam', 'dam>sire', 'sire>dam', 'sire>sire'], loc='upper left')
        plt.savefig("GenInt_PLOT.png")


class snpFiles:
    def __init__(self, AlphaSimDir):
        self.AlphaSimDir = AlphaSimDir
        self.chipFile = self.AlphaSimDir + '/SimulatedData/AllIndividualsSnpChips/Chip1Genotype.txt'

    def createBlupf90SNPFile(self):
        if os.path.isfile(
                        self.AlphaSimDir + 'GenoFile.txt'):  # if GenoFile.txt exists, only add the newIndividuals for genotypisation
            os.system(
                'grep -Fwaf IndForGeno_new.txt ' + self.chipFile + ' > ChosenInd.txt')  # only individuals chosen for genotypisation - ONLY NEW - LAST GEN!
            os.system("sed 's/^ *//' ChosenInd.txt > ChipFile.txt")  # Remove blank spaces at the beginning
            os.system("cut -f1 -d ' ' ChipFile.txt > Individuals.txt")  # obtain IDs
            os.system('''awk '{$1=""; print $0}' ChipFile.txt | sed 's/ //g' > Snps.txt''')  # obtain SNP genotypes
            os.system(
                r'''paste Individuals.txt Snps.txt | awk '{printf "%- 10s %+ 15s\n",$1,$2}' > GenoFile_new.txt''')  # obtain SNP genotypes of the last generation
            os.system(
                'grep -Fwf IndForGeno.txt GenoFile.txt  > GenoFile_Oldtmp && mv GenoFile_Oldtmp GenoFile.txt')  # here obtain updated old reference - removed one generation
            os.system("cat GenoFile.txt GenoFile_new.txt > GenoFileTmp && mv GenoFileTmp GenoFile.txt")
            os.system("less GenoFile.txt | sort -n | uniq > Genotmp && mv Genotmp GenoFile.txt")
        else:  # else create a new GenoFile containg all the individuals in the IndForGeno.txt
            os.system(
                'grep -Fwaf IndForGeno.txt ' + self.chipFile + ' > ChosenInd.txt')  # only individuals chosen for genotypisation - ALL
            os.system("sed 's/^ *//' ChosenInd.txt > ChipFile.txt")  # Remove blank spaces at the beginning
            os.system("cut -f1 -d ' ' ChipFile.txt > Individuals.txt")  # obtain IDs
            os.system('''awk '{$1=""; print $0}' ChipFile.txt | sed 's/ //g' > Snps.txt''')  # obtain SNP genotypes
            os.system(
                r'''paste Individuals.txt Snps.txt | awk '{printf "%- 10s %+ 15s\n",$1,$2}' > GenoFile.txt''')  # obtain SNP genotypes of the last generation
        pd.read_csv(self.AlphaSimDir + '/SimulatedData/Chip1SnpInformation.txt', sep='\s+')[
            ["SnpId", "ChromId", "SnpSeqIdOnChrom"]].to_csv(
            self.AlphaSimDir + 'SnpMap.txt', index=None, sep=" ", header=None)
