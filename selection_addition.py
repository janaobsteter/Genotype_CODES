# -*- coding: utf-8 -*-
from __future__ import division
import os
import sys
import shutil
from collections import defaultdict
import pandas as pd
import numpy as np
from collections import defaultdict
import random
from itertools import chain
from subprocess import call

# reload(selection)
# from selection import *
# -*- coding: utf-8 -*-

#############################################33
# tukaj naštimaj parametre
# 1) odstotki živali v kateogorijah
# nT = 20000 # velikost populacije
# perF = 0.9
# perM = 0.1 # odstotki ženskih in moških živali v aktivni populaciji
# odstotki kategorij znotraj spola
stNB = 6700

nrF = stNB * 0.5
telF = 0.966
pt = 0.85  # koliko jih pride do telic
# kraveRemont = 0.25 # ali pa število laktacij
# TUKAJ ŠE; KOLIKO JIH PREŽIVI PRVO LAKTACIJO!
bm = 0.0127  # TO JE OD KRAV!!!
# sum([nrF, ptel, tel, pt24, t24, k])


nrM = stNB * 0.5
potomciNP = 0.0135  # %OD NOVOROJENIH

vhlevljeni = 0.6  # TO JE ODSTOTEK OD POTOMCEV NAČRNIH PARJENJ

mladi = 0.30  # TO JE ZNOTRAJ POTOMCEV NAČRTNIH PARJENJ
pripust = 0.70  # to je znotraj vhlevljenih
# pripust = 0.0135 #TO JE ZNOTRAJ BIKOV 12!!!!
telM = 0.73  # KOLIKO BIKOV postane moška teleta

# cakajoci = 0.018 #to je znotraj bikov12!!!!!
bik12 = 0.12  # koliko bikov do 12 mesecev preživi še do 2. leta

pb = 0.5  # KOLIKO OD MLADIH POSTANE TESTIRANIH
# nrM + bik12 +  bik24 

# 2)števila let v uporabi
kraveUp = 4  # povprečno koliko let so krave v populaciji (koliko laktacij) = remont krav
bmUp = 3  # koliko let so v uporabi BM - REMONT!
cak = 3  # koliko časa so mladi biki v testu oz. koliko časa so čakajoči
pbUp = 5  # koliko let so povrpečno v uporabi biki v AI
pripustUp = 1.4  # koliko let so v uporabi biki v pripustu
genomUp = 1.3  # koliko let so povprečno v uporabi genomsko testirani biki
bmOdbira = 2

##številke doz letno
pripustDoz = 15
pozitivnoTestDoz = 220
mladiDoz = 250
##################################################################################################################3
##################################################################################################################3
# od tu naprej samo delo s parametri = spremenvljivkami
##################################################################################################################3

# številke
# ženske
nrFn = int(stNB * 0.5)

telFn = int(telF * nrFn)
ptn = int(pt * telFn)
bmn = int(round(ptn * kraveUp * bm))  # to je od vseh krav
# kn = int(k * nF)
# bmn = int(bm * k * nF)
# sum([nrFn, pteln, teln, pt24n, t24n, kn])

nM = 0.5
nrMn = int(round(stNB * 0.5 * (1 - potomciNP)))
telMn = int(round(telM * nrMn))
bik12n = int(round(bik12 * telMn))
potomciNPn = int(round(potomciNP * nrMn))
vhlevljenin = int(round(potomciNPn * vhlevljeni))

mladin = int(round(mladi * vhlevljenin))  ###Mladi: mladin = int(nrM * nM * mladi)

# cakajocin = int(cakajoci* nM)
pripustTn = int(round((vhlevljenin - mladin) * pripustUp))
pripust1n = vhlevljenin - mladin
pripust2n = pripustTn - pripust1n
pbn = int(pb * mladin)
# sum([nrMn, bik12n, bik24n, cakajocin, pripustn, pbn, mladin])
# sum([nrMn, bik12n, bik24n])
##########################################################
# ženske
# odberi novorojene
# tu se ne rabiš obremenjevat s tem, katere so matere in kateri očetje - to narediš na koncu, ko konstruiraš novorojene za naslednjo gene
# ti novorojeni so že prva generacija
# na koncu potegni skupaj ID in kategorije ter to uporabi za določitev staršev#
# tudi šele tam se ukvarjaj z [pb]

import selection

reload(selection)
from selection import *


##############################################################################################3
##############################################################################################3
##############################################################################################3    

# združi celo selekcijo v eno funkcijo

def selekcija_total(pedFile,):
    ped = pedigree(pedFile)
    
    #tukaj potem pridobi kategorije - če imaš samo eno burn-in, štartaš iz nule
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
        active = ped.save_active()
        ped = pedigree(pedFile) 
 
    elif max(ped.gens()) > 1:
        categories = create_categoriesDict('Categories_gen' + str(max(self.gens())) + 'DF.csv')  
        sex = create_sexDict('Sex_gen' + str(max(self.gens())) + 'DF.csv')  
        active = create_activeDict('Active_gen' + str(max(self.gens())) + 'DF.csv')  
       
    ped.set_sex_prevGen(sex)  # add sex information for individuals from prevGen
    ped.set_active_prevGen(active)  # add sex information for individuals from prevGen

    #remove category information from the ped itself
    for i in ped.gens():
        ped.set_cat_gen(i, "")

    #transfer culled (izlocene) category from prevGen
    ped.set_cat_old('izl', 'izl', categories)

    #compute age of the animals in the current selection year
    ped.compute_age()
    
    #################################################
    # FEMALES
    #################################################
    # age 0 - here you have newborn females (NB & potomkeNP) --> nekaj jih izloči, druge gredo naprej do ženskih telet
    ped.set_cat_sex_old("F", "potomciNP", "telF", categories) #potomke načrtnih parjenj gredo v telice
    izlF = nrFn - telFn  # number of culles NB females
    ped.izberi_poEBV_top("F", (nrFn - izlF), "nr", "telF", categories)  # izberi NB ženske, ki preživijo in postanejo telice
    ped.izloci_poEBV("F", izlF, "nr", categories)  # cull females (lowest on EBV) tukaj jih izloči, funkcija v modulu

    # age 1 - pri enem letu osemeni določeno število telic (% določen zgoraj), druge izloči
    if 'telF' in categories.keys():
        ped.izberi_poEBV_top("F", ptn, 'telF', 'pt', categories) #plemenske telice
        ped.izloci_poEBV("F", (len(categories['telF']) - ptn), 'telF', categories)  #preostale izloči

    # age > 2 - tukaj odbiraš in izločaš krave, odbiraš in izločaš BM
    # najprej dodaj nove krave, če jih že imaš v populaciji
    if ('pt' in categories.keys()): #če imaš v pedigreju plemenske telice
        ped.set_cat_old('pt', 'k', categories)  # osemenjene telice postanejo krave - predpostavimo, da vse
    # krave po 1., 2., 3. laktaciji prestavi naprej v krave - OZIROMA PODALJŠAJ STATUS!
    for i in range(2 + 1, (2 + kraveUp)):  # 2 + 1 - pri dveh letih prva laktacija, prestavljati začneš leto po tem
        ped.set_cat_age_old(i, 'k', 'k', categories)
    # potem izloči najstarejše krave - po 4. laktaciji
    if ('k' in categories.keys()) and ((kraveUp + 2) in ped.age()):  # izloči koliko laktacij + 2 leti
        ped.izloci_age_cat((kraveUp + 2), 'k', categories)


    # če imaš že dovolj stare krave, potem odberi BM
    # BM se odbira po drugi laktaciji - to je starost 3 - 4 (starost v pedigreju = 3, ker imaš tudi 0)
    if ('k' in categories.keys()) and ((1 + bmOdbira) in ped.age()):
        ped.izberi_poEBV_top_age("F", bmOdbira+1, int(bmn / bmUp), 'k', 'pBM', categories)  # izberi BM, ki jih osemeniš (plemenske BM = pBM) iz krav po 2. laktaciji
    # in izloči najastarejše BM, če jih imaš
    if ('bm' in categories.keys()):
        ped.izloci_cat('bm', categories)
    # ostale BM prestavi naprej - BM po 1. do izločitvene laktacije
    if 'pBM' in categories.keys():
        for i in range((1 + bmOdbira + 1), (
                1 + bmOdbira + bmUp)):  # 1 leto prva osemenitev, bm odbrane po 2. laktaciji, +1 da začneš prestavljat
            ped.set_cat_age_old(i, 'pBM', 'pBM', categories)
        # spremeni kategorijo iz plemenskih BM v bm v zadnji laktaciji 
        ped.set_cat_age_old((1 + bmOdbira + bmUp), 'pBM', 'bm',
                            categories)  

    #################################################################
    # MALES
    #################################################################
    # age 0: štartaš z NB in potomci NP --> odbereš vhlevljene iz potomcev NP in moška teleta in NB
    ped.izberi_poEBV_top("M", vhlevljenin, "potomciNP", "vhlevljeni",
                         categories)  # vhlevi najboljše potomceNP
    ped.izloci_poEBV("M", int(potomciNPn - vhlevljenin), 'potomciNP', categories) #druge potomceNP izloči
    ped.izberi_random("M", telMn, "nr", "telM", categories) #izberi moška teleta, ki preživijo (random)
    ped.izloci_random("M", int(nrMn - telMn), "nr", categories) #druga teleta izloči

    # age1: tukaj odbereš mlade iz vhlevljenih bikov in bike, ki preživijo do drugega leta
    if 'vhlevljeni' in categories.keys():
        ped.izberi_poEBV_top("M", mladin, "vhlevljeni", "mladi", categories)  # odberi mlade
        ped.izberi_poEBV_OdDo("M", mladin, vhlevljenin, "vhlevljeni", "pripust1", categories)  # preostali vhlevljeni gredo v pripust
    if 'telM' in categories.keys():
        ped.izberi_random("M", bik12n, 'telM', 'bik12', categories) #random odberi bike, ki preživijo do 2. leta
        ped.izloci_random("M", (len(categories['telM']) - bik12n), 'telM', categories) #izloči preostale

    # age > 2: tukaj mladi biki postanejo cakajoci in cakajo v testu
    #po koncanem testu odberes pozitivno testirane PB
    # mladi biki postanejo čakajoči (~1 leto, da se osemeni krave s semenom oz. malo po 2. letu)
    if 'mladi' in categories.keys():
        ped.set_cat_old('mladi', 'cak', categories) #mlade prestavi v cakajoce in jih izloci iz populacije
        ped.set_active_cat('mladi', 2, categories)

    if 'bik12' in categories.keys(): #izloci bike nad 2. leti
        ped.izloci_cat('bik12', categories)
    
    # povprečna doba v pripustu - glede na to odberi bike, ki preživijo še eno leto
    if 'pripust1' in categories.keys():
        ped.izberi_random("M", pripust2n, 'pripust1', 'pripust2', categories) #prestavi v 2. leto pripusta (ne vse - % glede na leta v UP)
        ped.izloci_random("M", (pripust1n - pripust2n), 'pripust1', categories) #preostale iz pripusta izloci

    if 'pripust2' in categories.keys(): #izloci po 2. letu v pripustu
        ped.izloci_cat('pripust2', categories)

    # čakajočim bikov podaljšaj status (do starosti 5 let oz. kolikor let v testu)
    # hkrati jim tudi nastavi status izl
    # ped.set_cat_age_old(2, 'cak', 'cak', categories)
    if 'cak' in categories.keys():
        for i in range((2 + 1), (2 + cak)):  # 1 leto, ko začnejo semenit in so mladi biki, 3 so čakajoči, +1 da začneš prestavlajt
            ped.set_cat_age_old(i, 'cak', 'cak', categories)



    # če že imaš bike dovolj dolgo v testu, odberi pozitivno testirane bike
    if ('cak' in categories.keys()) and ((cak + 2) in ped.age()):  # +2 - eno leto so teleta, eno leto mladi biki
        ped.izberi_poEBV_top_age("M", (cak + 2), int(mladin * 0.5), 'cak', 'pb', categories)
        ped.set_active_cat('cak', 2,
                           categories)  # tukaj moraš to nastaviti, zato ker fja izberi avtomatsko nastavi na active=1
        ped.izloci_poEBV_age("M", (cak + 2), int(mladin * 0.5), 'cak',
                             categories)  # TUKAJ MORA BITI ŠE STAROST!!!!!!!!!!!

    # plemenske bike prestavljaj naprej
    if 'pb' in categories.keys():
        ped.set_cat_old('pb', 'pb', categories)
    
    #########################################################
    #add new generation
    #########################################################
    #tukaj potem dodaj eno generacijo novorojenih    
    ped.add_new_gen_naive(stNB, potomciNPn*2)
    #določi starost glede na te novorojene
    ped.compute_age()
    #dodaj matere
    ped.doloci_matere(stNB, ptn, kraveUp)
    #preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
    ped.mother_nr_blank()
    #dodaj očete
    ped.doloci_ocete(stNB, potomciNPn, cak, pripustDoz, pbUp, mladiDoz, pozitivnoTestDoz)

    #preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
    ped.mother_nr_blank()   

    categories.clear()
    ped.save_cat_DF()
    ped.save_sex_DF()
    ped.save_active_DF()
    ped.write_ped("/home/jana/bin/AlphaSim1.05Linux/ExternalPedigree.txt")
    
    return ped, ped.save_cat(), ped.save_sex(), ped.save_active()



