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
from selection10 import *
#reload(selection)
#from selection import *
# -*- coding: utf-8 -*-

#############################################33
#tukaj naštimaj parametre
#1) odstotki živali v kateogorijah
#nT = 20000 # velikost populacije
#perF = 0.9
#perM = 0.1 # odstotki ženskih in moških živali v aktivni populaciji
#odstotki kategorij znotraj spola
stNB = 6700

nrF = stNB * 0.5
telF = 0.966
pt = 0.85 #koliko jih pride do telic
#kraveRemont = 0.25 # ali pa število laktacij
#TUKAJ ŠE; KOLIKO JIH PREŽIVI PRVO LAKTACIJO!
bm = 0.0127 #TO JE OD KRAV!!!
#sum([nrF, ptel, tel, pt24, t24, k])


nrM = stNB * 0.5
potomciNP = 0.0135 #%OD NOVOROJENIH

vhlevljeni = 0.6 #TO JE ODSTOTEK OD POTOMCEV NAČRNIH PARJENJ

mladi = 0.30 #TO JE ZNOTRAJ POTOMCEV NAČRTNIH PARJENJ
pripust = 0.70 # to je znotraj vhlevljenih
#pripust = 0.0135 #TO JE ZNOTRAJ BIKOV 12!!!!
telM = 0.73 #KOLIKO BIKOV postane moška teleta

#cakajoci = 0.018 #to je znotraj bikov12!!!!!
bik12 = 0.12 #koliko bikov do 12 mesecev preživi še do 2. leta


pb = 0.5 #KOLIKO OD MLADIH POSTANE TESTIRANIH
#nrM + bik12 +  bik24 

#2)števila let v uporabi
kraveUp = 4 #povprečno koliko let so krave v populaciji (koliko laktacij) = remont krav
bmUp = 3 # koliko let so v uporabi BM - REMONT!
cak = 3 #koliko časa so mladi biki v testu oz. koliko časa so čakajoči
pbUp = 5 #koliko let so povrpečno v uporabi biki v AI
pripustUp = 1.4 # koliko let so v uporabi biki v pripustu
genomUp = 1.3 #koliko let so povprečno v uporabi genomsko testirani biki
bmOdbira = 2



##številke doz letno
pripustDoz = 15
pozitivnoTestDoz = 220
mladiDoz = 250
##################################################################################################################3
##################################################################################################################3
#od tu naprej samo delo s parametri = spremenvljivkami
##################################################################################################################3

#številke
#ženske
nrFn = int(stNB * 0.5)

telFn = int(telF * nrFn)
ptn = int(pt* telFn)
bmn = int(round(ptn * kraveUp *bm)) #to je od vseh krav
#kn = int(k * nF)
#bmn = int(bm * k * nF)
#sum([nrFn, pteln, teln, pt24n, t24n, kn])

nM = 0.5
nrMn = int(round(stNB * 0.5 * (1 - potomciNP))) 
telMn = int(round(telM * nrMn))
bik12n = int(round(bik12 * telMn))
potomciNPn = int(round(potomciNP * nrMn))
vhlevljenin  = int(round(potomciNPn * vhlevljeni))

mladin = int(round(mladi * vhlevljenin)) ###Mladi: mladin = int(nrM * nM * mladi)

#cakajocin = int(cakajoci* nM)
pripustTn = int(round((vhlevljenin - mladin)*pripustUp))
pripust1n = vhlevljenin - mladin
pripust2n = pripustTn - pripust1n
pbn = int(pb * mladin)
#sum([nrMn, bik12n, bik24n, cakajocin, pripustn, pbn, mladin])
#sum([nrMn, bik12n, bik24n])
##########################################################
#ženske
#odberi novorojene
#tu se ne rabiš obremenjevat s tem, katere so matere in kateri očetje - to narediš na koncu, ko konstruiraš novorojene za naslednjo gene
#ti novorojeni so že prva generacija
#na koncu potegni skupaj ID in kategorije ter to uporabi za določitev staršev#
#tudi šele tam se ukvarjaj z [pb]

import selection
reload(selection)
from selection import *



     



##############################################################################################3
##############################################################################################3
##############################################################################################3    





#VEDNO NAJPREJ IZLOČI /ODBERI PO PV!!! - funckije za odbiro na random imajo pogoj, da je kateogrija prosta
def select_age_0_1(ped, categories): #tukaj odbereš iz novorojenih živali tel, ptel in mlade bike, pripust1
    #FEMALES
    ped.set_cat_sex_old("F", "potomciNP", "telF", categories)
    izlF = nrFn - telFn#koliko jih izločiš
    ped.izloci_poEBV("F", izlF, "nr", categories) #tukaj jih izloči, funkcija v modulu

    ped.izberi_poEBV_top("F", (nrFn - izlF), "nr", "telF", categories) #izberi telice, ki jih osemeniš --> krave
    
    
    #MALES
    ped.izberi_poEBV_top( "M", vhlevljenin, "potomciNP", "vhlevljeni", categories) #odberi mlade TO SAMO NA ZAČETKU; POTEM POTOMCI BM IN ELITE!
    ped.izloci_poEBV("M", int(potomciNPn - vhlevljenin), 'potomciNP', categories)
    ped.izberi_random("M", telMn, "nr", "telM", categories)
    ped.izloci_random("M", int(nrMn - telMn), "nr", categories)
    

def select_age_1_2(ped, categories): # tukaj odbereš nič pri kravah - razen, če so že bikovske matere, pripust 2, bike24
    #FEMALES
    ped.izberi_poEBV_top("F", ptn, 'telF', 'pt', categories)
    ped.izloci_poEBV("F", (len(categories['telF']) - ptn), 'telF', categories) #terlice postanejo
   
    
    #MALES
    ped.izberi_poEBV_top( "M", mladin, "vhlevljeni", "mladi", categories) #odberi mlade
    ped.izberi_poEBV_OdDo( "M", mladin, vhlevljenin, "vhlevljeni", "pripust1", categories) #odberi v pripustu
    ped.izberi_random( "M", bik12n, 'telM', 'bik12', categories)
    ped.izloci_random( "M", (len(categories['telM']) - bik12n), 'telM', categories)



#tukaj lahk daš vse v eno funkcijo - variabilno - koliko let krave, koliko let v testu
def select_age_2_3(ped, categories):
    #FEMALES
    #najprej dodaj nove krave
    ped.set_cat_old('pt', 'k', categories) #osemenjene telice postanejo krave - predpostavimo, da vse
    #potem izloči najstarejše krave - po 4. laktaciji
    if ('k' in categories.keys()) and ((kraveUp+2) in ped.age()): #izloči koliko laktacij + 2 leti
        ped.izloci_age_cat((kraveUp+2), 'k', categories)
    #ostale krave prestavi naprej v krave - OZIROMA PODALJŠAJ STATUS!
    ped.set_cat_age_old(3, 'k', 'k', categories)
    ped.set_cat_age_old(4, 'k', 'k', categories) 
    ped.set_cat_age_old(5, 'k', 'k', categories)    
              
    #če imaš že dovolj stare krave, potem odberi BM
    #BM se odbira po drugi laktaciji - to je starost 3 - 4 (starost v pedigreju = 3, ker imaš tudi 0)
    if ('k' in categories.keys()) and ((1 + bmOdbira) in ped.age()):
        ped.izberi_poEBV_top_age("F",3, int(bmn /bmUp), 'k', 'pBM', categories) #izberi bikovske matere
    #in izloči najastarejše BM, če jih imaš
    if ('bm' in categories.keys()):
        ped.izloci_cat('bm', categories)
    #ostale BM prestavi naprej
    for i in range((1 + bmOdbira + 1), (1 + bmOdbira + bmUp)): #1 leto prva osemenitev, bm odbrane po 2. laktaciji, +1 da začneš prestavljat
        ped.set_cat_age_old(i, 'pBM', 'pBM', categories)
    ped.set_cat_age_old((1 + bmOdbira + bmUp), 'pBM', 'bm', categories) #spremeni kategorijo iz plemenskih BM v bm v zadnji laktaciji
    
    #MALES
    #mladi biki postanejo čakajoči (~1 leto, da se osemeni krave s semenom oz. malo po 2. letu)
    ped.set_cat_old('mladi', 'cak', categories) 
    ped.set_active_cat('mladi', 2, categories)
    
    #čakajočim bikov podaljšaj status (do starosti 5 let)
    #hkrati jim tudi nastavi status izl
    #ped.set_cat_age_old(2, 'cak', 'cak', categories)
    for i in range((2 + 1), (2 + cak)): #1 leto, ko začnejo semenit in so mladi biki, 3 so čakajoči, +1 da začneš prestavlajt
        ped.set_cat_age_old(i, 'cak', 'cak', categories)
    
    #povprečna doba v pripustu - glede na to odberi bike, ki preživijo še eno leto
    if 'pripust1' in categories.keys():
        ped.izberi_random( "M", pripust2n, 'pripust1', 'pripust2', categories)
        ped.izloci_random( "M", (pripust1n - pripust2n), 'pripust1', categories)

    #plemenske bike prestavljaj naprej
    ped.set_cat_old('pb', 'pb', categories)
    ped.izloci_cat('bik12', categories)
    ped.izloci_cat('pripust2', categories)
    if ('cak' in categories.keys()) and ((cak+2) in ped.age()): #+2 - eno leto so teleta, eno leto mladi biki
        ped.izberi_poEBV_top_age("M", (cak +2), int(mladin * 0.5), 'cak', 'pb', categories)
        ped.set_active_cat('cak', 2, categories) #tukaj moraš to nastaviti, zato ker fja izberi avtomatsko nastavi na active=1
        ped.izloci_poEBV_age("M",(cak+2), int(mladin * 0.5), 'cak', categories) #TUKAJ MORA BITI ŠE STAROST!!!!!!!!!!!



def doloci_matere(ped):
    #MATERE
    sTbmMother = 90 if len(ped.catCurrent_indiv('pBM')) >= 90 else len(ped.catCurrent_indiv('pBM'))
    if sTbmMother != 0:
        bmMother = ped.select_mother_random('pBM', sTbmMother)
        ped.set_mother_catPotomca(bmMother, 'potomciNP')
    #
    
    if 'k' in ped.cat():#TUKAJ SO DOLOČENE SEDAJ VSE MATERE!!!
        mother = ped.select_mother_EBV_top('k', int(round(11000*0.7))) #tukaj odberi brez tistih, ki so za gospodarsko križanje
        if len(mother) >= (stNB - sTbmMother): # če že imaš dovolj krav, določi matere vsem novorojenim oz. odbiraš matere, saj jih imaš preveč!
            motherOther = random.sample(mother, (stNB - sTbmMother))
            ped.set_mother_catPotomca(motherOther, 'nr') #TUKAJ SO DOLOČENE SEDAJ VSE MATERE!!!
        elif len(mother) < (stNB - sTbmMother): # če jih še ni dovolj, ne odbiraš mater, ampak uporabiš vse MINUS gosp. križanmja
            ped.set_mother_catPotomca(mother, 'nr') 


def doloci_ocete(ped):
#OČETJE
    mladiOce = ped.catCurrent_indiv('mladi')
    pripustOce = ped.catCurrent_indiv('pripust1') + ped.catCurrent_indiv('pripust2') 
    testiraniOce = list(chain.from_iterable([ped.catCurrent_indiv_age('pb', (2 + cak + x)) for x in range(1, pbUp+1)])) # v času, ko določaš potomce, so že eno leto starjši!!!
    bmMother = 90 if len(ped.catCurrent_indiv('pBM')) >= 90 else len(ped.catCurrent_indiv('pBM'))
    if 'pb' in ped.cat():
        elita = np.random.choice(ped.catCurrent_indiv('pb'), bmMother, replace=True) #navidezna elita
#        pd.Series(elita).value_counts()#preveri zastopanost po bikih
        #naštimaj očete elite --> BM
        ped.set_father_catPotomca(elita, 'potomciNP')    

    ocetje = pripustOce*pripustDoz + testiraniOce*pozitivnoTestDoz + mladiOce*mladiDoz
    if len(ocetje) >= (stNB - potomciNPn*2): #če imaš dovolj DOZ za vse NB
        ocetjeNB = random.sample(ocetje, (stNB - potomciNPn*2)) #tukaj izbereš očete za vse krave  - razen BM!
        ped.set_father_catPotomca(ocetjeNB, 'nr')
    if len(ocetje) < (stNB - potomciNPn*2):
        ped.set_father_catPotomca(ocetje, 'nr')







#####################################################################
#tukaj je zdj funkcija, ki vse to dela!
##################################################################### 

#to je funkcija za odbiro in določanje staršev
#prvi pogoj if max gen = 1 je za primer, ko štartaš s praznim naivnim pedigrejem brez staršev - mam in očetov ni v pedigreju
#drugi pogoj,ko dodaš generacijo novorojenih in pelješ prejšnjo generacijo naprej
#tretji krog so združene vse selekcijske odločitve po tem - počasi dobiš bm in pb, če jih ni, se pač ti starši ne določajo
def selekcija_ena_gen(pedFile, categories = None, sex = None, active = None):
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
        select_age_0_1(ped, categories)
        ped.add_new_gen_naive(stNB, potomciNPn*2)
        
        ped.compute_age()
        #dodaj matere
        doloci_matere(ped)
        #preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
        ped.mother_nr_blank()
        #dodaj očete
        doloci_ocete(ped)
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
        select_age_0_1(ped, categories)
        select_age_1_2(ped, categories)

        ped.add_new_gen_naive(stNB, potomciNPn*2)
        ped.compute_age()
        
        #dodaj matere
        doloci_matere(ped)
        #preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
        ped.mother_nr_blank()
        #dodaj očete
        doloci_ocete(ped)
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
        select_age_0_1(ped, categories)
        select_age_1_2(ped, categories)

        select_age_2_3(ped, categories)

        ped.add_new_gen_naive(stNB, potomciNPn*2)  
        ped.compute_age()
        
        #dodaj matere
        doloci_matere(ped)
        #preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
        ped.mother_nr_blank()
        #dodaj očete
        doloci_ocete(ped)
        #preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
        ped.mother_nr_blank()
                
        categories.clear() #sprazni slovar od prejšnjega leta
        ped.write_ped("/home/jana/bin/AlphaSim1.05Linux/ExternalPedigree.txt")

    return ped, ped.save_cat(), ped.save_sex(), ped.save_active()


#######################################################
#TO JE, ČE ŠTARTAŠ S POLNO AKTIVNO POPULACIJO IN DOLOČIŠ KATEGORIJE
#######################################################
def nastavi_cat(PedFile):
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
    doloci_matere(ped)
    #preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
    ped.mother_nr_blank()
    #dodaj očete
    doloci_ocete(ped)
    #preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
    ped.mother_nr_blank()   
    
    ped.write_ped("/home/jana/bin/AlphaSim1.05Linux/ExternalPedigree.txt")
    
    return ped, ped.save_cat(), ped.save_sex(), ped.save_active()

###########################################################################
########################################################################

#ped = pedigree("~/Documents/PhD/Simulaton/Pedigrees/PedPython.txt")

#TUKAJ PA JE SEDAJ PROGRAM
#Najprej določi, ali štartaš od začetka in počasi polniš populacijo ali štartaš z polnim pedigrejem 
OPTION = raw_input("1 - Polnjenje populacije; 2 - Start z polnim pedigrejem ")
#PedFile = raw_input("Vnesi pot do pedigreja")
StBurnInGen = input("Vnesi stevilo burn in generacij: ")
StSelGen = input("Vnesi stevilo krogov oz. generacij: ")
AlphaSimDir = '/home/jana/bin/AlphaSim1.05Linux'
AlphaSimPed = raw_input("Vnesi pot do output AlphaSim pedigrejev im ime file")
AlphaSimPed = "/home/jana/Documents/PhD/Simulaton/Pedigrees/Pedigree_10burnIn_10gen.txt"
AlphaSimPed = '/home/jana/bin/AlphaSim1.05Linux/SimulatedData/PedigreeAndGeneticValues.txt'

"""
if OPTION == 1:
    for krog in StKrogov:
        #PRERAČUNAŠ EBV v Ru in ZAPIŠEŠ PEDIGRE
        shutil.copy ("Rcorr_PedEBV.R", "Rcorr_PedEBV_ThisGen.R")
        os.system('sed -i "s|AlphaSimPed|' + AlphaSimPed + '|g" Rcorr_PedEBV_ThisGen.R')
        call('Rscript Rcorr_PedEBV_ThisGen.R', shell=True)
        selekcija_ena_gen('GenPed_EBV.txt') #to ti določi kategorije in starše
        #prestavi se v AlphaSim Dir
        os.chdir(AlphaSimDir)
        #TUKAJ POTEM POPRAVIŠ AlphaSimSpec
        #tukaj poženeš prvič po burn inu
        os.system('sed -i "s|StartStopGeneration                               ,1,' + str(StBurnInGen) + '|StartStopGeneration                               ,' + str(StBurnInGen+1) + ',' + str(StBurnInGen+1) + '|g" Rcorr_PedEBV_ThisGen.R')
        os.system('sed -i "s|Internal|ExternalPedigree_NextGen.txt|g" AlphaSimSpec.txt')
        #POŽENEŠ ALPHASIM        
        os.system('./AlphaSim1.05')
"""

 
if OPTION == 2:
    BurnInYN = raw_input("Do you already have a burn in population? [Y/N] ")
    if BurnInYN == 'N':  
        for roundNo in range(StSelGen+1):
            if roundNo == 0: #do burn in
                #prestavi se v AlphaSim Dir
                os.chdir(AlphaSimDir)
                shutil.copy('/home/jana/Genotipi/Genotipi_CODES/AlphaSimSpec.txt', AlphaSimDir)
                os.system('sed -i "s|PedigreeType|Internal|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|EnterBurnInGenerationNumber|' + str(StBurnInGen) + '|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|EnterSelectionGenerationNumber|' + str(StSelGen) + '|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|EnterNumberOfSires|' + str(NumberOfSires) + '|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|EnterNumberOfDams|' + str(NumberOfDams) + '|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|TurnOnGenFlex|On|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|StartFlexGen,StopFlexGen|1,' + str(StBurnInGen + 1) + '|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|TurnOnSelFlex|On|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|TheImportedGenerationPed|' +  str(StBurnInGen + 1) + '|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|TBVComputation|1|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|EnterIndividualInPopulation|' +str(stNB)+ '|g" AlphaSimSpec.txt') 
                #POŽENEŠ ALPHASIM        
                os.system('./AlphaSim1.05')
    
    
            elif roundNo == 1:
                os.chdir('/home/jana/Genotipi/Genotipi_CODES/')
                #PRERAČUNAŠ EBV v Ru in ZAPIŠEŠ PEDIGRE
                shutil.copy ("Rcorr_PedEBV.R", "Rcorr_PedEBV_ThisGen.R")
                os.system('sed -i "s|AlphaSimPed|' + AlphaSimPed + '|g" Rcorr_PedEBV_ThisGen.R')
                call('Rscript Rcorr_PedEBV_ThisGen.R', shell=True)
                #tukaj nastvaiš začetne kategorije
                global ped, categories, sex, active
                ped, categories, sex, active = nastavi_cat('GenPed_EBV.txt')
                #prestavi se v AlphaSim Dir
                os.chdir(AlphaSimDir)
                #kopiraj pedigre v selection folder
                shutil.copy('ExternalPedigree.txt', AlphaSimDir + '/Selection/SelectionFolder' + str(roundNo) + '/')
                #TUKAJ POTEM POPRAVIŠ AlphaSimSpec
                #PRVIČ PO BURN IN-U
                shutil.copy('/home/jana/Genotipi/Genotipi_CODES/AlphaSimSpec.txt', AlphaSimDir)
                os.system('sed -i "s|PedigreeType|ExternalPedigree.txt|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|EnterBurnInGenerationNumber|' + str(StBurnInGen) + '|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|EnterSelectionGenerationNumber|' + str(StSelGen) + '|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|EnterNumberOfSires|0|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|EnterNumberOfDams|0|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|TurnOnGenFlex|On|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|StartFlexGen,StopFlexGen|' +str(StBurnInGen + roundNo)+ ','+ str(StBurnInGen + roundNo) + '|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|TurnOnSelFlex|On|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|TheImportedGenerationPed|' +  str(StBurnInGen + roundNo) + '|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|TBVComputation|2|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|EnterIndividualInPopulation|' +str(stNB)+ '|g" AlphaSimSpec.txt') 

                #POŽENEŠ ALPHASIM        
                os.system('./AlphaSim1.05')
                    
    
            elif roundNo > 1:
                os.chdir('/home/jana/Genotipi/Genotipi_CODES/')
                #PRERAČUNAŠ EBV v Ru in ZAPIŠEŠ PEDIGRE
                shutil.copy ("Rcorr_PedEBV.R", "Rcorr_PedEBV_ThisGen.R")
                os.system('sed -i "s|AlphaSimPed|' + AlphaSimPed + '|g" Rcorr_PedEBV_ThisGen.R')
                call('Rscript Rcorr_PedEBV_ThisGen.R', shell=True)
                #tukaj nastvaiš začetne kategorije
                global ped, categories, sex, active
                ped, categories, sex, active = selekcija_ena_gen('GenPed_EBV.txt')
                #prestavi se v AlphaSim Dir
                os.chdir(AlphaSimDir)
                #kopiraj pedigre v selection folder
                shutil.copy('ExternalPedigree.txt', AlphaSimDir + '/Selection/SelectionFolder' + str(roundNo) + '/')
                #TUKAJ POTEM POPRAVIŠ AlphaSimSpec
                #PRVIČ PO BURN IN-U
                shutil.copy('/home/jana/Genotipi/Genotipi_CODES/AlphaSimSpec.txt', AlphaSimDir)
                os.system('sed -i "s|PedigreeType|ExternalPedigree.txt|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|EnterBurnInGenerationNumber|' + str(StBurnInGen) + '|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|EnterSelectionGenerationNumber|' + str(StSelGen) + '|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|EnterNumberOfSires|0|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|EnterNumberOfDams|0|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|TurnOnGenFlex|On|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|StartFlexGen,StopFlexGen|' +str(StBurnInGen + roundNo)+ ','+ str(StBurnInGen + roundNo) + '|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|TurnOnSelFlex|On|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|TheImportedGenerationPed|' +  str(StBurnInGen + roundNo) + '|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|TBVComputation|2|g" AlphaSimSpec.txt') 
                os.system('sed -i "s|EnterIndividualInPopulation|' +str(stNB)+ '|g" AlphaSimSpec.txt') 

                #POŽENEŠ ALPHASIM        
                os.system('./AlphaSim1.05')
                
                
        if BurnInYN == 'Y':  
            for roundNo in range(6, (StSelGen+1)):
                if roundNo == 1:
                    os.chdir('/home/jana/Genotipi/Genotipi_CODES/')
                    #PRERAČUNAŠ EBV v Ru in ZAPIŠEŠ PEDIGRE
                    shutil.copy ("Rcorr_PedEBV.R", "Rcorr_PedEBV_ThisGen.R")
                    os.system('sed -i "s|AlphaSimPed|' + AlphaSimPed + '|g" Rcorr_PedEBV_ThisGen.R')
                    call('Rscript Rcorr_PedEBV_ThisGen.R', shell=True)
                    #tukaj nastvaiš začetne kategorije
                    global ped, categories, sex, active
                    ped, categories, sex, active = nastavi_cat('GenPed_EBV.txt')

                    #prestavi se v AlphaSim Dir
                    os.chdir(AlphaSimDir)
                    #kopiraj pedigre v selection folder
                    shutil.copy('ExternalPedigree.txt', AlphaSimDir + '/Selection/SelectionFolder' + str(roundNo) + '/')
                    #TUKAJ POTEM POPRAVIŠ AlphaSimSpec
                    #PRVIČ PO BURN IN-U
                    shutil.copy('/home/jana/Genotipi/Genotipi_CODES/AlphaSimSpec.txt', AlphaSimDir)
                    os.system('sed -i "s|PedigreeType|ExternalPedigree.txt|g" AlphaSimSpec.txt') 
                    os.system('sed -i "s|EnterBurnInGenerationNumber|' + str(StBurnInGen) + '|g" AlphaSimSpec.txt') 
                    os.system('sed -i "s|EnterSelectionGenerationNumber|' + str(StSelGen) + '|g" AlphaSimSpec.txt') 
                    os.system('sed -i "s|EnterNumberOfSires|0|g" AlphaSimSpec.txt') 
                    os.system('sed -i "s|EnterNumberOfDams|0|g" AlphaSimSpec.txt') 
                    os.system('sed -i "s|TurnOnGenFlex|On|g" AlphaSimSpec.txt') 
                    os.system('sed -i "s|StartFlexGen,StopFlexGen|' +str(StBurnInGen + roundNo)+ ','+ str(StBurnInGen + roundNo) + '|g" AlphaSimSpec.txt') 
                    os.system('sed -i "s|TurnOnSelFlex|On|g" AlphaSimSpec.txt') 
                    os.system('sed -i "s|TheImportedGenerationPed|' +  str(StBurnInGen + roundNo) + '|g" AlphaSimSpec.txt') 
                    os.system('sed -i "s|TBVComputation|2|g" AlphaSimSpec.txt') 
                    os.system('sed -i "s|EnterIndividualInPopulation|' +str(stNB)+ '|g" AlphaSimSpec.txt') 
    
                
    
                    #POŽENEŠ ALPHASIM        
                    os.system('./AlphaSim1.05')
                        
        
                elif roundNo > 1:
                    os.chdir('/home/jana/Genotipi/Genotipi_CODES/')
                    #PRERAČUNAŠ EBV v Ru in ZAPIŠEŠ PEDIGRE
                    shutil.copy ("Rcorr_PedEBV.R", "Rcorr_PedEBV_ThisGen.R")
                    os.system('sed -i "s|AlphaSimPed|' + AlphaSimPed + '|g" Rcorr_PedEBV_ThisGen.R')
                    call('Rscript Rcorr_PedEBV_ThisGen.R', shell=True)
                    #global ped, categories, sex, active
                    ped, categories, sex, active = selekcija_ena_gen('GenPed_EBV.txt', categories=categories, sex=sex, active=active)
                    #prestavi se v AlphaSim Dir
                    os.chdir(AlphaSimDir)
                    #kopiraj pedigre v selection folder
                    os.system('mkdir ' +  AlphaSimDir + '/Selection/SelectionFolder' + str(roundNo))
                    shutil.copy('ExternalPedigree.txt', AlphaSimDir + '/Selection/SelectionFolder' + str(roundNo) + '/')
                    #TUKAJ POTEM POPRAVIŠ AlphaSimSpec
                    #PRVIČ PO BURN IN-U
                    shutil.copy('/home/jana/Genotipi/Genotipi_CODES/AlphaSimSpec.txt', AlphaSimDir)
                    os.system('sed -i "s|PedigreeType|ExternalPedigree.txt|g" AlphaSimSpec.txt') 
                    os.system('sed -i "s|EnterBurnInGenerationNumber|' + str(StBurnInGen) + '|g" AlphaSimSpec.txt') 
                    os.system('sed -i "s|EnterSelectionGenerationNumber|' + str(StSelGen) + '|g" AlphaSimSpec.txt') 
                    os.system('sed -i "s|EnterNumberOfSires|0|g" AlphaSimSpec.txt') 
                    os.system('sed -i "s|EnterNumberOfDams|0|g" AlphaSimSpec.txt') 
                    os.system('sed -i "s|TurnOnGenFlex|On|g" AlphaSimSpec.txt') 
                    os.system('sed -i "s|StartFlexGen,StopFlexGen|' +str(StBurnInGen + roundNo)+ ','+ str(StBurnInGen + roundNo) + '|g" AlphaSimSpec.txt') 
                    os.system('sed -i "s|TurnOnSelFlex|On|g" AlphaSimSpec.txt') 
                    os.system('sed -i "s|TheImportedGenerationPed|' +  str(StBurnInGen + roundNo) + '|g" AlphaSimSpec.txt') 
                    os.system('sed -i "s|TBVComputation|2|g" AlphaSimSpec.txt') 
                    os.system('sed -i "s|EnterIndividualInPopulation|' +str(stNB)+ '|g" AlphaSimSpec.txt') 
    
                    #POŽENEŠ ALPHASIM        
                    os.system('./AlphaSim1.05')





###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
#plot the results
#class TBVGenTable (SelectionTbvTest.txt)
from scipy import stats

TBVmeans.clear()
TBVmeans = defaultdict(list)
for roundNo in range(1,rounds+1):
    
    TBVt = TBVGenTable(AlphaSimDir + '/Selection/SelectionFolder' + str(roundNo) + '/SelectionTbvTest.txt')
    TBVmeans[roundNo] = TBVt.TBVmean
    
TBV  = TBVPed()
#plt.errorbar(x = TBV.gens, y = TBV.means, yerr = TBV.vars)

plt.plot( TBV.gens, TBV.means,  label = 'Mean Gen TBV')
plt.xlabel('Selected Generation')
plt.ylabel('Mean Generation TBV')
pylab.legend(loc='upper left')
plt.show()
plt.plot(TBV.gens, TBV.vars, label = 'TBV Var')
pylab.legend(loc='upper left')
plt.xlabel('Selected Generation')
plt.ylabel('Generation TBV variance')
plt.show()

