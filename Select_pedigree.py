# -*- coding: utf-8 -*-
from __future__ import division
from collections import defaultdict
import pandas as pd
import numpy as np
from collections import defaultdict
import random
from itertools import chain
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

teln = int(tel * nrFn)
ptn = int(pt* teln)
bmn = int(round(ptn * kraveUp *bm)) #to je od vseh krav
#kn = int(k * nF)
#bmn = int(bm * k * nF)
#sum([nrFn, pteln, teln, pt24n, t24n, kn])

nM = 0.5
nrMn = int(round(stNB * 0.5 * (1 - potomciNP))) 
telMn = int(round(telM * nrMn))
bik24n = int(round(bik24 * telMn))
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
ped = pedigree("~/Documents/PhD/Simulaton/Pedigrees/PedPython.txt")    

###################################
#loop za vajo
#################################

 #tretja odbira
 #TUKAJ NE MOREŠ ODBIRATI PO EBV DOKLER DA NIMAŠ ZGENERIRANIH! - zato random za vajo
stevilo_krogov = 50

for krog in (range(0,stevilo_krogov)):
    if krog == 0:
        ped.set_cat_gen(max(ped.gen), "nr") #to je samo na prvem loopu
        ped.set_sex_list([x for x in range(0,ped.rows()) if x%2==0], "F")
        ped.set_sex_list([x for x in range(0,ped.rows()) if x%2!=0], "M")
        ped.izberi_poEBV_top_catCurrent("F", int(potomciNPn/2), 'nr', 'potomciNP')
        ped.izberi_poEBV_top_catCurrent("M", int(potomciNPn/2), 'nr', 'potomciNP')
        categories = ped.save_cat()
        ped = pedigree("~/Documents/PhD/Simulaton/Pedigrees/PedPython.txt")
        ped.set_sex_list([x for x in range(0,ped.rows()) if x%2==0], "F")
        ped.set_sex_list([x for x in range(0,ped.rows()) if x%2!=0], "M")
        #prva odbira
        ped.compute_age()        
        select_age_0_1(ped)
        ped.add_new_gen_naive(stNB, stpotomciNPn)
            
        categories.clear()
        categories = ped.save_cat()
        sex = ped.save_sex()
    
    if krog == 1:
        #SETSEX!!!
        #druga odbira
        ped.set_cat_gen(1, "")
        ped.set_cat_gen(2, "")
        ped.set_cat_old('izl', 'izl', categories)
        
        ped.compute_age()        
        select_age_0_1(ped) 
        select_age_1_2(ped)
        
        ped.add_new_gen_naive(stNB, stpotomciNPn)  
        
        categories.clear()
        categories = ped.save_cat()
        sex = ped.save_sex()
        active = ped.save_active()
        age = ped.save_age() 
    
    if krog >= 2:
        for i in ped.gens():
            ped.set_cat_gen(i, "")
            
        ped.set_cat_old('izl', 'izl', categories) 

        ped.compute_age()        
        select_age_0_1(ped)
        select_age_1_2(ped)

        select_age_2_3(ped)

        ped.add_new_gen_naive(stNB, stpotomciNPn) 
        ped.compute_age() 
        doloci_matere(ped)
        doloci_ocete(ped)
        #ped.set_cat_mother_catCurrent('bm', 'potomciNP') #TO DAJ V FUNKCIJO!
        
        categories.clear()
        categories = ped.save_cat()
        sex = ped.save_sex()
        active = ped.save_active()
        age = ped.save_age() 
        ped.compute_age() #drugače imaš negativne vrednosti!
        


     



##############################################################################################3
##############################################################################################3
##############################################################################################3    
    
#VEDNO NAJPREJ IZLOČI /ODBERI PO PV!!! - funckije za odbiro na random imajo pogoj, da je kateogrija prosta
def select_age_0_1(ped): #tukaj odbereš iz novorojenih živali tel, ptel in mlade bike, pripust1
    #FEMALES
    ped.set_cat_sex_old("F", "potomciNP", "tel", categories)
    izlF = nrFn - teln#koliko jih izločiš
    ped.izloci_poEBV("F", izlF, "nr", categories) #tukaj jih izloči, funkcija v modulu

    ped.izberi_poEBV_top("F", (nrFn - izlF), "nr", "tel", categories) #izberi telice, ki jih osemeniš --> krave
    
    
    #MALES
    ped.izberi_poEBV_top( "M", vhlevljenin, "potomciNP", "vhlevljeni", categories) #odberi mlade TO SAMO NA ZAČETKU; POTEM POTOMCI BM IN ELITE!
    ped.izloci_poEBV("M", int((len(categories['potomciNP'])/2) - vhlevljenin), 'potomciNP', categories)
    ped.izberi_random( "M", bik12n, "nr", "bik12", categories)
    ped.izloci_random( "M", int(nrMn - bik12n - (len(categories['potomciNP'])/2)),"nr", categories)
    

def select_age_1_2(ped): # tukaj odbereš nič pri kravah - razen, če so že bikovske matere, pripust 2, bike24
    #FEMALES
    ped.izberi_poEBV_top("F", ptn, 'tel', 'pt', categories)
    ped.izloci_poEBV("F", (len(categories['tel']) - ptn),'tel', categories) #terlice postanejo
   
    
    #MALES
    ped.izberi_poEBV_top( "M", mladin, "vhlevljeni", "mladi", categories) #odberi mlade
    ped.izberi_poEBV_OdDo( "M", mladin, (mladin + (vhlevljenin - mladin)), "vhlevljeni", "pripust1", categories) #odberi v pripustu
    ped.izberi_random( "M", bik24n, 'bik12', 'bik24', categories)
    ped.izloci_random( "M", (len(categories['bik12']) - bik24n), 'bik12', categories)



#tukaj lahk daš vse v eno funkcijo - variabilno - koliko let krave, koliko let v testu
def select_age_2_3(ped):
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
    if ('k' in categories.keys()) and (3 in ped.age()):
        ped.izberi_poEBV_top_age("F",3, int(bmn /bmUp), 'k', 'bm', categories) #izberi bikovske matere
    #in izloči najastarejše BM, če jih imaš
    if ('bm' in categories.keys()) and ((bmUp+3) in ped.age()):
        ped.izloci_age_cat((bmUp+3), 'bm', categories)
    #ostale BM prestavi naprej
    ped.set_cat_age_old(4, 'bm', 'bm', categories)
    ped.set_cat_age_old(5, 'bm', 'bm', categories)
    
    #MALES
    #mladi biki postanejo čakajoči (~1 leto, da se osemeni krave s semenom oz. malo po 2. letu)
    ped.set_cat_old('mladi', 'cak', categories) 
    ped.set_active_cat('mladi', 2, categories)
    
    #čakajočim bikov podaljšaj status (do starosti 5 let)
    #hkrati jim tudi nastavi status izl
    #ped.set_cat_age_old(2, 'cak', 'cak', categories)
    ped.set_cat_age_old(3, 'cak', 'cak', categories)
    ped.set_cat_age_old(4, 'cak', 'cak', categories)
    
    #povprečna doba v pripustu - glede na to odberi bike, ki preživijo še eno leto
    if 'pripust1' in categories.keys():
        ped.izberi_random( "M", pripust2n, 'pripust1', 'pripust2', categories)
        ped.izloci_random( "M", (pripust1n - pripust2n), 'pripust1', categories)

    #plemenske bike prestavljaj naprej
    ped.set_cat_old('pb', 'pb', categories)
    ped.izloci_cat('bik24', categories)
    ped.izloci_cat('pripust2', categories)
    if ('cak' in categories.keys()) and ((cak+1) in ped.age()):
        ped.izberi_poEBV_top_age("M", (cak +1), int(mladin * 0.5), 'cak', 'pb', categories)
        ped.set_active_cat('cak', 2, categories) #tukaj moraš to nastaviti, zato ker fja izberi avtomatsko nastavi na active=1
        ped.izloci_poEBV_age("M",(cak+1), int(mladin * 0.5), 'cak', categories) #TUKAJ MORA BITI ŠE STAROST!!!!!!!!!!!



def doloci_matere(ped):
    #MATERE
    sTbmMother = 90 if len(ped.catCurrent_indiv('bm')) >= 90 else len(ped.catCurrent_indiv('bm'))
    if sTbmMother != 0:
        bmMother = ped.select_mother_random('bm', sTbmMother)
        ped.set_mother_catPotomca(bmMother, 'potomciNP')
    #
    
    if 'k' in ped.cat():#TUKAJ SO DOLOČENE SEDAJ VSE MATERE!!!
        mother = ped.select_mother_EBV_top('k', int(round(11000*0.7))) #tukaj odberi brez tistih, ki so za gospodarsko križanje
        if len(mother) >= (stNB - sTbmMother): # če že imaš dovolj krav, določi matere vsem novorojenim oz. odbiraš matere, saj jih imaš preveč!
            motherOther = random.sample(mother, (stNB - len(bmMother)))
            ped.set_mother_catPotomca(motherOther, 'nr') #TUKAJ SO DOLOČENE SEDAJ VSE MATERE!!!
        elif len(mother) < (stNB - sTbmMother): # če jih še ni dovolj, ne odbiraš mater, ampak uporabiš vse - gosp. križanmja
            ped.set_mother_catPotomca(mother, 'nr') 


def doloci_ocete(ped):
#OČETJE
    mladiOce = ped.catCurrent_indiv('mladi')
    pripustOce = ped.catCurrent_indiv('pripust1') + ped.catCurrent_indiv('pripust2') 
    testiraniOce = list(chain.from_iterable([ped.catCurrent_indiv_age('pb', (1 + cak + x)) for x in range(1, pbUp+1)]))
    bmMother = 90 if len(ped.catCurrent_indiv('bm')) >= 90 else len(ped.catCurrent_indiv('bm'))
    elita = np.random.choice(range(-20,-1), bmMother, replace=True) #navidezna elita
    pd.Series(elita).value_counts()#preveri zastopanost po bikih
    
    ocetje = pripustOce*pripustDoz + testiraniOce*pozitivnoTestDoz + mladiOce*mladiDoz
    if len(ocetje) >= (stNB - potomciNPn*2): #če imaš dovolj DOZ za vse NB
        ocetjeNB = random.sample(ocetje, (stNB - potomciNPn*2)) #tukaj izbereš očete za vse krave  - razen BM!
        ped.set_father_catPotomca(ocetjeNB, 'nr')
    if len(ocetje) < (stNB - potomciNPn*2):
        ped.set_father_catPotomca(ocetje, 'nr')

    #naštimaj očete elite --> BM
    ped.set_father_catPotomca(elita, 'potomciNP')





#####################################################################
#tukaj je zdj funkcija, ki vse to dela!
##################################################################### 

#to je funkcija za odbiro in določanje staršev
#prvi pogoj if max gen = 1 je za primer, ko štartaš s praznim naivnim pedigrejem brez staršev - mam in očetov ni v pedigreju
#drugi pogoj,ko dodaš generacijo novorojenih in pelješ prejšnjo generacijo naprej
#tretji krog so združene vse selekcijske odločitve po tem - počasi dobiš bm in pb, če jih ni, se pač ti starši ne določajo
def selekcija_ena_gen(pedFile):
    ped = pedigree(pedFile) 
    
    if max(ped.gen) == 1:
        ped.set_cat_gen(max(ped.gen), "nr") #to je samo na prvem loopu
        ped.set_sex_list([x for x in range(0,ped.rows()) if x%2==0], "F")
        ped.set_sex_list([x for x in range(0,ped.rows()) if x%2!=0], "M")
        ped.izberi_poEBV_top_catCurrent("F", int(potomciNPn/2), 'nr', 'potomciNP')
        ped.izberi_poEBV_top_catCurrent("M", int(potomciNPn/2), 'nr', 'potomciNP')
        
        global categories #to moraš dat global samo v prvenm loopu, drugje dobiš return
        categories = ped.save_cat()
        ped = pedigree(pedFile)
        ped.set_sex_list([x for x in range(0,ped.rows()) if x%2==0], "F")
        ped.set_sex_list([x for x in range(0,ped.rows()) if x%2!=0], "M")
        #prva odbira
        ped.compute_age()        
        select_age_0_1(ped)
        ped.add_new_gen_naive(stNB, stpotomciNPn)
        
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
    
    
    if max(ped.gen) == 2:
        print "DRUGA SELEKCIJA"
        #druga odbira
        ped.set_cat_gen(1, "")
        ped.set_cat_gen(2, "")
        ped.set_cat_old('izl', 'izl', categories)
        
        ped.compute_age()
        select_age_0_1(ped) 
        select_age_1_2(ped)
        
        ped.add_new_gen_naive(stNB, stpotomciNPn)  
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
        print "TRETJA SELEKCIJA"
        for i in ped.gens():
            ped.set_cat_gen(i, "")
            
        ped.set_cat_old('izl', 'izl', categories) 
        
        ped.compute_age()        
        select_age_0_1(ped)
        select_age_1_2(ped)
        select_age_2_3(ped)
        
        ped.add_new_gen_naive(stNB)  
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
    
    return ped.save_cat(), ped.save_sex(), ped.save_active()



      