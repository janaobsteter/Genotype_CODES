# -*- coding: utf-8 -*-
from __future__ import division
from collections import defaultdict
import pandas as pd
import numpy as np
from collections import defaultdict
import random
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
t12 = 0.966
pt = 0.85 #koliko jih pride do telic
kraveRemont = 0.25
bm = 0.0127 #TO JE OD KRAV!!!
sum([nrF, pt12, t12, pt24, t24, k])


nrM = 0.55
mladi = 0.0073 #TO JE ZNOTRAJ NOVOROJENIH BIKOV!!!! (drugače procent 0.4%)
bik12 = 0.40
pripust = 0.0135 #TO JE ZNOTRAJ BIKOV 12!!!!
#cakajoci = 0.018 #to je znotraj bikov12!!!!!
bik24 = 0.05
potomciNP = 0.041 #%OD NOVOROJENIH
vhlevljeni = 0.6 #TO JE ODSTOTEK OD POTOMCEV NAČRNIH PARJENJ


pb = 0.009
nrM + bik12 +  bik24 

#2)števila let v uporabi
kraveUp = 4 #povprečno koliko let so krave v populaciji (koliko laktacij)
bmUp = 3 # koliko let so v uporabi BM - REMONT!
cak = 4 #koliko časa so mladi biki v testu oz. koliko časa so čakajoči
pbUp = 5 #koliko let so povrpečno v uporabi biki v AI
pripustUp = 1.4 # koliko let so v uporabi biki v pripustu
genomUp = 1.3 #koliko let so povprečno v uporabi genomsko testirani biki

##################################################################################################################3
##################################################################################################################3
#od tu naprej samo delo s parametri = spremenvljivkami
##################################################################################################################3

#številke
#ženske
nrFn = int(stNB * 0.5)

t12n = int(t12 * nrFn)
ptn = int(pt* t12n)
t24n = int(t24 * nF)
kn = int(k * nF)
bmn = int(bm * k * nF)
sum([nrFn, pt12n, t12n, pt24n, t24n, kn])

nM = 0.5
nrMn = int(stNB * nM)
bik12n = int(bik12 * nM)
bik24n = int(bik24 * nM)
potomciNPn = int(potomciNP * nM * nrM)
mladin = int(mladi * nrMn) ###Mladi: mladin = int(nrM * nM * mladi)
vhlevljenin  = int(potomciNP * nM * nrM * vhlevljeni)
#cakajocin = int(cakajoci* nM)
pripustTn = int(pripust *  nM)
pripust1n = vhlevljenin - mladin
pripust2n = pripustTn - pripust1n
pbn = int(pb * nM)
#sum([nrMn, bik12n, bik24n, cakajocin, pripustn, pbn, mladin])
sum([nrMn, bik12n, bik24n])
##########################################################
#ženske
#odberi novorojene
#tu se ne rabiš obremenjevat s tem, katere so matere in kateri očetje - to narediš na koncu, ko konstruiraš novorojene za naslednjo gene
#ti novorojeni so že prva generacija
#na koncu potegni skupaj ID in kategorije ter to uporabi za določitev staršev#
#tudi šele tam se ukvarjaj z [pb]
reload(selection)
import selection
from selection import *
ped = pedigree("~/Documents/PhD/Simulaton/Pedigrees/PedPython.txt")



ped.set_cat_gen(max(ped.gen), "nr") #to je samo na prvem loopu

categories = ped.save_cat()

ped = pedigree("~/Documents/PhD/Simulaton/Pedigrees/PedPython.txt")

ped.set_sex(0, nrFn, "F") #choose female new borns from generation before

#prva odbira
select_age_0_1(ped)
ped.add_new_gen_blank(nrFn, nrMn)
     
categories.clear()
categories = ped.save_cat()
sex = ped.save_sex() 


#druga odbira
ped.set_cat_gen(1, "")
ped.set_cat_gen(2, "")
ped.set_cat_old('izl', 'izl', categories)

select_age_0_1(ped) 
select_age_1_2(ped)
 
ped.add_new_gen_blank(nrFn, nrMn)  

categories.clear()
categories = ped.save_cat()
sex = ped.save_sex()
active = ped.save_active()
age = ped.save_age() 
    
    

#################################

 #tretja odbira
ped.set_cat_gen(1, "")
ped.set_cat_gen(2, "")
ped.set_cat_gen(3, "")  
ped.set_cat_old('izl', 'izl', categories)


#TUKAJ NE MOREŠ ODBIRATI PO EBV DOKLER DA NIMAŠ ZGENERIRANIH!
select_age_0_1(ped)
select_age_1_2(ped)
ped.compute_age()
select_age_2_3(ped)

ped.add_new_gen_blank(nrFn, nrMn)  

categories.clear()
categories = ped.save_cat()
sex = ped.save_sex()
active = ped.save_active()
age = ped.save_age() 
    
#četrti krog
ped.set_cat_gen(1, "")
ped.set_cat_gen(2, "")
ped.set_cat_gen(3, "") 
ped.set_cat_gen(4, "")   
ped.set_cat_old('izl', 'izl', categories) 

select_age_0_1(ped)
select_age_1_2(ped)
ped.compute_age()
select_age_2_3(ped)

ped.add_new_gen_blank(nrFn, nrMn)  

categories.clear()
categories = ped.save_cat()
sex = ped.save_sex()
active = ped.save_active()
age = ped.save_age() 


#PETI KROG
ped.set_cat_gen(1, "")
ped.set_cat_gen(2, "")
ped.set_cat_gen(3, "") 
ped.set_cat_gen(4, "")   
ped.set_cat_gen(5, "")
ped.set_cat_old('izl', 'izl', categories) 

select_age_0_1(ped)
select_age_1_2(ped)
ped.compute_age()
select_age_2_3(ped)

ped.add_new_gen_blank(nrFn, nrMn)  

categories.clear()
categories = ped.save_cat()
sex = ped.save_sex()
active = ped.save_active()
age = ped.save_age() 

#ŠESTI KROG
ped.set_cat_gen(1, "")
ped.set_cat_gen(2, "")
ped.set_cat_gen(3, "") 
ped.set_cat_gen(4, "")   
ped.set_cat_gen(5, "")
ped.set_cat_gen(6, "")
ped.set_cat_old('izl', 'izl', categories) 

select_age_0_1(ped)
select_age_1_2(ped)
ped.compute_age()
select_age_2_3(ped)

ped.add_new_gen_blank(nrFn, nrMn)  

categories.clear()
categories = ped.save_cat()
sex = ped.save_sex()
active = ped.save_active()
age = ped.save_age() 

#SEDMI KROG
ped.set_cat_gen(1, "")
ped.set_cat_gen(2, "")
ped.set_cat_gen(3, "") 
ped.set_cat_gen(4, "")   
ped.set_cat_gen(5, "")
ped.set_cat_gen(6, "")
ped.set_cat_gen(7, "")
ped.set_cat_old('izl', 'izl', categories) 

select_age_0_1(ped)
select_age_1_2(ped)
ped.compute_age()
select_age_2_3(ped)

ped.add_new_gen_blank(nrFn, nrMn)  

categories.clear()
categories = ped.save_cat()
sex = ped.save_sex()
active = ped.save_active()
age = ped.save_age() 

#OSMI KROG
ped.set_cat_gen(1, "")
ped.set_cat_gen(2, "")
ped.set_cat_gen(3, "") 
ped.set_cat_gen(4, "")   
ped.set_cat_gen(5, "")
ped.set_cat_gen(6, "")
ped.set_cat_gen(7, "")
ped.set_cat_gen(8, "")
ped.set_cat_old('izl', 'izl', categories) 

select_age_0_1(ped)
select_age_1_2(ped)
ped.compute_age()
select_age_2_3(ped)

ped.add_new_gen_blank(nrFn, nrMn)  

categories.clear()
categories = ped.save_cat()
sex = ped.save_sex()
active = ped.save_active()
age = ped.save_age() 

#DEVETI KROG
ped.set_cat_gen(1, "")
ped.set_cat_gen(2, "")
ped.set_cat_gen(3, "") 
ped.set_cat_gen(4, "")   
ped.set_cat_gen(5, "")
ped.set_cat_gen(6, "")
ped.set_cat_gen(7, "")
ped.set_cat_gen(8, "")
ped.set_cat_gen(9, "")
ped.set_cat_old('izl', 'izl', categories) 

select_age_0_1(ped)
select_age_1_2(ped)
ped.compute_age()
select_age_2_3(ped)

ped.add_new_gen_blank(nrFn, nrMn)  

categories.clear()
categories = ped.save_cat()
sex = ped.save_sex()
active = ped.save_active()
age = ped.save_age() 

#DESETI KROG
ped.set_cat_gen(1, "")
ped.set_cat_gen(2, "")
ped.set_cat_gen(3, "") 
ped.set_cat_gen(4, "")   
ped.set_cat_gen(5, "")
ped.set_cat_gen(6, "")
ped.set_cat_gen(7, "")
ped.set_cat_gen(8, "")
ped.set_cat_gen(9, "")
ped.set_cat_gen(10, "")
ped.set_cat_old('izl', 'izl', categories) 

select_age_0_1(ped)
select_age_1_2(ped)
ped.compute_age()
select_age_2_3(ped)

ped.add_new_gen_blank(nrFn, nrMn)  

categories.clear()
categories = ped.save_cat()
sex = ped.save_sex()
active = ped.save_active()
age = ped.save_age() 
##############################################################################################3
##############################################################################################3
##############################################################################################3    
    
#VEDNO NAJPREJ IZLOČI /ODBERI PO PV!!! - funckije za odbiro na random imajo pogoj, da je kateogrija prosta
def select_age_0_1(ped): #tukaj odbereš iz novorojenih živali t12, pt12 in mlade bike, pripust1
    #FEMALES
    izlF = nrFn - t12n#koliko jih izločiš
    ped.izloci_poEBV("F", izlF, "nr", categories) #tukaj jih izloči, funkcija v modulu
    ped.izberi_poEBV_top("F", (nrFn - izlF), "nr", "t12", categories) #izberi telice, ki jih osemeniš --> krave
    
    
    #MALES
    ped.set_sex(nrFn, (nrFn + nrMn), "M") #nastavi sex = M od zadnje novorojene teličke do
    izlM = nrMn - bik12n #choose female new borns from generation before
    ped.izberi_poEBV_top( "M", mladin, "nr", "mladi", categories) #odberi mlade
    ped.izberi_poEBV_OdDo(max(ped.gens()), "M", mladin, (mladin + (vhlevljenin - mladin)), "nr", "pripust1", categories) #odberi v pripustu
    ped.izberi_random(max(ped.gens()), "M", bik12n, "nr", "bik12", categories)
    ped.izloci_random(max(ped.gens()), "M", (nrMn - bik12n - vhlevljenin),"nr", categories)
    

def select_age_1_2(ped): # tukaj odbereš nič pri kravah - razen, če so že bikovske matere, pripust 2, bike24

    currentGen = max(ped.gens())-1
    #FEMALES
    ped.izberi_poEBV_top("F", ptn, 't12', 'pt', categories)
    ped.izloci_poEBV("F", (t12n - ptn),'t12', categories) #terlice postanejo
   
    
    #MALES
    ped.izberi_random(currentGen, "M", bik24n, 'bik12', 'bik24', categories)
    ped.izberi_random(currentGen, "M", (bik12n - bik24n), 'bik12', 'izl', categories)
    ped.set_cat_old('mladi', 'cak', categories) #mladih ne prestavljaj, ampak samo daj v pb, ko dosežejo starost
    ped.izberi_random(currentGen, "M", pripust2n, 'pripust1', 'pripust2', categories)
    ped.izloci_random(currentGen, "M", (pripust1n - pripust2n), 'pripust1', categories)


#tukaj lahk daš vse v eno funkcijo - variabilno - koliko let krave, koliko let v testu
def select_age_2_3(ped):
    #FEMALES
    #najprej dodaj nove krave
    ped.set_cat_old('pt', 'k', categories) #osemenjene telice postanejo krave - predpostavimo, da vse
    #potem izloči najstarejše krave - po 4. laktaciji
    if ('k' in categories.keys()) and ((kraveUp+1) in ped.age()):
        ped.izloci_age_cat((kraveUp+1), 'k', categories)
    #ostale krave prestavi naprej v krave
    ped.set_cat_age_old(3, 'k', 'k', categories)
    ped.set_cat_age_old(4, 'k', 'k', categories) 
    ped.set_cat_age_old(5, 'k', 'k', categories)    
              
    #če imaš že dovolj stare krave, potem odberi BM
    #BM se odbira po drugi laktaciji - to je starost 3 - 4 (starost v pedigreju = 3, ker imaš tudi 0)
    if ('k' in categories.keys()) and (3 in ped.age()):
        ped.izberi_poEBV_top_age("F",3, int(bmn /bmUp), 'k', 'bm', categories) #izberi bikovske matere
    #in izloči najastarejše BM, če jih imaš
    if ('bm' in categories.keys()) and (6 in ped.age()):
        ped.izloci_age_cat(6, 'bm', categories)
    #ostale BM prestavi naprej
    ped.set_cat_age_old(4, 'bm', 'bm', categories)
    ped.set_cat_age_old(5, 'bm', 'bm', categories)
    
    #MALES
    #mlade bike prestavi naprej (starost 1,2,3)
    ped.set_cat_age_old(2, 'cak', 'cak', categories)
    ped.set_cat_age_old(3, 'cak', 'cak', categories)
    ped.set_cat_age_old(4, 'cak', 'cak', categories)
    #plemenske bike prestavljaj naprej
    ped.set_cat_old('pb', 'pb', categories)
    ped.izloci_cat('bik24', categories)
    ped.izloci_cat('pripust2', categories)
    if ('cak' in categories.keys()) and ((cak+1) in ped.age()):
        ped.izberi_poEBV_top_age("M", (cak +1), int(mladin * 0.5), 'cak', 'pb', categories)
        ped.izloci_poEBV_age("M",(cak+1), int(mladin * 0.5), 'cak', categories) #TUKAJ MORA BITI ŠE STAROST!!!!!!!!!!!

 
 
 """   
def select_age_3_4(ped):
    #FEMALES
    ped.set_cat_old('k2', 'k3', categories)
    
    
    #MALES
    ped.set_cat_old('mladi3', 'mladi4', categories)
 
def select_age_4_5(ped):
    #FEMALES
    #ped.set_cat_old('bm1', 'bm2', categories)
    ped.set_cat_old('k3', 'k4', categories)
    
    #MALES
    ped.set_cat_old('mladi4', 'mladi5', categories)
 """  
def select_age_5_6(ped):
    ped.set_cat_old('k4', 'k5', categories)
    #ped.set_cat_old('bm2', 'bm3', categories)
    


def osemeni_bm_elita(ped):
      