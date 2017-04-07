# -*- coding: utf-8 -*-
from __future__ import division
from collections import defaultdict
import pandas as pd
import numpy as np
from collections import defaultdict
import random
reload(selection)
from selection import *
# -*- coding: utf-8 -*-

#############################################33
#tukaj naštimaj parametre
#1) odstotki živali v kateogorijah
nT = 20000 # velikost populacije
perF = 0.9
perM = 0.1 # odstotki ženskih in moških živali v aktivni populaciji
#odstotki kategorij znotraj spola
nrF = 0.15
pt12 = 0.11
t12 = 0.035
pt24 = 0.05
t24 = 0.03
k = 0.625
sum([nrF, pt12, t12, pt24, t24, k])


nrM = 0.55
mladi = 0.0073 #TO JE ZNOTRAJ NOVOROJENIH BIKOV!!!! (drugače procent 0.4%)
bik12 = 0.40
pripust = 0.041 #TO JE ZNOTRAJ BIKOV 12!!!!
#cakajoci = 0.018 #to je znotraj bikov12!!!!!
bik24 = 0.05
potomciNP = 0.041 #%OD NOVOROJENIH
vhlevljeni = 0.6 #TO JE ODSTOTEK OD POTOMCEV NAČRNIH PARJENJ


pb = 0.009
nrM + bik12 +  bik24 

#2)števila let v uporabi
kraveUp = 5 #povprečno koliko let so krave v populaciji (koliko laktacij)
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
nF = nT * perF
nrFn = int(nrF * nF)
pt12n = int(pt12 * nF)
t12n = int(t12 * nF)
pt24n = int(pt24* nF)
t24n = int(t24 * nF)
kn = int(k * nF)
sum([nrFn, pt12n, t12n, pt24n, t24n, kn])

nM = nT * perM
nrMn = int(nrM * nM)
bik12n = int(bik12 * nM)
bik24n = int(bik24 * nM)
potomciNPn = int(potomciNP * nM * nrM)
mladin = int(mladi * nrMn) ###Mladi: mladin = int(nrM * nM * mladi)
vhlevljenin  = int(potomciNP * nM * nrM * vhlevljeni)
#cakajocin = int(cakajoci* nM)
pripustn = int(pripust * bik12 * nM)
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
from selection import *
ped = pedigree("~/Documents/PhD/Simulaton/Pedigrees/PedPython.txt")

#VEDNO NAJPREJ IZLOČI /ODBERI PO PV!!! - funckije za odbiro na random imajo pogoj, da je kateogrija prosta
def select_age_0_1(ped):
    #FEMALES
    ped.set_sex(0, nrFn, "F") #choose female new borns from generation before
    izlF = nrFn - pt12n - t12n #koliko jih izločiš
    ped.izloci_poEBV(1, "F", izlF) #tukaj jih izloči, funkcija v modulu
    ped.izberi_random(1, "F", t12n, "t12") #na random izberi telice, ki jih preneseš naprej
    ped.izberi_random(1, "F", pt12n, "pt12") #na random izberi osemenjene telice

    #MALES
    ped.set_sex(nrFn, (nrFn + nrMn), "M") #nastavi sex = M od zadnje novorojene teličke do
    izlM = nrMn - bik12n #choose female new borns from generation before
    ped.izberi_poEBV_top(1, "M", mladin, "mladi1") #odberi mlade
    ped.izberi_poEBV_OdDo(1, "M", mladin, (mladin + (vhlevljenin - mladin)), "pripust1") #odberi v pripustu
    ped.izberi_random(1, "M", bik12n, "bik12")
    ped.izloci_random(1, "M", (nrMn - bik12n - vhlevljenin))
    categories = ped.save_cat()
    sex = ped.save_sex()
    active = ped.save_active()
    