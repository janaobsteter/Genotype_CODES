# -*- coding: utf-8 -*-
import ast
import pandas as pd
from collections import defaultdict
import os
from selection10 import nastavi_cat_TGV
from selection10 import selekcija_total_TGV
from selection10 import *

WorkingDir = "/home/jana/"
os.chdir(WorkingDir)
scenario = "Class"

parhome = pd.read_csv(os.getcwd() + "/Essentials/SelectionParam_Class.csv", header=None, names=["Keys", "Vals"])
parhome.to_dict()
selParhome = defaultdict()
for key, val in zip(parhome.Keys, parhome.Vals):
    if key not in ['BurnInYN', 'EBV', 'gEBV', 'PA', 'AlphaSimDir', 'genotyped', 'EliteDamsPTBulls',
                   'EliteDamsPABulls', 'UpdateGenRef', 'sexToUpdate', 'EliteDamsGenBulls', 'gpb_pb',
                   'genTest_mladi', 'genTest_gpb', 'genFemale']:
        try:
            selParhome[key] = int(val)
        except:
            selParhome[key] = float(val)
    if key in ['BurnInYN', 'EBV', 'gEBV', 'PA', 'AlphaSimDir', 'EliteDamsPTBulls',
               'EliteDamsPABulls', 'UpdateGenRef', 'sexToUpdate', 'EliteDamsGenBulls', 'gpb_pb',
               'genTest_mladi', 'genTest_gpb', 'genFemale']:
        if val in ['False', 'True']:
            selParhome[key] = bool(val == 'True')
        else:
            selParhome[key] = val
    if key == 'genotyped':
        selParhome[key] = ast.literal_eval(val)

if selParhome['EBV']:
    seltype = 'class'
if selParhome['gEBV']:
    seltype = 'gen'

#     # tukaj pa še parametri za "large" population
# parimport = pd.read_csv(os.getcwd() + "/Essentials/SelectionParam_Class_LargePop.csv", header=None, names=["Keys", "Vals"])
# parimport.to_dict()
# selParimport = defaultdict()
# for key, val in zip(parimport.Keys, parimport.Vals):
#     if key not in ['BurnInYN', 'EBV', 'gEBV', 'PA', 'AlphaSimDir', 'genotyped', 'EliteDamsPTBulls',
#                    'EliteDamsPABulls', 'UpdateGenRef', 'sexToUpdate', 'EliteDamsGenBulls', 'gpb_pb',
#                    'genTest_mladi', 'genTest_gpb', 'genFemale']:
#         try:
#             selParimport[key] = int(val)
#         except:
#             selParimport[key] = float(val)
#     if key in ['BurnInYN', 'EBV', 'gEBV', 'PA', 'AlphaSimDir', 'EliteDamsPTBulls',
#                'EliteDamsPABulls', 'UpdateGenRef', 'sexToUpdate', 'EliteDamsGenBulls', 'gpb_pb',
#                'genTest_mladi', 'genTest_gpb', 'genFemale']:
#         if val in ['False', 'True']:
#             selParimport[key] = bool(val == 'True')
#         else:
#             selParimport[key] = val
#     if key == 'genotyped':
#         selParimport[key] = ast.literal_eval(val)

# BurnInYN = "False"  # ali izvedeš tudi BurnIn
# SelYN = "True"  # ali izvedeš tudi BurnIn
# StNB = 17280
# StBurnInGen = 20
# StSelGen = 40
# StartSelGen = 21
# StopSelGen = 40
# NumberOfSires = 30
# NumberOfDams = 8640
# AlphaSimDir = os.getcwd() + '/'
# selParimport['AlphaSimDir'] = os.getcwd()
# selParhome['AlphaSimDir'] = os.getcwd()
# if selParimport['EBV']:
#     seltype = 'class'
# if selParimport['gEBV']:
#     seltype = 'gen'

BurnInYN = "False"  # ali izvedeš tudi BurnIn
SelYN = "True"  # ali izvedeš tudi BurnIn
StNB = 8640
StBurnInGen = 20
StSelGen = 40
StartSelGen = 21
StopSelGen = 40
NumberOfSires = 6
NumberOfDams = 4320
AlphaSimDir = os.getcwd() + '/'
#selParimport['AlphaSimDir'] = os.getcwd()
selParhome['AlphaSimDir'] = os.getcwd()
# if selParimport['EBV']:
#     seltype = 'class'
# if selParimport['gEBV']:
#     seltype = 'gen'
if selParhome['EBV']:
    seltype = 'class'
if selParhome['gEBV']:
    seltype = 'gen'




#selPar["stNBn"] = 8640*2
# os.chdir(AlphaSimDir)
#
# selPar['pbn'] = 1
#
# selPar['genpbn'] = 1
# selPar['pozitivnoTestDoz'] = 2000
#
# selPar['pbUp'] = 1

IndCat = pd.DataFrame()
ped0, c0, s0, a0 = nastavi_cat("/home/jana/GenPed_EBV.txt", **selParhome)
createHerds = Herds(AlphaSimDir)  # to ne naredi nič, samo prebere datotek
createHerds.create_herds()  # ustvari črede, zapiši fajle
PedCat = OrigPed(AlphaSimDir, WorkingDir + '/CodeDir')
PedCat.addInfo()  # to ti zapiše PedigreeAndGeneticValues_cat.txt v AlphaSim/SimualatedData
# ped0, c0, s0, a0 = selekcija_total("/home/jana//GenPed_EBV.txt", **selPar)
blupNextGen = estimateBV(AlphaSimDir, WorkingDir + "/CodeDir", way='milk', sel=seltype)
blupNextGen.computeEBV_permEnv_herd(setVar=True, varPE=varPE, varE=varE, varH=varH,
                                    repeats=repeats)

inds = ped0.catCurrent_indiv('')
Cats = []
for ind in inds:
    for (key, value) in categories.iteritems():
        if ind in value:
            print(key)
            



Ind_['Indiv'] = ped.ped.Indiv
IndCat['catBurnIN'] = ped.ped.cat

krogov = 10

for krog in range(krogov): #ponavljaj kolikor krogov selekcije hočeš
    ped, c, s, a = selekcija_total("/home/jana/PedTotal.txt", **selParhome)
    IndCat[str('cat' + str(krog))] = ped.ped.cat





#
# set_group_two(AlphaSimDir, "home", "import", int(selParimport["stNBn"]))
# splitGenPed("PopulationSplit.txt")
# ped0, c0, s0, a0 = nastavi_cat('GenPed_EBVhome.txt', externalPedName = "ExternalPedigreehome", group=True, groupNumber=0, noGroups=2, **selParhome) #nastavi cat za prvo skupino
# ped1, c1, s1, a1 = nastavi_cat_TGV('GenPed_EBVimport.txt', externalPedName = "ExternalPedigreeimport", group=True, groupNumber=1, noGroups=2,  **selParimport) #nastavi cat za drugo skupino
# ped0.write_pedTotal("/home/jana/PedTotalhome.txt")
# ped1.write_pedTotal("/home/jana/PedTotalimport.txt")
# os.system("Rscript /home/jana/Genotipi/Genotipi_CODES/Combine_PedTotals.R")
# joinExternalPeds(["ExternalPedigreehome", "ExternalPedigreeimport"], AlphaSimDir)
# record_groups(["home", "import"], "PopulationSplit.txt")
#
# krogov = 2
# for krog in range(krogov): #ponavljaj kolikor krogov selekcije hočeš
#     splitGenPed("PopulationSplit.txt")
#     pedH, cH, sH, aH = selekcija_total('GenPed_EBVhome.txt', externalPedName = "ExternalPedigreehome", group=True, groupNumber=0, noGroups=2,  **selParhome)
#     pedI, cI, sI, aI = selekcija_total_TGV('GenPed_EBVimport.txt', externalPedName = "ExternalPedigreeimport", group=True, groupNumber=1, noGroups=2, **selParimport)
#     pedH.write_pedTotal("/home/jana/PedTotalhome.txt")
#     pedI.write_pedTotal("/home/jana/PedTotalimport.txt")
#     os.system("Rscript /home/jana/Genotipi/Genotipi_CODES/Combine_PedTotals.R")
#     joinExternalPeds(["ExternalPedigreehome", "ExternalPedigreeimport"], AlphaSimDir)
#     record_groups(["home", "import"], "PopulationSplit.txt")
#
#
# #test import
#     splitGenPed("PopulationSplit.txt")
#     pedI, cI, sI, aI = selekcija_total_TGV('GenPed_EBVimport.txt', externalPedName="ExternalPedigreeimport", group=True, groupNumber=1, noGroups=2, **selParimport)
#     if selParimport['EBV']:
#         Oce_import = pedI.izberi_ocete_PT(selParimport["pbUp"]) #tukaj so PT testirani očetje
#     if selParimport['gEBV']:
#         Oce_import = pedI.izberi_ocete_gen(selParimport["pbUp"])  # tukaj so genomsko testirani očetje
#
#     #odberi starše domače populacije
#     pedH, cH, sH, aH = selekcija_importOcetov('GenPed_EBVhome.txt', externalPedName="ExternalPedigreehome", group=True, groupNumber=0, noGroups=2,
#                                            importBool=True, importGroup="bm", FatherList=Oce_import, **selParhome)
#     pedH.write_pedTotal("/home/jana/PedTotalhome.txt")
#     pedI.write_pedTotal("/home/jana/PedTotalimport.txt")
#     os.system("Rscript /home/jana/Genotipi/Genotipi_CODES/Combine_PedTotals.R")
#
#     joinExternalPeds(["ExternalPedigreehome", "ExternalPedigreeimport"], AlphaSimDir)
#     record_groups(["home", "import"], "PopulationSplit.txt")