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


#parhome = pd.read_csv(os.getcwd() + "/Essentials/SelectionParam_ClassTest.csv", header=None, names=["Keys", "Vals"])
parhome = pd.read_csv(os.getcwd() + "/Essentials/SelectionParam_Class.csv", header=None, names=["Keys", "Vals"])
parhome.to_dict()
selParhome = defaultdict()
for key, val in zip(parhome.Keys, parhome.Vals):
    if key not in ['importPer', 'BurnInYN', 'EBV', 'gEBV', 'PA', 'AlphaSimDir', 'genotyped', 'EliteDamsPTBulls',
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
    if key in ['genotyped', 'importPer']:
        selParhome[key] = ast.literal_eval(val)


if selParhome['EBV']:
    seltype = 'class'
if selParhome['gEBV']:
    seltype = 'gen'

    # tukaj pa še parametri za "large" population
#parimport = pd.read_csv(os.getcwd() + "/Essentials/SelectionParam_Class_LargePop_Test.csv", header=None, names=["Keys", "Vals"])
parimport = pd.read_csv(os.getcwd() + "/Essentials/SelectionParam_Class_LargePop.csv", header=None, names=["Keys", "Vals"])
parimport.to_dict()
selParimport = defaultdict()
for key, val in zip(parimport.Keys, parimport.Vals):
    if key not in ['BurnInYN', 'EBV', 'gEBV', 'PA', 'AlphaSimDir', 'genotyped', 'EliteDamsPTBulls',
                   'EliteDamsPABulls', 'UpdateGenRef', 'sexToUpdate', 'EliteDamsGenBulls', 'gpb_pb',
                   'genTest_mladi', 'genTest_gpb', 'genFemale']:
        try:
            selParimport[key] = int(val)
        except:
            selParimport[key] = float(val)
    if key in ['BurnInYN', 'EBV', 'gEBV', 'PA', 'AlphaSimDir', 'EliteDamsPTBulls',
               'EliteDamsPABulls', 'UpdateGenRef', 'sexToUpdate', 'EliteDamsGenBulls', 'gpb_pb',
               'genTest_mladi', 'genTest_gpb', 'genFemale']:
        if val in ['False', 'True']:
            selParimport[key] = bool(val == 'True')
        else:
            selParimport[key] = val
    if key == 'genotyped':
        selParimport[key] = ast.literal_eval(val)

BurnInYN = "False"  # ali izvedeš tudi BurnIn
SelYN = "True"  # ali izvedeš tudi BurnIn
StNB = 17280
StBurnInGen = 20
StSelGen = 40
StartSelGen = 21
StopSelGen = 40
NumberOfSires = 20
NumberOfDams = 8640
AlphaSimDir = os.getcwd() + '/'
selParhome['AlphaSimDir'] = os.getcwd()
AlphaSimPed = selParhome['AlphaSimDir'] + '/SimulatedData/PedigreeAndGeneticValues.txt'
if selParhome['EBV']:
    seltype = 'class'
if selParhome['gEBV']:
    seltype = 'gen'

# BurnInYN = "False"  # ali izvedeš tudi BurnIn
# SelYN = "True"  # ali izvedeš tudi BurnIn
# StNB = 200
# StBurnInGen = 10
# StSelGen = 40
# StartSelGen = 21
# StopSelGen = 40
# NumberOfSires = 12
# NumberOfDams = 50
# AlphaSimDir = os.getcwd() + '/'
# selParhome['AlphaSimDir'] = os.getcwd()
# AlphaSimPed = selParhome['AlphaSimDir'] + '/SimulatedData/PedigreeAndGeneticValues.txt'
# if selParhome['EBV']:
#     seltype = 'class'
# if selParhome['gEBV']:
#     seltype = 'gen'


BurnInYN = "False"  # ali izvedeš tudi BurnIn
SelYN = "True"  # ali izvedeš tudi BurnIn
StNB = 17280
StBurnInGen = 20
StSelGen = 40
StartSelGen = 21
StopSelGen = 40
NumberOfSires = 30
NumberOfDams = 8640
AlphaSimDir = os.getcwd() + '/'
selParimport['AlphaSimDir'] = os.getcwd()
AlphaSimPed = selParimport['AlphaSimDir'] + '/SimulatedData/PedigreeAndGeneticValues.txt'
if selParimport['EBV']:
    seltype = 'class'
if selParimport['gEBV']:
    seltype = 'gen'



#selPar["stNBn"] = 8640*2
os.chdir(AlphaSimDir)
#
# selPar['pbn'] = 1
#
# selPar['genpbn'] = 1
# selPar['pozitivnoTestDoz'] = 2000
#
# selPar['pbUp'] = 1

IndCat = pd.DataFrame()

# inds = ped0.catCurrent_indiv('')
# Cats = []
# for ind in inds:
#     for (key, value) in categories.iteritems():
#         if ind in value:
#             print(key)
#
# Ind_['Indiv'] = ped.ped.Indiv
# IndCat['catBurnIN'] = ped.ped.cat


set_group_two(AlphaSimDir, "home", "import", int(selParimport["stNBn"]))
splitGenPed("PopulationSplit.txt")
ped0, c0, s0, a0 = nastavi_cat('GenPed_EBVhome.txt', externalPedName = "ExternalPedigreehome", group=True, groupNumber=0, noGroups=2, **selParhome) #nastavi cat za prvo skupino
ped1, c1, s1, a1 = nastavi_cat_TGV('GenPed_EBVimport.txt', externalPedName = "ExternalPedigreeimport", group=True, groupNumber=1, noGroups=2,  **selParimport) #nastavi cat za drugo skupino
ped0.write_pedTotal("/home/jana/PedTotalhome.txt")
ped1.write_pedTotal("/home/jana/PedTotalimport.txt")
os.system("Rscript /home/jana/Genotipi/Genotipi_CODES/Combine_PedTotals.R")
joinExternalPeds(["ExternalPedigreehome", "ExternalPedigreeimport"], AlphaSimDir)
record_groups(["home", "import"], "PopulationSplit.txt")

krogov = 2
for krog in range(krogov): #ponavljaj kolikor krogov selekcije hočeš
    # tukaj razdeli populacijo na domačo in tujo
    splitGenPed("PopulationSplit.txt")
    # tukaj izvedi celotno selekcijo v tuji populaciji --> naknadno shrani še očete z izberi_ocete_PT
    # v domači odberi in nastavi matere --> očete (za bm) uvoziš
    pedI, cI, sI, aI = selekcija_total_TGV('GenPed_EBVimport.txt', externalPedName="ExternalPedigreeimport", group=True,
                                           groupNumber=1, noGroups=2, **selParimport)
    if selParimport['EBV']:
        Oce_import = pedI.izberi_ocete_PT(selParimport["pbUp"])  # tukaj so PT testirani očetje
    if selParimport['gEBV']:
        Oce_import = pedI.izberi_ocete_gen(selParimport["pbUp"])  # tukaj so genomsko testirani očetje

    # odberi starše domače populacije
    pedH, cH, sH, aH = selekcija_importOcetov('GenPed_EBVhome.txt', externalPedName="ExternalPedigreehome", group=True,
                                              groupNumber=0, noGroups=2,
                                              importBool=True, importGroup="bm", FatherList=Oce_import, **selParhome)
    joinExternalPeds(["ExternalPedigreehome", "ExternalPedigreeimport"], AlphaSimDir)
    record_groups(["home", "import"], "PopulationSplit.txt")

    # kopiraj pedigre v selection folder
    if not os.path.exists(AlphaSimDir + '/Selection/SelectionFolder' + str(roundNo) + '/'):
        os.makedirs(AlphaSimDir + '/Selection/SelectionFolder' + str(roundNo) + '/')
    shutil.copy(AlphaSimDir + '/ExternalPedigree.txt',
                AlphaSimDir + '/Selection/SelectionFolder' + str(roundNo) + '/')
    # TUKAJ POTEM popravis AlphaSimSpec
    # PRVIc PO BURN IN-U
    SpecFile = AlphaSimSpec(os.getcwd(),
                            WorkingDir + "/CodeDir")  # AlphaSimSpec je class iz selection, ki omogoča nastavljanje parametrov AlphaSimSpec fila

    SpecFile.setPedType("ExternalPedigree.txt")
    SpecFile.setBurnInGen(StBurnInGen)
    SpecFile.setSelGen(StSelGen)
    SpecFile.setNoSires(NumberOfSires)
    SpecFile.setNoDams(NumberOfDams)
    SpecFile.turnOnGenFlex()
    SpecFile.setFlexGenToFrom((StBurnInGen + roundNo), (StBurnInGen + roundNo))
    SpecFile.turnOnSelFlex()
    SpecFile.setExtPedForGen(StBurnInGen + roundNo)
    SpecFile.setTBVComp(2)
    SpecFile.setNB(StNB)
    # pozenes ALPHASIM
    os.system('./AlphaSim1.08')
    # tukaj odstrani chip2 genotype file in izračunaj heterozigotnost na nevtralnih lokusih (chip2 - chip1)
    os.system(
        "/exports/cmvm/eddie/eb/groups/tier2_hickey_external/R-3.4.2/bin/Rscript MeanHetChip2_NeutralMarker.R " + str(
            roundNo + 20) + " " + str(rep) + " " + str(scenario) + " " + str(strategy))
    os.system("bash ChangeChip2Geno_IDs.sh")

    # tukaj dodaj kategorije k PedigreeAndGeneticValues (AlphaSim File)
    PedCat = OrigPed(AlphaSimDir, WorkingDir + '/CodeDir')
    PedCat.addInfo()  # to ti zapiše PedigreeAndGeneticValues_cat.txt v AlphaSim/SimualatedData

    # tukaj pridobi podatke za generacijske intervale
    GenInt = genInterval(AlphaSimDir)  # tukaj preberi celoten pedigre
    if seltype == 'class':
        GenInt.prepareGenInts(
            ['vhlevljeni', 'pt'])  # pri klasični so izrbrani potomci vhlevljeni (test in pripust) in plemenske telice
    if seltype == 'gen':
        GenInt.prepareGenInts(['genTest',
                               'pt'])  # pri klasični so izbrani potomci vsi genomsko testirani (pozTest in pripust) in plemenske telice
    blupNextGen = estimateBV(AlphaSimDir, WorkingDir + "/CodeDir", way='milk', sel=seltype)
    blupNextGen.computeEBV()  # estimate EBV with added phenotypes only of animals of certain category (here = milk)
    # Acc.saveAcc()
    # GenTrends.saveTrends()
    # zdaj za vsako zapiši, ker vsakič na novo prebereš
    # Acc.writeAcc()
    # GenTrends.writeTrends()
    
    
    
    splitGenPed("PopulationSplit.txt")
    pedH, cH, sH, aH = selekcija_total('GenPed_EBVhome.txt', externalPedName = "ExternalPedigreehome", group=True, groupNumber=0, noGroups=2,  **selParhome)
    pedI, cI, sI, aI = selekcija_total_TGV('GenPed_EBVimport.txt', externalPedName = "ExternalPedigreeimport", group=True, groupNumber=1, noGroups=2, **selParimport)
    pedH.write_pedTotal("/home/jana/PedTotalhome.txt")
    pedI.write_pedTotal("/home/jana/PedTotalimport.txt")
    os.system("Rscript /home/jana/Genotipi/Genotipi_CODES/Combine_PedTotals.R")
    joinExternalPeds(["ExternalPedigreehome", "ExternalPedigreeimport"], AlphaSimDir)
    record_groups(["home", "import"], "PopulationSplit.txt")


#test import
    splitGenPed("PopulationSplit.txt")
    pedI, cI, sI, aI = selekcija_total_TGV('GenPed_EBVimport.txt', externalPedName="ExternalPedigreeimport", group=True, groupNumber=1, noGroups=2, **selParimport)
    if selParimport['EBV']:
        Oce_import = pedI.izberi_ocete_PT(selParimport["pbUp"]) #tukaj so PT testirani očetje
    if selParimport['gEBV']:
        Oce_import = pedI.izberi_ocete_gen(selParimport["pbUp"])  # tukaj so genomsko testirani očetje

    #odberi starše domače populacije
    pedH, cH, sH, aH = selekcija_importOcetov('GenPed_EBVhome.txt', externalPedName="ExternalPedigreehome", group=True, groupNumber=0, noGroups=2,
                                           importBool=True, importGroup="bm", FatherList=Oce_import, **selParhome)
    pedH.write_pedTotal("/home/jana/PedTotalhome.txt")
    pedI.write_pedTotal("/home/jana/PedTotalimport.txt")
    os.system("Rscript /home/jana/Genotipi/Genotipi_CODES/Combine_PedTotals.R")

    joinExternalPeds(["ExternalPedigreehome", "ExternalPedigreeimport"], AlphaSimDir)
    record_groups(["home", "import"], "PopulationSplit.txt")