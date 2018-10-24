# -*- coding: utf-8 -*-
from __future__ import division
import sys
import os
import math
from selection10 import *
from selection10 import nastavi_cat, selekcija_total
from collections import defaultdict
import shutil
import pandas as pd
import numpy as np
import resource
import ast
WorkingDir = os.getcwd()

#sys.argv: 1 = rep, 2 = scenario, 3 = strategy, 4 = reference size

class estimateBV:
    def __init__(self, AlphaSimDir, codeDir, way, sel):
        self.AlphaSimDir = AlphaSimDir
        self.way = way
        self.sel = sel
        self.codeDir = codeDir

    def computeEBV(self):
        # pripravi fajle za blupf90
        blupFiles = blupf90(self.AlphaSimDir, self.codeDir, way=self.way)
        # listUnphenotyped = ['potomciNP', 'nr', 'telF', 'telM', 'pt', 'mladi', 'vhlevljeni', 'cak'] #list of unphenotyped categories (better ages?)
        # blupFiles.preparePedDat_cat(listUnphenotyped) #pripravi ped, dat file za blup #skopiraj generičen paramfile v AlphaSim Directory
        if self.way == 'milk':
            blupFiles.makeDat_removePhen_milk()  # odstrani genotip moškim živali in pripravi dat file - repetabaility model
        if self.way == 'burnin_milk':  # če je takoj po burn inu - nimaš še kategorij (za prvih n burningeneracij brze dodane naslednje)
            blupFiles.makeDat_sex(2)
        # skopiraj paramFile za renumf90
        if self.sel == 'gen':
            shutil.copy(blupFiles.blupgenParamFile,
                        blupFiles.AlphaSimDir + 'renumf90.par')  # skopiraj template blupparam file
        if self.sel == 'class':
            shutil.copy(blupFiles.blupgenParamFile_Clas,
                        blupFiles.AlphaSimDir + 'renumf90.par')  # skopiraj template blupparam file

        # uredi blupparam file
        # get variance components from AlphaSim Output Files
        OutputFiles = AlphaSim_OutputFile(self.AlphaSimDir)
        genvar = OutputFiles.getAddVar()  # dobi additivno varianco
        resvar = OutputFiles.getResVar()  # dobi varianco za ostanek

        blupFiles.prepareParamFiles(genvar, resvar,
                                    self.AlphaSimDir + '/renumf90.par')  # set levels of random aniaml effect, add var and res var
        # the paramfile is now set
        blupFiles.makePed_gen()  # make ped file for blup, no Code!
        if self.sel == 'gen':
            GenFiles = snpFiles(self.AlphaSimDir)
            GenFiles.createBlupf90SNPFile()

        os.system("./renumf90 < renumParam")  # run renumf90

        # if self.sel == 'class': #ZDJ TEGA NI, KER JE POSEBEJ FILE!
        # if self.sel == 'class': #ZDJ TEGA NI, KER JE POSEBEJ FILE!
        # os.system("head -n-3 renf90.par > tmp && mv tmp renf90.par")
        # os.system("./blupf90 blupf90_Selection")  # run blupf90

        resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
        # os.system('./preGSf90 renf90.par')
        os.system('./blupf90 renf90.par')
        # os.system('./postGSf90 renf90.par')


        # renumber the solutions
        # copy the solution in a file that does not get overwritten
        os.system("bash Match_AFTERRenum.sh")
        shutil.copy('renumbered_Solutions', 'renumbered_Solutions_' + str(blupFiles.gen))
        # shutil.copy('solutions', 'renumbered_Solutions_' + str(blupFiles.gen))

        blupFiles.prepareSelPed()  # obtain solution and add them to
        # AlphaPed PedigreeAndGeneticValues files --> Write them to GenPed_EBV.txt, which is read by module selection


######################################################################################
######################################################################################
#os.chdir("/home/jana/bin/AlphaSim1.05Linux/EddieScripts/")
#argument 0 is the name of the script
rep = sys.argv[1]
scenario = sys.argv[2]
strategy = sys.argv[3]
refSize = sys.argv[4]

os.chdir(strategy + "/")

print("Creating directory " + scenario + str(rep))
if not os.path.isdir(scenario + str(rep)):
    os.makedirs(scenario + str(rep))
SelectionDir = scenario + str(rep) + "/"





######################################################################################################
#now run all the scenarios
######################################################################################################
#potem se prestavi nazaj v working directory
os.chdir(SelectionDir)

print("Copying files to " + SelectionDir)
os.system('cp -r ' + WorkingDir + '/FillInBurnIn' + str(rep) + '/* .')
os.system('cp -r ' + WorkingDir + '/Essentials/* .')
os.system('cp -r ' + WorkingDir + '/CodeDir/* .')
os.system('mv IndForGeno_' + refSize + '.txt IndForGeno.txt')

os.system("chmod a+x AlphaSim1.08")
os.system("chmod a+x renumf90")
os.system("chmod a+x blupf90")

parhome = pd.read_csv(os.getcwd() + "/SelectionParam_Class.csv", header=None, names=["Keys", "Vals"])
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

# tukaj pa še parametri za "large" population
parimport = pd.read_csv(os.getcwd() + "/SelectionParam_" + scenario + "_LargePop.csv", header=None, names=["Keys", "Vals"])
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
NumberOfSires = 30
NumberOfDams = 8640
AlphaSimDir = os.getcwd() + '/'
selParimport['AlphaSimDir'] = os.getcwd()
if selParimport['EBV']:
    seltypeimport = 'class'
if selParimport['gEBV']:
    seltype = 'gen'

##############################################################################
#SELEKCIJA
##############################################################################
print(AlphaSimDir)
for roundNo in range(21,41): #za vsak krog selekcije
    # prestavi se v AlphaSim Dir
    if not os.path.isfile(AlphaSimDir + 'ReferenceSize.txt') and os.path.isfile(AlphaSimDir + "IndForGeno.txt"):
        os.system("less IndForGeno.txt | wc -l > ReferenceSize.txt")

    # Štartaj že po 20 gen kalsične selekcije
    Acc = accuracies(AlphaSimDir)
    #GenTrends = TBVCat(AlphaSimDir)
    # izvedi selekcijo, doloci kategorije zivali, dodaj novo generacijo in dodeli starse
    # pedigre se zapise v AlphaSimDir/SelectionFolder/ExternalPedigree.txt

    #tukaj razdeli populacijo na domačo in tujo
    splitGenPed("PopulationSplit.txt")
    #tukaj izvedi celotno selekcijo v tuji populaciji --> naknadno shrani še očete z izberi_ocete_PT
    #v domači odberi in nastavi matere --> očete (za bm) uvoziš
    pedI, cI, sI, aI = selekcija_total_TGV('GenPed_EBVimport.txt', externalPedName="ExternalPedigreeimport", group=True, groupNumber=1,
                        **selParimport)
    if selParimport['EBV']:
        Oce_import = pedI.izberi_ocete_PT(selParimport["pbUp"]) #tukaj so PT testirani očetje
    if selParimport['gEBV']:
        Oce_import = pedI.izberi_ocete_gen(selParimport["pbUp"])  # tukaj so genomsko testirani očetje

    #odberi starše domače populacije
    pedH, cH, sH, aH = selekcija_importOcetov('GenPed_EBVhome.txt', externalPedName="ExternalPedigreehome", group=True, groupNumber=0,
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
    SpecFile = AlphaSimSpec(os.getcwd(),WorkingDir + "/CodeDir")  # AlphaSimSpec je class iz selection, ki omogoča nastavljanje parametrov AlphaSimSpec fila

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
    #tukaj odstrani chip2 genotype file in izračunaj heterozigotnost na nevtralnih lokusih (chip2 - chip1)
    os.system("/exports/cmvm/eddie/eb/groups/tier2_hickey_external/R-3.4.2/bin/Rscript MeanHetChip2_NeutralMarker.R " + str(roundNo+20) + " " + str(rep) + " " + str(scenario) + " " + str(strategy))			
    os.system("bash ChangeChip2Geno_IDs.sh")    	
	
    # tukaj dodaj kategorije k PedigreeAndGeneticValues (AlphaSim File)
    PedCat = OrigPed(AlphaSimDir, WorkingDir + '/CodeDir')
    PedCat.addInfo() #to ti zapiše PedigreeAndGeneticValues_cat.txt v AlphaSim/SimualatedData

    #tukaj pridobi podatke za generacijske intervale
    GenInt = genInterval(AlphaSimDir) #tukaj preberi celoten pedigre
    if seltype == 'class':
        GenInt.prepareGenInts(['vhlevljeni', 'pt']) #pri klasični so izrbrani potomci vhlevljeni (test in pripust) in plemenske telice
    if seltype == 'gen':
        GenInt.prepareGenInts(['genTest', 'pt']) #pri klasični so izbrani potomci vsi genomsko testirani (pozTest in pripust) in plemenske telice
    blupNextGen = estimateBV(AlphaSimDir, WorkingDir + "/CodeDir",  way='milk', sel=seltype)
    blupNextGen.computeEBV() #estimate EBV with added phenotypes only of animals of certain category (here = milk)
    #Acc.saveAcc()
    #GenTrends.saveTrends()
    #zdaj za vsako zapiši, ker vsakič na novo prebereš
    #Acc.writeAcc()
    #GenTrends.writeTrends()


os.system('rm -rf Chromosomes Selection && cp * ' + scenario + str(rep))
#os.system('rm SimulatedData/UnrestrictedQtnIndivGenotypes.txt')
os.system('rm SimulatedData/RestrictedQtnIndivGenotypes.txt')

