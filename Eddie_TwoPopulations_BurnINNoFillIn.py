# -*- coding: utf-8 -*-
from __future__ import division
import sys
import os
import math
from selection10 import *
from selection10 import nastavi_cat_TGV, selekcija_total_TGV
from collections import defaultdict
import shutil
import pandas as pd
import numpy as np
import resource
import ast
from random import randint

WorkingDir = "/home/v1jobste/jobsteter/"


class estimateBV:
    def __init__(self, AlphaSimDir, codeDir, way, sel):
        self.AlphaSimDir = AlphaSimDir
        self.way = way
        self.sel = sel
        self.codeDir = codeDir

    def computeEBV(self, group = None, dataGroup = False, prepareSelPed = True, traitEBV = 1, multipleTraitsTBV = None):
        """
        A function to prepare input files, estimate breeding values with blupf90 and prepare output files
        :param group: The group you are estimating the breeding values for
        :param dataGroup: Use only group data for the estimation
        :param prepareSelPed: Do you want to prepare selection ped = GenPed_EBV in this function
        :return:
        """
        # pripravi fajle za blupf90
        blupFiles = blupf90(self.AlphaSimDir, self.codeDir, way=self.way, trait = traitEBV)
        # listUnphenotyped = ['potomciNP', 'nr', 'telF', 'telM', 'pt', 'mladi', 'vhlevljeni', 'cak'] #list of unphenotyped categories (better ages?)
        # blupFiles.preparePedDat_cat(listUnphenotyped) #pripravi ped, dat file za blup #skopiraj generičen paramfile v AlphaSim Directory
        if self.way == 'milk':
            blupFiles.makeDat_removePhen_milk()  # odstrani genotip moškim živali in pripravi dat file - repetabaility model
        if self.way == 'burnin_milk':  # če je takoj po burn inu - nimaš še kategorij (za prvih n burningeneracij brze dodane naslednje)
            blupFiles.makeDat_sex(2)

        #make ped files
        blupFiles.makePed_gen()  # make ped file for blup, no Code!

        # skopiraj paramFile za renumf90
        if not dataGroup:
            if self.sel == 'gen':
                shutil.copy(blupFiles.blupgenParamFile,
                            blupFiles.AlphaSimDir + 'renumf90.par')  # skopiraj template blupparam file
                #prepare genoFile
                os.system(
                    "sed -i 's/Blupf90.dat/Blupf90_trait" + str(traitEBV) + ".dat/g' renumf90.par")
                GenFiles = snpFiles(self.AlphaSimDir)
                GenFiles.createBlupf90SNPFile()
            elif self.sel == 'class':
                shutil.copy(blupFiles.blupgenParamFile_Clas,
                            blupFiles.AlphaSimDir + 'renumf90.par')  # skopiraj template blupparam file
                os.system(
                    "sed -i 's/Blupf90.dat/Blupf90_trait" + str(traitEBV) + ".dat/g' renumf90.par")

        elif dataGroup:
            #first split the .dat file (for the group)
            blupFiles.popSplitDat('Blupf90_trait' + str(traitEBV) + '.dat', "PopulationSplit.txt")
            if self.sel == 'class':
                shutil.copy(blupFiles.blupgenParamFile_Clas_group,
                            blupFiles.AlphaSimDir + 'renumf90.par')  # skopiraj template blupparam file
                os.system("sed -i 's/Blupf90_group.dat/Blupf90_trait" + str(traitEBV) + "_" + group + ".dat/g' renumf90.par")
            if self.sel == 'gen':
                shutil.copy(blupFiles.blupgenParamFile_group,
                            blupFiles.AlphaSimDir + 'renumf90.par')  # skopiraj template blupparam file
                #prepare genoFile
                GenFiles = snpFiles(self.AlphaSimDir)
                GenFiles.createBlupf90SNPFile(group = group)
                os.system("sed -i 's/Blupf90_group.dat/Blupf90_trait" + str(traitEBV) + "_" + group + ".dat/g' renumf90.par")
                os.system("sed -i 's/GenoFile_group.txt/GenoFile" + group +".txt/g' renumf90.par")


        # uredi blupparam file
        # get variance components from AlphaSim Output Files
        OutputFiles = AlphaSim_OutputFile(self.AlphaSimDir)
        genvar = OutputFiles.getAddVar(trait=traitEBV)  # dobi additivno varianco
        resvar = OutputFiles.getResVar(trait=traitEBV)  # dobi varianco za ostanek

        blupFiles.prepareParamFiles(genvar, resvar,
                                    self.AlphaSimDir + '/renumf90.par')  # set levels of random aniaml effect, add var and res var

        #This is just to track wqhether everything is ok with the parameters files and dat files
        if group:
            os.system("cp renumf90.par renumf90_" + group + "_trait_" + str(traitEBV) + ".par")

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
        if not group:
            shutil.copy('renumbered_Solutions', 'renumbered_Solutions_' + str(blupFiles.gen))
        if group:
            shutil.copy('renumbered_Solutions', 'renumbered_Solutions_' + group + '_' + str(blupFiles.gen))
        # shutil.copy('solutions', 'renumbered_Solutions_' + str(blupFiles.gen))

        if prepareSelPed:
            if not group:
                blupFiles.prepareSelPed(multipleTraits=multipleTraitsTBV)  # obtain solution and add them to
                # AlphaPed PedigreeAndGeneticValues files --> Write them to GenPed_EBV.txt, which is read by module selection
            elif group:
                blupFiles.prepareSelPed_group("PopulationSplit.txt", multipleTraits=multipleTraitsTBV)


######################################################################################
######################################################################################

rep = sys.argv[1]
scenario = sys.argv[2]
strategy = sys.argv[3]
refSize = sys.argv[4]
traitHome = int(sys.argv[5])
traitImport = int(sys.argv[6])

#####################################################################################################
#####################################################################################################
# FILL IN - 2x newborns
#####################################################################################################
if not os.path.isdir("BurnIn_TwoPop_" + str(rep)+ "_" + str(traitHome) + str(traitImport)):
    os.makedirs("BurnIn_TwoPop_" + str(rep) + "_" +  str(traitHome) + str(traitImport))
os.chdir("BurnIn_TwoPop_" + str(rep) + "_" +  str(traitHome) + str(traitImport))  # prestavi se v FillInBurnin za ta replikat
os.system('cp -r ' + WorkingDir + '/FillIn_TwoPop_' + str(rep) + '/* .')  # skopiraj vse iz Esentials
os.system('cp -r ' + WorkingDir + '/Essentials/* .')  # skopiraj vse iz Esentials
#os.system('cp -r ' + WorkingDir + '/CodeDir/* .')  # skopiraj vse iz CodeDir


#####################################################################################################
#####################################################################################################
# THEN MAKE A BURN IN - classical selection!
parhome = pd.read_csv(WorkingDir + "/SelPar/10K/SU55SelPar/SelectionParam_" + scenario + ".csv", header=None, names=["Keys", "Vals"])
parhome.to_dict()
selParhome = defaultdict()
for key, val in zip(parhome.Keys, parhome.Vals):
    if key not in ['BurnInYN', 'EBV', 'gEBV', 'PA', 'AlphaSimDir', 'genotyped', 'EliteDamsPTBulls',
                   'EliteDamsPABulls', 'UpdateGenRef', 'sexToUpdate', 'EliteDamsGenBulls', 'gpb_pb',
                   'genTest_mladi', 'genTest_gpb', 'importPer', 'genFemale']:
        try:
            selParhome[key] = int(val)
        except:
            selParhome[key] = float(val)
    if key in ['BurnInYN', 'EBV', 'gEBV', 'PA', 'AlphaSimDir', 'EliteDamsPTBulls',
               'EliteDamsPABulls', 'UpdateGenRef', 'sexToUpdate', 'EliteDamsGenBulls', 'gpb_pb',
               'genTest_mladi', 'genTest_gpb']:
        if val in ['False', 'True']:
            selParhome[key] = bool(val == 'True')
        else:
            selParhome[key] = val
    if key == 'genotyped':
        selParhome[key] = ast.literal_eval(val)

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

# tukaj pa še parametri za "large" population
parimport = pd.read_csv(WorkingDir + "/SelPar/10K/SU55SelPar/SelectionParam_" + scenario + "_LargePop.csv", header=None, names=["Keys", "Vals"])
parimport.to_dict()
selParimport = defaultdict()
for key, val in zip(parimport.Keys, parimport.Vals):
    if key not in ['BurnInYN', 'EBV', 'gEBV', 'PA', 'AlphaSimDir', 'genotyped', 'EliteDamsPTBulls',
                   'EliteDamsPABulls', 'UpdateGenRef', 'sexToUpdate', 'EliteDamsGenBulls', 'gpb_pb',
                   'genTest_mladi', 'genTest_gpb', 'importPer', 'genFemale']:
        try:
            selParimport[key] = int(val)
        except:
            selParimport[key] = float(val)
    if key in ['BurnInYN', 'EBV', 'gEBV', 'PA', 'AlphaSimDir', 'EliteDamsPTBulls',
               'EliteDamsPABulls', 'UpdateGenRef', 'sexToUpdate', 'EliteDamsGenBulls', 'gpb_pb',
               'genTest_mladi', 'genTest_gpb']:
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
AlphaSimPed = selParimport['AlphaSimDir'] + '/SimulatedData/PedigreeAndGeneticValues.txt'
if selParimport['EBV']:
    seltype = 'class'
if selParimport['gEBV']:
    seltype = 'gen'

##############################################################################
# SELEKCIJA - 20 krogov klasične selekcije
##############################################################################

for roundNo in range(1, 21):
    if roundNo == 1:  # če je to prvi krog - nimaš še kategorij od prej, nimaš niti EBV-jev
        # odstrani Blupf90 fajle iz prejšnjih runov - ker se merge-a
        # enako tudi za generacijski interval in file z genotipi
        if os.path.isfile(AlphaSimDir + 'Blupf90.dat'):
            os.remove(AlphaSimDir + 'Blupf90.dat')
        if os.path.isfile(AlphaSimDir + 'GenInts.txt'):
            os.remove(AlphaSimDir + 'GenInts.txt')
        if os.path.isfile(AlphaSimDir + 'GenoFile.txt'):
            os.remove(AlphaSimDir + 'GenoFile.txt')
        if os.path.isfile(AlphaSimDir + 'IndForGeno.txt'):
            os.remove(AlphaSimDir + 'IndForGeno.txt')
        if os.path.isfile(AlphaSimDir + 'GenTrends_gen.csv'):
            os.remove(AlphaSimDir + 'GenTrends_gen.csv')
        if os.path.isfile(AlphaSimDir + 'GenTrends_cat.csv'):
            os.remove(AlphaSimDir + 'GenTrends_cat.csv')
        if os.path.isfile(AlphaSimDir + 'Accuracies_Cat.csv'):
            os.remove(AlphaSimDir + 'Accuracies_Cat.csv')
        if os.path.isfile(AlphaSimDir + 'Accuracies_Gen.csv'):
            os.remove(AlphaSimDir + 'Accuracies_Gen.csv')
            # if os.path.isfile(self.AlphaSimDir + 'AccuraciesBV.csv'):
            # os.remove(self.AlphaSimDir + 'AccuraciesBV.csv')

        AccHome = accuracies(AlphaSimDir, group='home', trait = traitHome)
        AccImport = accuracies(AlphaSimDir, group='import', trait = traitImport)
        # GenTrends = TBVCat(AlphaSimDir)
        # nimaš GenPed_EBV.txt
        print("Compute initial EBVs")
        blups = estimateBV(AlphaSimDir, WorkingDir + "/CodeDir", way='burnin_milk', sel='class')
        # estimate EBVs for home population with only domestic data
        blups.computeEBV(traitEBV = traitHome, multipleTraitsTBV = [traitHome,traitImport])
        # Acc.saveAcc()
        # tukaj sedaj določi populacije
        set_group_two(AlphaSimDir, "home", "import", int(selParimport["stNBn"]))
        splitGenPed("PopulationSplit.txt")
        nastavi_cat('GenPed_EBVhome.txt', externalPedName="ExternalPedigreehome", group=True, groupNumber=0, noGroups = 2,
                    **selParhome)  # nastavi cat za prvo skupino
        nastavi_cat_TGV('GenPed_EBVimport.txt', externalPedName="ExternalPedigreeimport", group=True, groupNumber=1, noGroups = 2, trait = traitImport,
                        **selParimport)  # nastavi cat za drugo skupino
        joinExternalPeds(["ExternalPedigreehome", "ExternalPedigreeimport"], AlphaSimDir)
        record_groups(['home', 'import'],"PopulationSplit.txt")

        #os.system("rm Blup*dat")

    else:
        AccHome = accuracies(AlphaSimDir, group='home', trait=traitHome)
        AccImport = accuracies(AlphaSimDir, group='import', trait=traitImport)
        GenTrends = TBVCat(AlphaSimDir)
        # izvedi selekcijo, doloci kategorije zivali, dodaj novo generacijo in dodeli starse
        # pedigre se zapise v AlphaSimDir/SelectionFolder/ExternalPedigree.txt
        #splitGenPed("PopulationSplit.txt") This is no longer needed, since you have prepareSelPed_group
        selekcija_total('GenPed_EBVhome.txt', externalPedName="ExternalPedigreehome", group=True, groupNumber=0, noGroups = 2,
                        **selParhome)
        selekcija_total_TGV('GenPed_EBVimport.txt', externalPedName="ExternalPedigreeimport", group=True, groupNumber=1, noGroups = 2, trait = traitImport,
                            **selParimport)
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
                            WorkingDir + "/CodeDir", type="Multitrait")  # AlphaSimSpec je class iz selection, ki omogoča nastavljanje parametrov AlphaSimSpec fila
    SpecFile.setPedType("ExternalPedigree.txt")
    SpecFile.setBurnInGen(StBurnInGen)
    SpecFile.setSelGen(StSelGen)
    SpecFile.setNoSires(20)
    SpecFile.setNoDams(8640)
    SpecFile.turnOnGenFlex()
    SpecFile.setFlexGenToFrom((StBurnInGen + roundNo), (StBurnInGen + roundNo))
    SpecFile.turnOnSelFlex()
    SpecFile.setExtPedForGen(StBurnInGen + roundNo)
    SpecFile.setTBVComp(2)
    SpecFile.setNB(StNB)
    # pozenes ALPHASIM
    os.system('./AlphaSim1.08')
    #tukaj odstrani chip2 genotipe (pusti ID-je) in izračunaj heterozigotnost na nevtralnih lokusih (chip2 - chip1)
    os.system("/exports/cmvm/eddie/eb/groups/tier2_hickey_external/R-3.4.2/bin/Rscript MeanHetMarker_Neutral_QTN_import.R " + str(roundNo) + " " + str(rep) + " " + str(scenario) + str(scenario) + " " + str(strategy))
    os.system("bash ChangeChip2Geno_IDs.sh")

    # tukaj dodaj kategorije k PedigreeAndGeneticValues (AlphaSim File)
    PedCat = OrigPed(AlphaSimDir, WorkingDir + '/CodeDir')
    PedCat.addInfo()  # to ti zapiše PedigreeAndGeneticValues_cat.txt v AlphaSim/SimualatedData

    # tukaj pridobi podatke za generacijske intervale
    GenInt = genInterval(AlphaSimDir)  # tukaj preberi celoten pedigre
    if seltype == 'class':
        GenInt.prepareGenInts(['vhlevljeni',
                               'pt'])  # pri klasični so izrbrani potomci vhlevljeni (test in pripust) in plemenske telice
    if seltype == 'gen':
        GenInt.prepareGenInts(['genTest',
                               'pt'])  # pri klasični so izbrani potomci vsi genomsko testirani (pozTest in pripust) in plemenske telice
    print("Compute selection EBVs")
    blupNextGen = estimateBV(AlphaSimDir, WorkingDir + "/CodeDir", way='milk', sel=seltype)
    #estimate EBVs for home population with only domestic data
    blupNextGen.computeEBV(group = "home", dataGroup = True, prepareSelPed = False, traitEBV = traitHome, multipleTraitsTBV=[traitHome,traitImport])
    #estimate EBVs for import population, use all data, create GenPed_EBVs.txt for both groups (only once, sinve it is one populationsplit file)
    blupNextGen.computeEBV(group = "import", dataGroup = True, prepareSelPed = True, traitEBV = traitImport, multipleTraitsTBV=[traitHome,traitImport])
    AccHome.saveAcc()
    AccImport.saveAcc()
    #GenTrends.saveTrends()
    #zdaj za vsako zapiši, ker vsakič na novo prebereš
    AccHome.writeAcc()
    AccImport.writeAcc()
    #GenTrends.writeTrends()
