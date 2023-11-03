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


WorkingDir = "/home/v1jobste/JanaO/"

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
#scenarios = ['Class_LargePop'] #, 'GenSLO', 'BmGen', 'OtherCowsGen', 'Gen']
REP = sys.argv[1]



for rep in [REP]:
    #####################################################################################################
    #####################################################################################################
    #FILL IN - 2x newborns
    #####################################################################################################
    if not os.path.isdir("FillInBurnIn_TwoPop_" + str(rep)):
        os.makedirs("FillInBurnIn_TwoPop_" + str(rep))
    os.chdir("FillInBurnIn" + str(rep)) #prestavi se v FillInBurnin za ta replikat
    os.system('cp -r ' + WorkingDir + '/Essentials/* .') # skopiraj vse iz Esentials
    os.system('cp -r ' + WorkingDir + '/CodeDir/* .') # skopiraj vse iz CodeDir
    seed =  randint(-100000000, -1)
    os.system("echo " + str(seed) + " > Seed.txt")

    #first make a FILLIN
    #nastavi AlphaSimSpec
    print(os.getcwd())

    SpecFile = AlphaSimSpec(os.getcwd(), WorkingDir + "/CodeDir")  # AlphaSimSpec je class iz selection, ki omogoča nastavljanje parametrov AlphaSimSpec fila
    SpecFile.setPedType("Internal")  # pedigree je za burn in internal
    SpecFile.setNB(17280)  # stevilo novorojenih
    SpecFile.setBurnInGen(20)  # stevilo burnINGen
    SpecFile.setSelGen(40)  # st selection gen
    SpecFile.setNoSires(20)
    SpecFile.setNoDams(8640)
    SpecFile.turnOnGenFlex()
    SpecFile.setFlexGenToFrom(1, 21)  # pozeni od generacije 1 do burnin+1
    SpecFile.turnOnSelFlex()
    SpecFile.setExtPedForGen(21)  # za katero generacijo uvozi external pedigre  - ena po burn in
    SpecFile.setTBVComp(1)  # skomputiraj TBV
    # pozenes ALPHASIM
    os.system('./AlphaSim1.08')


    #####################################################################################################
    #####################################################################################################
        #THEN MAKE A BURN IN - classical selection!
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

     #tukaj pa še parametri za "large" population
    parimport = pd.read_csv(os.getcwd() + "/SelectionParam_Class_LargePop.csv", header=None, names=["Keys", "Vals"])
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
                   'genTest_mladi', 'genTest_gpb','genFemale']:
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
        seltype = 'class'
    if selParimport['gEBV']:
        seltype = 'gen'




    ##############################################################################
    # SELEKCIJA - 20 krogov klasične selekcije
    ##############################################################################

    for roundNo in range(1,21):  # za vsak krog selekcije
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

            Acc = accuracies(AlphaSimDir)  # nastavi
            #GenTrends = TBVCat(AlphaSimDir)
            # nimaš GenPed_EBV.txt
            blups = estimateBV(AlphaSimDir, WorkingDir + "/CodeDir", way='burnin_milk', sel='class')
            blups.computeEBV()  # tukaj izbriši samo fenotipe moških - ne morš po kategorijah, ker jih nimaš
            # Acc.saveAcc()
            #tukaj sedaj določi populacije
            set_group_two(AlphaSimDir, "home", "import", int(selParhome["stNBn"]))
            splitGenPed("PopulationSplit.txt")
            nastavi_cat('GenPed_EBVhome.txt', externalPedName = "ExternalPedigreehome", group=True, groupNumber=0, numberOfGroups = 1, **selParhome) #nastavi cat za prvo skupino
            nastavi_cat_TGV('GenPed_EBVimport.txt', externalPedName = "ExternalPedigreeimport", group=True, groupNumber=1, numberOfGroups = 2, **selParimport) #nastavi cat za drugo skupino
            joinExternalPeds(["ExternalPedigreehome", "ExternalPedigreeimport"], AlphaSimDir)
            record_groups(["home", "import"],"PopulationSplit.txt")

        else:
            Acc = accuracies(AlphaSimDir)
            GenTrends = TBVCat(AlphaSimDir)
            # izvedi selekcijo, doloci kategorije zivali, dodaj novo generacijo in dodeli starse
            # pedigre se zapise v AlphaSimDir/SelectionFolder/ExternalPedigree.txt
            splitGenPed("PopulationSplit.txt")
            selekcija_total('GenPed_EBVhome.txt', externalPedName = "ExternalPedigreehome", group=True, groupNumber=0, numberOfGroups = 2, **selParhome)
            selekcija_total_TGV('GenPed_EBVimport.txt', externalPedName = "ExternalPedigreeimport", group=True, groupNumber=1, numberOfGroups = 2, **selParimport)
            joinExternalPeds(["ExternalPedigreehome", "ExternalPedigreeimport"], AlphaSimDir)
            record_groups(["home", "import"],"PopulationSplit.txt")

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
        SpecFile.setNoSires(35)
        SpecFile.setNoDams(4320)
        SpecFile.turnOnGenFlex()
        SpecFile.setFlexGenToFrom((StBurnInGen + roundNo), (StBurnInGen + roundNo))
        SpecFile.turnOnSelFlex()
        SpecFile.setExtPedForGen(StBurnInGen + roundNo)
        SpecFile.setTBVComp(2)
        SpecFile.setNB(StNB)
        # pozenes ALPHASIM
        os.system('./AlphaSim1.08')
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
        blupNextGen = estimateBV(AlphaSimDir, WorkingDir + "/CodeDir", way='milk', sel=seltype)
        blupNextGen.computeEBV()  # estimate EBV with added phenotypes only of animals of certain category (here = milk)
        Acc.saveAcc()
        #GenTrends.saveTrends()
        # zdaj za vsako zapiši, ker vsakič na novo prebereš
        Acc.writeAcc()
        #GenTrends.writeTrends()



