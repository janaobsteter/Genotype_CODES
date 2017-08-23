# -*- coding: utf-8 -*-
from __future__ import division
import sys
import os
from PyQt4 import QtGui, QtCore, uic
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import math
import selection10
from selection10 import *
from selection10 import nastavi_cat, selekcija_total
from collections import defaultdict
import shutil
import pandas as pd
import numpy as np
import resource


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
os.chdir("/home/jana/bin/AlphaSim1.05Linux/EddieScripts/")
scenarios = ['Clas', 'GenSLO', 'BmGen', 'OtherCowsGen', 'Gen']
replicates = 1

for scenario in scenarios:
    par = pd.read_csv("~/SelectionParam_" + scenario + ".csv", header=None, names=["Keys", "Vals"])
    par.to_dict()
    selPar = defaultdict()
    for key, val in zip(par.Keys, par.Vals):
        if key not in ['BurnInYN','EBV','gEBV','PA', 'AlphaSimDir', 'genotyped', 'EliteDamsPTBulls', 'EliteDamsPABulls', 'UpdateGenRef', 'sexToUpdate', 'EliteDamsGenBulls', 'gpb_pb', 'genTest_mladi', 'genTest_gpb']:
            try:
                selPar[key] = int(val)
            except:
                selPar[key] = float(val)
        if key in  ['BurnInYN','EBV','gEBV','PA', 'AlphaSimDir', 'genotyped', 'EliteDamsPTBulls', 'EliteDamsPABulls', 'UpdateGenRef', 'sexToUpdate', 'EliteDamsGenBulls', 'gpb_pb', 'genTest_mladi', 'genTest_gpb']:
            if val in ['False', 'True']:
                selPar[key] = bool(val == 'True')
            else:
                selPar[key] = val




    for rep in range(replicates):
        os.makedirs(scenario + str(rep))
        BurnInYN = "False" #ali izvedeš tudi BurnIn
        SelYN = "True" #ali izvedeš tudi BurnIn
        StNB = 8640
        StBurnInGen = 20
        StSelGen = 40
        StartSelGen = 21
        StopSelGen = 40
        NumberOfSires = 4
        NumberOfDams = 3500
        selPar['AlphaSimDir'] = os.getcwd() + '/'
        AlphaSimDir = os.getcwd() + '/'
        AlphaSimPed = selPar['AlphaSimDir'] + '/SimulatedData/PedigreeAndGeneticValues.txt'
        if selPar['EBV']:
            seltype = 'class'
        if selPar['gEBV']:
            seltype = 'gen'


        os.system('cp -r Essentials/* . && cp -r REALFillIn20BurnIn20/* .')
        os.system('mv IndForGeno_Gen.txt IndForGen.txt')

        ##############################################################################
        #SELEKCIJA
        ##############################################################################

        for roundNo in range(1): #za vsak krog selekcije
            # prestavi se v AlphaSim Dir
            os.chdir(AlphaSimDir)
            if not os.path.isfile(AlphaSimDir + 'ReferenceSize.txt') and os.path.isfile(AlphaSimDir + "IndForGeno.txt"):
                os.system("less IndForGeno.txt | wc -l > ReferenceSize.txt")

            # Štartaj že po 20 gen kalsične selekcije
            Acc = accuracies(AlphaSimDir)
            GenTrends = TBVCat(AlphaSimDir)
            # izvedi selekcijo, doloci kategorije zivali, dodaj novo generacijo in dodeli starse
            # pedigre se zapise v AlphaSimDir/SelectionFolder/ExternalPedigree.txt
            selekcija_total('GenPed_EBV.txt', **selPar)

            # kopiraj pedigre v selection folder
            if not os.path.exists(AlphaSimDir + '/Selection/SelectionFolder' + str(roundNo) + '/'):
                os.makedirs(AlphaSimDir + '/Selection/SelectionFolder' + str(roundNo) + '/')
            shutil.copy(AlphaSimDir + '/ExternalPedigree.txt',
                        AlphaSimDir + '/Selection/SelectionFolder' + str(roundNo) + '/')
            # TUKAJ POTEM popravis AlphaSimSpec
            # PRVIc PO BURN IN-U
            SpecFile = AlphaSimSpec(AlphaSimDir, "./CodeDir")  # AlphaSimSpec je class iz selection, ki omogoča nastavljanje parametrov AlphaSimSpec fila
            shutil.copy(SpecFile.genSpecFile,
                        AlphaSimDir)  # skopiraj generično ALphaSimSpec datoteko v AlphaSimDir
            SpecFile.setPedType("ExternalPedigree.txt")
            SpecFile.setBurnInGen(StBurnInGen)
            SpecFile.setSelGen(StSelGen)
            SpecFile.setNoSires(0)
            SpecFile.setNoDams(0)
            SpecFile.turnOnGenFlex()
            SpecFile.setFlexGenToFrom((StBurnInGen + roundNo), (StBurnInGen + roundNo))
            SpecFile.turnOnSelFlex()
            SpecFile.setExtPedForGen(StBurnInGen + roundNo)
            SpecFile.setTBVComp(2)
            SpecFile.setNB(StNB)
            # pozenes ALPHASIM
            os.system('./AlphaSim1.05')

            # tukaj dodaj kategorije k PedigreeAndGeneticValues (AlphaSim File)
            PedCat = OrigPed(AlphaSimDir, './CodeDir')
            PedCat.addInfo() #to ti zapiše PedigreeAndGeneticValues_cat.txt v AlphaSim/SimualatedData

            #tukaj pridobi podatke za generacijske intervale
            GenInt = genInterval(AlphaSimDir) #tukaj preberi celoten pedigre
            if seltype == 'class':
                GenInt.prepareGenInts(['vhlevljeni', 'pt']) #pri klasični so izrbrani potomci vhlevljeni (test in pripust) in plemenske telice
            if seltype == 'gen':
                GenInt.prepareGenInts(['genTest', 'pt']) #pri klasični so izbrani potomci vsi genomsko testirani (pozTest in pripust) in plemenske telice
            blupNextGen = estimateBV(AlphaSimDir, "./CodeDir",  way='milk', sel=seltype)
            blupNextGen.computeEBV() #estimate EBV with added phenotypes only of animals of certain category (here = milk)
            Acc.saveAcc()
            GenTrends.saveTrends()
            #zdaj za vsako zapiši, ker vsakič na novo prebereš
            Acc.writeAcc()
            GenTrends.writeTrends()

        os.system('rm -rf Chromosomes Selection && cp * ' + scenario + str(rep))
        os.system('mv SimulatedData ' + scenario + str(rep))
        os.system('ls |  grep -v Eddie_SelectionRep.py | xargs rm')