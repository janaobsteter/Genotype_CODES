# -*- coding: utf-8 -*-
from __future__ import division
import sys
import os
import math
import selection10
from selection10 import *
from selection10 import nastavi_cat, selekcija_total
from collections import defaultdict
import shutil
import pandas as pd
import numpy as np
import resource
import ast
from random import shuffle

WorkingDir = os.getcwd()


reload(selection10)
#sys.argv: 1 = rep, 2 = scenario

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
degree = sys.argv[3]



print("Creating directory " + scenario + str(rep) +"_" + degree + "OCS")
if not os.path.isdir(scenario + str(rep) +"_" + degree + "OCS"):
    os.makedirs(scenario + str(rep) +"_" + degree + "OCS")
SelectionDir = scenario + str(rep) +"_" + degree + "OCS/"





######################################################################################################
#now run all the scenarios
######################################################################################################
#potem se prestavi nazaj v working directory
os.chdir(SelectionDir)

print("Copying files to " + SelectionDir)
os.system('cp -r ' + WorkingDir + '/Essentials/* .')
os.system('cp -r ' + WorkingDir + '/FillInBurnIn' + str(rep) + '/* .')
os.system('cp -r ' + WorkingDir + '/CodeDir/* .')
os.system('cp ' + WorkingDir + '/AlphaMate .')
os.system('mv IndForGeno_10000.txt IndForGeno.txt')

par = pd.read_csv(WorkingDir + "/Essentials/SelectionParam_" + scenario + ".csv", header=None, names=["Keys", "Vals"])
par.to_dict()
selPar = defaultdict()
for key, val in zip(par.Keys, par.Vals):
    if key not in ['BurnInYN', 'EBV', 'gEBV', 'PA', 'AlphaSimDir', 'genotyped', 'EliteDamsPTBulls',
                   'EliteDamsPABulls', 'UpdateGenRef', 'sexToUpdate', 'EliteDamsGenBulls', 'gpb_pb',
                   'genTest_mladi', 'genTest_gpb']:
        try:
            selPar[key] = int(val)
        except:
            selPar[key] = float(val)
    if key in ['BurnInYN', 'EBV', 'gEBV', 'PA', 'AlphaSimDir', 'EliteDamsPTBulls',
               'EliteDamsPABulls', 'UpdateGenRef', 'sexToUpdate', 'EliteDamsGenBulls', 'gpb_pb',
               'genTest_mladi', 'genTest_gpb']:
        if val in ['False', 'True']:
            selPar[key] = bool(val == 'True')
        else:
            selPar[key] = val
    if key == 'genotyped':
        selPar[key] = ast.literal_eval(val)


BurnInYN = "False" #ali izvedeš tudi BurnIn
SelYN = "True" #ali izvedeš tudi BurnIn
StNB = 8640
StBurnInGen = 20
StFillInBurnIn = 40
StSelGen = 40
StartSelGen = 21
StopSelGen = 40
NumberOfSires = 12
NumberOfDams = 3500
selPar['AlphaSimDir'] = os.getcwd() + '/'
AlphaSimDir = os.getcwd() + '/'
AlphaSimPed = selPar['AlphaSimDir'] + '/SimulatedData/PedigreeAndGeneticValues.txt'
if selPar['EBV']:
    seltype = 'class'
if selPar['gEBV']:
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
    #TO JE SEDAJ OCS KORAK!!!!!!
    ped = odberiStarse_OCSgen('GenPed_EBV.txt', AlphaSimDir, **selPar)

    # prepare pedigree matrix for selected individuals
    pedA = AlphaRelate(AlphaSimDir, AlphaSimDir)
    pedA.preparePedigree()
    pedA.runAlphaRelate()

    mate = AlphaMate(AlphaSimDir, AlphaSimDir, roundNo+19)
    mate.prepareGender()
    mate.prepareCriterionFile()
  
    os.system("Rscript optiSel.R")
    Ocetje = list(pd.read_table("Ocetje.txt", header=None).loc[:, 0])
    #len(Ocetje)

    # dodaj novo generacijo
    ped.add_new_gen_naive(selPar['stNBn'], selPar['potomciNPn'] * 2)
    ped.compute_age()
    # dodaj matere
    ped.doloci_matere(selPar.get('stNBn'), selPar.get('potomciNPn'), selPar.get('ptn'), selPar.get('kraveUp'))
    # preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
    ped.mother_nr_blank()
    # dodaj očete
    print(len(Ocetje))
    if len(Ocetje) == selPar.get('stNBn'):
    	ped.ped.loc[ped.ped.cat.isin(['nr', 'potomciNP']), 'Father'] = Ocetje
	ped.ped.Father = ped.ped.Father.astype(int)
    else:
        print("Not enough fathers!!!!")

    ped.ped.loc[ped.ped.cat == "kandidati", 'cat'] = "izl"
    ped.ped.loc[ped.ped.Indiv.isin(list(set(Ocetje))), "cat"] = "gpb"


    # preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
    ped.mother_nr_blank()

    # ped.UpdateIndCat('/home/jana/')
    ped.save_cat_DF()
    ped.save_sex_DF()
    ped.save_active_DF()
    if selPar.get('gEBV'):
        if selPar.get('UpdateGenRef'):
            ped.updateAndSaveIndForGeno(selPar.get('genotyped'), selPar.get('NbUpdatedGen'), selPar.get('sexToUpdate'),
                                        selPar.get('AlphaSimDir'))
        if not selPar.get('UpdateGenRef'):
            ped.saveIndForGeno(selPar.get('genotyped'))
        os.system(
            'less IndForGeno.txt | wc -l > ReferenceSize_new.txt && cat ReferenceSize_new.txt ReferenceSize.txt > Reftmp && mv Reftmp ReferenceSize.txt')
    ped.write_ped(AlphaSimDir + "/ExternalPedigree.txt")
    ped.write_pedTotal(AlphaSimDir + "/ExternalPedigreeTotal.txt")

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
    os.system(AlphaSimDir + '/AlphaSim1.08')

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
os.system('rm SimulatedData/UnrestrictedQtnIndivGenotypes.txt')
os.system('rm SimulatedData/RestrictedQtnIndivGenotypes.txt')

