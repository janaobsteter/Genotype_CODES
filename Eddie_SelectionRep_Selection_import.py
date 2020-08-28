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
#os.chdir("/home/jana/bin/AlphaSim1.05Linux/EddieScripts/")
#argument 0 is the name of the script
rep = sys.argv[1]
scenarioHome = sys.argv[2]
scenarioImport = sys.argv[3]
strategy = sys.argv[4]
refSize = sys.argv[5]
traitHome = int(sys.argv[6])
traitImport = int(sys.argv[7])
percentageImport_k = sys.argv[8] if len(sys.argv) > 7 else 0
percentageImport_bm = sys.argv[9] if len(sys.argv) > 8 else 0


os.chdir(refSize + "/" + strategy + "_import/")

print("Creating directory " + scenarioHome + scenarioImport + str(rep) + "_" + str(percentageImport_bm))
if not os.path.isdir(scenarioHome + scenarioImport + str(rep) + "_" + str(percentageImport_bm) + str(traitHome) + str(traitImport)):
    os.makedirs(scenarioHome + scenarioImport + str(rep) + "_" + str(percentageImport_bm) + str(traitHome) + str(traitImport))
SelectionDir = scenarioHome + scenarioImport + str(rep) + "_" + str(percentageImport_bm) + str(traitHome) + str(traitImport) + "/"





######################################################################################################
#now run all the scenarios
######################################################################################################
#potem se prestavi nazaj v working directory
os.chdir(SelectionDir)

print("Copying files to " + SelectionDir)
os.system('cp -r ' + WorkingDir + '/BurnIn_TwoPop_' + str(rep) + '_' +  str(traitHome) + str(traitImport) + '/*' + ' .')
#os.system('cp -r ' + WorkingDir + '/Essentials/* .')
#os.system('cp -r ' + WorkingDir + '/CodeDir/* .')


os.system("chmod a+x AlphaSim1.08")
os.system("chmod a+x renumf90")
os.system("chmod a+x blupf90")

parhome = pd.read_csv(WorkingDir + "/SelPar/10K/SU55SelPar/SelectionParam_" + scenarioHome + ".csv", header=None, names=["Keys", "Vals"])
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

if selParhome['EBV']:
    seltype = 'class'
if selParhome['gEBV']:
    seltype = 'gen'

if percentageImport_k != 0 and percentageImport_bm != 0:
    selParhome['importPer'] = {'k': int(percentageImport_k), 'bm': int(percentageImport_bm)}

print(selParhome)
print("This is the import per: ", str(selParhome['importPer']))

# tukaj pa še parametri za "large" population
parimport = pd.read_csv(WorkingDir + "/SelPar/10K/SU55SelPar/SelectionParam_" + scenarioImport + "_LargePop.csv", header=None, names=["Keys", "Vals"])
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
StSelGen = 60
StartSelGen = 21
StopSelGen = 80
NumberOfSires = 30
NumberOfDams = 8640
AlphaSimDir = os.getcwd() + '/'
selParimport['AlphaSimDir'] = os.getcwd()
selParhome['AlphaSimDir'] = os.getcwd()
if selParimport['EBV']:
    seltypeimport = 'class'
if selParimport['gEBV']:
    seltype = 'gen'

##############################################################################
#SELEKCIJA
##############################################################################
print(AlphaSimDir)
for roundNo in range(21,41): #za vsak krog selekcije
    if roundNo == 21:
        print("Creating and initial training population of all active cows and PT bulls.")
        ped = pd.read_csv('./SimulatedData/PedigreeAndGeneticValues_cat.txt', sep='\s+')
        popSplit = pd.read_csv("PopulationSplit.txt")
        #create a reference from both populations
        pd.DataFrame({"ID" : ped.loc[ped.cat.isin(["k", "pb"]), 'Indiv']}).to_csv('IndForGeno.txt', index=None, header=None, sep='\n')
        #create a "home" reference
        pd.DataFrame({"ID" : ped.loc[ped.cat.isin(["k", "pb"]) & (ped.Indiv.isin(popSplit.ID[popSplit.Group == "home"])), 'Indiv']}).\
            to_csv('IndForGenohome.txt', index=None, header=None, sep='\n')
        #create a "import" reference
        pd.DataFrame({"ID" : ped.loc[ped.cat.isin(["k", "pb"]) & (ped.Indiv.isin(popSplit.ID[popSplit.Group == "import"])), 'Indiv']}).\
            to_csv('IndForGenoimport.txt', index=None, header=None, sep='\n')


        # prestavi se v AlphaSim Dir
        #if not os.path.isfile(AlphaSimDir + 'ReferenceSize.txt') and os.path.isfile(AlphaSimDir + "IndForGeno.txt"):
        os.system("less IndForGeno.txt | wc -l > ReferenceSize.txt")
        #if not os.path.isfile(AlphaSimDir + 'ReferenceSize_home.txt') and os.path.isfile(AlphaSimDir + "IndForGenohome.txt"):
        print("Creating new ReferenceSizeHome file " + str(roundNo))
        os.system("less IndForGenohome.txt | wc -l > ReferenceSizehome.txt")
        #if not os.path.isfile(AlphaSimDir + 'ReferenceSize_import.txt') and os.path.isfile(AlphaSimDir + "IndForGenoimport.txt"):
        print("Creating new ReferenceSizeImport file " + str(roundNo))
        os.system("less IndForGenoimport.txt | wc -l > ReferenceSizeimport.txt")

    # Štartaj že po 20 gen kalsične selekcije
    AccHome = accuracies(AlphaSimDir, group = 'home', trait = traitHome)
    AccImport = accuracies(AlphaSimDir, group = 'import', trait = traitImport)
    #GenTrends = TBVCat(AlphaSimDir)
    # izvedi selekcijo, doloci kategorije zivali, dodaj novo generacijo in dodeli starse
    # pedigre se zapise v AlphaSimDir/SelectionFolder/ExternalPedigree.txt

    #tukaj razdeli populacijo na domačo in tujo
    #THis is now no longer needed, since prepareSelPed_groups does this
    #splitGenPed("PopulationSplit.txt")
    #tukaj izvedi celotno selekcijo v tuji populaciji --> naknadno shrani še očete z izberi_ocete_PT
    #v domači odberi in nastavi matere --> očete (za bm) uvoziš
    pedI, cI, sI, aI = selekcija_total('GenPed_EBVimport.txt', externalPedName="ExternalPedigreeimport", group=True, groupNumber=1,
                                           groupName = "import", noGroups=2, **selParimport)
    if selParimport['EBV']:
        Oce_import = pedI.izberi_ocete_PT(selParimport["pbUp"]) #tukaj so PT testirani očetje
    if selParimport['gEBV']:
        Oce_import = pedI.izberi_ocete_gen(selParimport["pbUp"])  # tukaj so genomsko testirani očetje

    #odberi starše domače populacije
    pedH, cH, sH, aH = selekcija_importOcetov('GenPed_EBVhome.txt', externalPedName="ExternalPedigreehome", group=True, groupNumber=0,
                                              groupName = "home", noGroups=2,
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
    SpecFile = AlphaSimSpec(os.getcwd(),WorkingDir + "/CodeDir", type="Multitrait" )  # AlphaSimSpec je class iz selection, ki omogoča nastavljanje parametrov AlphaSimSpec fila

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
    os.system("/exports/cmvm/eddie/eb/groups/tier2_hickey_external/R-3.4.2/bin/Rscript MeanHetMarker_Neutral_QTN_import.R " +
              str(roundNo+20) + " " + str(rep) + " " + str(scenarioHome) + " " + str(scenarioImport) + " " +
              str(traitHome) + " " + str(traitImport) + " " + str(strategy))
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
    #estimate EBVs for home population with only domestic data
    blupNextGen.computeEBV(group="home", dataGroup=True, prepareSelPed=False, traitEBV=traitHome,
                       multipleTraitsTBV=[traitHome, traitImport])
    #estimate EBVs for import population, use all data, create GenPed_EBVs.txt for both groups (only once, sinve it is one populationsplit file)
    blupNextGen.computeEBV(group="import", dataGroup=True, prepareSelPed=True, traitEBV=traitImport,
                       multipleTraitsTBV=[traitHome, traitImport])
    AccHome.saveAcc()
    AccImport.saveAcc()
    #GenTrends.saveTrends()
    #zdaj za vsako zapiši, ker vsakič na novo prebereš
    AccHome.writeAcc()
    AccImport.writeAcc()
    #GenTrends.writeTrends()


#os.system('rm -rf Chromosomes Selection && cp * ' + scenarioHome + scenarioImport +  str(rep))
#os.system('rm SimulatedData/UnrestrictedQtnIndivGenotypes.txt')
#os.system('rm SimulatedData/RestrictedQtnIndivGenotypes.txt')

