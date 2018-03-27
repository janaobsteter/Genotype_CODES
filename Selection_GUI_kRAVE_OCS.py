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



#nalozi GUI za selekcijo
#qtCreatorFile = '/home/jana/Genotipi/Genotipi_CODES/SelectionParameters.ui' # Enter file here.
qtCreatorFile = '/home/jana/Genotipi/Genotipi_CODES/SelectionParameters.ui' # Enter file here.
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)


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
        #skopiraj paramFile za renumf90
        if self.sel == 'gen':
            shutil.copy(blupFiles.blupgenParamFile, blupFiles.AlphaSimDir + 'renumf90.par')  # skopiraj template blupparam file
        if self.sel == 'class':
            shutil.copy(blupFiles.blupgenParamFile_Clas, blupFiles.AlphaSimDir + 'renumf90.par')  # skopiraj template blupparam file

        # uredi blupparam file
        # get variance components from AlphaSim Output Files
        OutputFiles = AlphaSim_OutputFile(self.AlphaSimDir)
        genvar = OutputFiles.getAddVar()  # dobi additivno varianco
        resvar = OutputFiles.getResVar()  # dobi varianco za ostanek

        blupFiles.prepareParamFiles(genvar, resvar, self.AlphaSimDir + '/renumf90.par')  # set levels of random aniaml effect, add var and res var
        # the paramfile is now set
        blupFiles.makePed_gen()  # make ped file for blup, no Code!
        if self.sel == 'gen':
            GenFiles = snpFiles(self.AlphaSimDir)
            GenFiles.createBlupf90SNPFile()

        os.system("./renumf90 < renumParam")  # run renumf90

        #if self.sel == 'class': #ZDJ TEGA NI, KER JE POSEBEJ FILE!
        #if self.sel == 'class': #ZDJ TEGA NI, KER JE POSEBEJ FILE!
            #os.system("head -n-3 renf90.par > tmp && mv tmp renf90.par")
            #os.system("./blupf90 blupf90_Selection")  # run blupf90

        resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
        #os.system('./preGSf90 renf90.par')
        os.system('./blupf90 renf90.par')
        # os.system('./postGSf90 renf90.par')


        #renumber the solutions
        # copy the solution in a file that does not get overwritten
        os.system("bash Match_AFTERRenum.sh")
        shutil.copy('renumbered_Solutions', 'renumbered_Solutions_' + str(blupFiles.gen))
        #shutil.copy('solutions', 'renumbered_Solutions_' + str(blupFiles.gen))

        blupFiles.prepareSelPed()  # obtain solution and add them to
        # AlphaPed PedigreeAndGeneticValues files --> Write them to GenPed_EBV.txt, which is read by module selection


#SelParam je class za okno
class SelParam(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)
        self.AlphaSimDir.clicked.connect(self.choose_dir)
        self.DoMagic.clicked.connect(self.selekcija)
        self.DoMagic.clicked.connect(self.setSelParam)
        self.SpecFile = AlphaSimSpec('/home/jana/bin/AlphaSim1.05Linux/', '/home/jana/Genotipi/Genotipi_CODES/')
        # AlphaSimSpec je class iz selection, ki omogoča nastavljanje parametrov AlphaSimSpec fila
        self.setParamDict = defaultdict()
        #self.DoMagic.clicked.connect(self.setSelParam)
        #self.DoMagic.clicked.connect(self.setText)

        self.gEBV_YN.stateChanged.connect(self.disableEBV)
        self.EBV_YN.stateChanged.connect(self.disablegEBV)
        self.PA_YN.stateChanged.connect(self.disablegEBVEBV)
        self.genTest_mladi.stateChanged.connect(self.disable_genTest_mladi)
        self.genTest_gpb.stateChanged.connect(self.disable_genTest_gpb)

        criteria = ['random', 'PA', 'EBV']
        for i in [self.potomciNP_M_genC, self.potomciNP_F_genC, self.vhlevljeni_genC, self.mladi_genC, self.cak_genC, self.nrM_genC, self.nrF_genC, self.telF_genC, self.k_genC, self.bm_genC, self.pb_genC]:
            i.addItems(criteria)


    #poveži tipko za AlphaSimDir z QFileDialog (posploši za vse take tipke!)
    def choose_dir(self):
        AlphaSimDirPath = QtGui.QFileDialog.getExistingDirectory(self, 'Save Directory')
        if AlphaSimDirPath:
            self.AlphaSimDirShow.setText(AlphaSimDirPath)

    def disable_genTest_mladi(self):
        if self.genTest_mladi.isChecked():
            self.genTest_gpb.setEnabled(False)
            self.gpb_pb.setEnabled(False)
        if not self.genTest_mladi.isChecked():
            self.genTest_gpb.setEnabled(True)
            self.gpb_pb.setEnabled(True)




    def disable_genTest_gpb(self):
        if self.genTest_gpb.isChecked():
            self.genTest_mladi.setEnabled(False)
        if not self.genTest_gpb.isChecked():
            self.genTest_mladi.setEnabled(True)

    def disableEBV(self):
        if self.gEBV_YN.isChecked():
            #self.vhlevljeni.setEnabled(False)
            #self.mladi.setEnabled(False)
            #self.cakE.setEnabled(False)
            #self.mladiDozE.setEnabled(False)
            self.EBV_YN.setEnabled(False)
            self.PA_YN.setEnabled(False)
        if not self.gEBV_YN.isChecked():
            self.vhlevljeni.setEnabled(True)
            self.mladi.setEnabled(True)
            self.cakE.setEnabled(True)
            self.mladiDozE.setEnabled(True)
            self.EBV_YN.setEnabled(True)
            self.PA_YN.setEnabled(True)

    def disablegEBV(self):
        if self.EBV_YN.isChecked():
            self.gEBV_YN.setEnabled(False)
            self.PA_YN.setEnabled(False)
            self.genpb.setEnabled(False)
            self.genpripust1.setEnabled(False)
            self.GenotypedInd.setEnabled(False)
        if not self.EBV_YN.isChecked():
            self.gEBV_YN.setEnabled(True)
            self.PA_YN.setEnabled(True)
            self.GenotypedInd.setEnabled(True)
            #self.genmladi.setEnabled(True)

    def disablegEBVEBV(self):
        if self.PA_YN.isChecked():
            self.gEBV_YN.setEnabled(False)
            self.EBV_YN.setEnabled(False)
            self.GenotypedInd.setEnabled(False)
        if not self.PA_YN.isChecked():
            self.gEBV_YN.setEnabled(True)
            self.EBV_YN.setEnabled(True)
            self.GenotypedInd.setEnabled(True)
            #self.genmladi.setEnabled(True)

    #funkcija, ki ustvari dict in vseh vnešenih parametrov okna SelParam
    #ta dict potem daš funkciji, ki dela selekcijo oz. nastavlja kategorije
    def setSelParam(self):
        self.setParamDict['genotyped'] = [(str(x.text()), float(xP.text()), str(xC.currentText()), str(sex)) for (x, xP, xC, sex) in
                                          [(self.potomciNP_M_gen, self.potomciNP_M_genP, self.potomciNP_M_genC, 'M'),
                                           (self.potomciNP_F_gen, self.potomciNP_F_genP, self.potomciNP_F_genC, 'F'),
                                           (self.nrF_gen, self.nrF_genP, self.nrF_genC, 'F'),
                                           (self.nrM_gen, self.nrM_genP, self.nrM_genC, 'M'),
                                           (self.telF_gen, self.telF_genP, self.telF_genC, 'F'),
                                           (self.k_gen, self.k_genP, self.k_genC, 'F'),
                                           (self.bm_gen, self.bm_genP, self.bm_genC, 'F'),
                                           (self.cak_gen, self.cak_genP, self.cak_genC, 'M'),
                                           (self.vhlevljeni_gen, self.vhlevljeni_genP, self.vhlevljeni_genC, 'M'),
                                           (self.mladi_gen, self.mladi_genP, self.mladi_genC, 'M'),
                                           (self.pb_gen, self.pb_genP, self.pb_genC, 'M')] if x.isChecked()]
        self.setParamDict['EBV'] = self.EBV_YN.isChecked()
        self.setParamDict['gEBV'] = self.gEBV_YN.isChecked()
        self.setParamDict['stNBn'] = int(self.NoNB.text()) if not self.NoNB.text().isEmpty() else 0
        self.setParamDict['kraveUp'] = int(self.kraveUpE.text()) if not self.kraveUpE.text().isEmpty() else 0
        self.setParamDict['bmOdbira'] = int(self.bmOdbiraE.text()) if not self.bmOdbiraE.text().isEmpty() else 0
        self.setParamDict['bmUp'] = int(self.bmUpE.text()) if not self.bmUpE.text().isEmpty() else 0
        self.setParamDict['cak'] = int(self.cakE.text()) if not self.cakE.text().isEmpty() else 0
        self.setParamDict['pripustUp'] = float((self.pripustUpE.text())) if not self.pripustUpE.text().isEmpty() else 0
        self.setParamDict['pbUp'] = int(self.pbUpE.text()) if not self.pbUpE.text().isEmpty() else 0

        self.setParamDict['pripustDoz'] = int(self.pripustDozE.text()) if not self.pripustDozE.text().isEmpty() else 0
        self.setParamDict['mladiDoz'] = int(self.mladiDozE.text()) if not self.mladiDozE.text().isEmpty() else 0
        self.setParamDict['pozitivnoTestDoz'] = int(
            self.pozitivnoTestDozE.text()) if not self.pozitivnoTestDozE.text().isEmpty() else 0
        self.setParamDict['BurnInYN'] = self.BurnInYNE.isChecked()
        self.setParamDict['StBurnInGen'] = int(self.StBurnInGenE.text()) if not self.StBurnInGenE.text().isEmpty() else 0
        self.setParamDict['StSelGen'] = int(self.SelToGen.text()) if not self.SelToGen.text().isEmpty() else 0
        self.setParamDict['NumberOfSires'] = int(self.NoSires.text()) if not self.NoSires.text().isEmpty() else 0
        self.setParamDict['NumberOfDams'] = int(self.NoDams.text()) if not self.NoDams.text().isEmpty() else 0

        # self.setParamDict['AlphaSimDir'] = choose_dir()

        self.setParamDict['nrMn'] = int(self.setParamDict['stNBn'] * 0.5)
        self.setParamDict['potomciNPn'] = int(
            float(self.potomciNP.text()) * self.setParamDict['nrMn']) if not self.potomciNP.text().isEmpty() else 0

        self.setParamDict['nrFn'] = int(self.setParamDict['stNBn'] * 0.5) - self.setParamDict['potomciNPn']
        self.setParamDict['telFn'] = int(float(self.telF.text()) * self.setParamDict['nrFn']) if not self.telF.text().isEmpty() else 0
        self.setParamDict['telFnTotal'] = int(float(self.telF.text()) * self.setParamDict['nrFn'] + self.setParamDict['potomciNPn']) \
            if not self.telF.text().isEmpty() else 0 #tu prištej, da dobiš končno število telFn (pridejo po dveh poteh, ene preko potomciNP)
        self.setParamDict['ptn'] = int(float(self.pt.text()) * self.setParamDict['telFnTotal']) if not self.pt.text().isEmpty() else 0
        self.setParamDict['MinusDamLact'] = int((1 - float(self.cowsLactLact.text())) * self.setParamDict['ptn']) if not self.cowsLactLact.text().isEmpty() else 0
        self.setParamDict['kn'] = int(self.setParamDict['ptn'] *
                                     self.setParamDict['kraveUp'] - (self.setParamDict['MinusDamLact'] * sum(range((self.setParamDict['kraveUp'])))))
        self.setParamDict['bmn'] = int(
            float(self.bm.text()) * (self.setParamDict['ptn'] *
                                     self.setParamDict['kraveUp'] - (self.setParamDict['MinusDamLact'] * sum(range((self.setParamDict['kraveUp'])))))) if not self.bm.text().isEmpty() else 0

        self.setParamDict['telMn'] = int(float(self.telM.text()) * self.setParamDict['nrMn']) if not self.telF.text().isEmpty() else 0
        self.setParamDict['potomciNPn'] = int(
            float(self.potomciNP.text()) * self.setParamDict['nrMn']) if not self.potomciNP.text().isEmpty() else 0

      #  if self.EBV_YN.isChecked():
        self.setParamDict['vhlevljenin'] = int(
            float(self.vhlevljeni.text()) * self.setParamDict['potomciNPn']) if not self.vhlevljeni.text().isEmpty() else 0
        self.setParamDict['mladin'] = int(
            float(self.mladi.text()) * self.setParamDict['vhlevljenin']) if not self.mladi.text().isEmpty() else 0
        self.setParamDict['pbn'] = int(float(self.pb.text()) * self.setParamDict['mladin']) if not self.pb.text().isEmpty() else 0
        self.setParamDict['pripust1n'] = int(round(
            float(self.pripust1.text()) * self.setParamDict['vhlevljenin'])) if not self.pripust1.text().isEmpty() else 0

        #if self.gEBV_YN.isChecked():
        self.setParamDict['genpbn'] = int(round(float(self.genpb.text()) * self.setParamDict['potomciNPn'])) if not self.pb.text().isEmpty() else 0
        self.setParamDict['genpripust1n'] = int(round(round(
            float(self.genpripust1.text()) * self.setParamDict['potomciNPn'], 1))) if not self.pripust1.text().isEmpty() else 0

        self.setParamDict['bik12n'] = int(float(self.bik12.text()) * self.setParamDict['telMn']) if not self.bik12.text().isEmpty() else 0
        self.setParamDict['pripust2n'] = float(float(self.setParamDict['pripustUp']) - 1)
        self.setParamDict['EliteDamsPTBulls'] = self.EliteDamsPTBulls.isChecked()
        self.setParamDict['EliteDamsGenBulls'] = self.EliteDamsGenBulls.isChecked()
        self.setParamDict['EliteDamsPABulls'] = self.EliteDamsPABulls.isChecked()
        self.setParamDict['CowsGenBulls_Per'] = int((float(self.CowsGenBullsPer.text()) / 100) * self.setParamDict['stNBn']) \
            if not self.CowsGenBullsPer.text().isEmpty() else 0
        self.setParamDict['genTest_mladi'] = self.genTest_mladi.isChecked()
        self.setParamDict['genTest_gpb'] = self.genTest_gpb.isChecked()
        self.setParamDict['gpb_pb'] = self.gpb_pb.isChecked()
        self.setParamDict['UpdateGenRef'] = self.UpdateGenRef.isChecked()
        self.setParamDict['sexToUpdate'] = str(self.sexToUpdate.text()) if not self.sexToUpdate.text().isEmpty() else ''
        self.setParamDict['NbUpdatedGen'] = int(self.NbUpdatedGen.text()) if not self.NbUpdatedGen.text().isEmpty() else 0
        self.setParamDict['AlphaSimDir'] = str(self.AlphaSimDirShow.text())
        pd.DataFrame.from_dict(self.setParamDict, orient='index').to_csv('/home/jana/SelectionParam.csv', header=False)
        pd.DataFrame.from_dict(self.setParamDict, orient='index').to_csv(self.setParamDict.get('AlphaSimDir') + '/SelectionParam.csv', header=False)
        return self.setParamDict

    #funkcija selekcija
    def selekcija(self):
        self.BurnInYN = self.BurnInYNE.isChecked() #ali izvedeš tudi BurnIn
        self.SelYN = self.PerformSel.isChecked() #ali izvedeš tudi BurnIn
        self.StNB = int(self.NoNB.text())
        self.StBurnInGen = int(self.StBurnInGenE.text())
        self.StSelGen = int(self.NoSelGen.text())
        self.StartSelGen = int(self.SelFromGen.text())
        self.StopSelGen = int(self.SelToGen.text())
        self.NumberOfSires = int(self.NoSires.text()) if not self.NoSires.text().isEmpty() else 0#number of sires in the population
        self.NumberOfDams = int(self.NoDams.text()) if not self.NoDams.text().isEmpty() else 0 # to je za burn in - NoDams in NoSires
        self.AlphaSimDir = str(self.AlphaSimDirShow.text())
        self.codeDir = '/home/jana/Genotipi/Genotipi_CODES/'
        self.AlphaSimPed = str(self.AlphaSimDir).strip('/.') + '/SimulatedData/PedigreeAndGeneticValues.txt'
        if self.EBV_YN.isChecked():
            self.seltype = 'class'
        if self.gEBV_YN.isChecked():
            self.seltype = 'gen'
        if self.PA_YN.isChecked():
            self.seltype = 'PA'


        ##############################################################################
#burnin        #ce NAJ se izvede tudi burn in - ce ga se nimas
##############################################################################
        if self.BurnInYN:
            # prestavi se v AlphaSim Dir
            os.chdir(self.AlphaSimDir)
            shutil.copy(self.SpecFile.genSpecFile, self.AlphaSimDir) #skopiraj generično ALphaSimSpec datoteko v AlphaSimDir
            self.SpecFile.setPedType("Internal") #pedigree je za burn in internal
            self.SpecFile.setNB(self.StNB) #stevilo novorojenih
            self.SpecFile.setBurnInGen(self.StBurnInGen) #stevilo burnINGen
            self.SpecFile.setSelGen(self.StSelGen) #st selection gen
            self.SpecFile.setNoSires(self.NumberOfSires)
            self.SpecFile.setNoDams(self.NumberOfDams)
            self.SpecFile.turnOnGenFlex()
            self.SpecFile.setFlexGenToFrom(1, (self.StBurnInGen + 1)) #pozeni od generacije 1 do burnin+1
            self.SpecFile.turnOnSelFlex()
            self.SpecFile.setExtPedForGen(self.StBurnInGen + 1) #za katero generacijo uvozi external pedigre  - ena po burn in
            self.SpecFile.setTBVComp(1) #skomputiraj TBV
            # pozenes ALPHASIM
            os.system('./AlphaSim1.08')
##############################################################################
        if self.SelYN: #če naj izvedem tudi selekcijo
            #ko si ustvaril burn in - ali pa ga imaš od prej - imaš torej PEdigreeAndGeneticValues
                #za VSAK KROG selekcije
            #1) vzami AlphaSim sproduciran PedigreeAndGeneticValues
            #2) preračunaj TBV --> EBV, doloci zeljeno korelacijo
            #3) nastavi kategorije (nastavi kategorije glede na EBV) ali izvedi selekcijo (določi kategorije glede na prejšenje leto)
            #obe določita starše novim živalim - torej dodata novo generacijo plus starše
            #obe ustvarita external pedigree za naslednjo generacijo


            for roundNo in range(self.StartSelGen, (self.StopSelGen + 1)): #za vsak krog selekcije
                # prestavi se v AlphaSim Dir
                os.chdir(self.AlphaSimDir)
                if not os.path.isfile(self.AlphaSimDir + 'ReferenceSize.txt') and os.path.isfile(self.AlphaSimDir + "IndForGeno.txt"):
                    os.system("less IndForGeno.txt | wc -l > ReferenceSize.txt")

                # USTVARI EXTERNAL PEDIGREE
                # doloci kategorije zivalim v pedigreju - če je to prvi krog, nastavi kategorije,
                # ce pa je to eden od naslednjih krogov, pa preberi kategorije iz prejsnje generacije
                # selekcija_total zapise kategorije, sex in active za vsako generacijo
                # nastavi_cat in selekcija_total ti zapišeta ExternalPedigree.txt
                if roundNo == 1:  # če je to prvi krog - nimaš še kategorij od prej, nimaš niti EBV-jev
                    # odstrani Blupf90 fajle iz prejšnjih runov - ker se merge-a
                    # enako tudi za generacijski interval in file z genotipi
                    if os.path.isfile(self.AlphaSimDir + 'Blupf90.dat'):
                        os.remove(self.AlphaSimDir + 'Blupf90.dat')
                    if os.path.isfile(self.AlphaSimDir + 'GenInts.txt'):
                        os.remove(self.AlphaSimDir + 'GenInts.txt')
                    if os.path.isfile(self.AlphaSimDir + 'GenoFile.txt'):
                        os.remove(self.AlphaSimDir + 'GenoFile.txt')
                    if os.path.isfile(self.AlphaSimDir + 'IndForGeno.txt'):
                        os.remove(self.AlphaSimDir + 'IndForGeno.txt')
                    if os.path.isfile(self.AlphaSimDir + 'GenTrends_gen.csv'):
                        os.remove(self.AlphaSimDir + 'GenTrends_gen.csv')
                    if os.path.isfile(self.AlphaSimDir + 'GenTrends_cat.csv'):
                        os.remove(self.AlphaSimDir + 'GenTrends_cat.csv')
                    if os.path.isfile(self.AlphaSimDir + 'Accuracies_Cat.csv'):
                        os.remove(self.AlphaSimDir + 'Accuracies_Cat.csv')
                    if os.path.isfile(self.AlphaSimDir + 'Accuracies_Gen.csv'):
                        os.remove(self.AlphaSimDir + 'Accuracies_Gen.csv')
                    #if os.path.isfile(self.AlphaSimDir + 'AccuraciesBV.csv'):
                        #os.remove(self.AlphaSimDir + 'AccuraciesBV.csv')

                    Acc = accuracies(self.AlphaSimDir) #nastavi
                    GenTrends = TBVCat(self.AlphaSimDir)
                    #nimaš GenPed_EBV.txt
                    blups = estimateBV(self.AlphaSimDir, self.codeDir, way='burnin_milk', sel=self.seltype)
                    blups.computeEBV(self.StNB) #tukaj izbriši samo fenotipe moških - ne morš po kategorijah, ker jih nimaš
                    #Acc.saveAcc()
                    nastavi_cat('GenPed_EBV.txt', **self.setSelParam())

                elif roundNo > 1:
                    Acc = accuracies(self.AlphaSimDir)
                    GenTrends = TBVCat(self.AlphaSimDir)
                    # izvedi selekcijo, doloci kategorije zivali, dodaj novo generacijo in dodeli starse
                    # pedigre se zapise v AlphaSimDir/SelectionFolder/ExternalPedigree.txt
                    selekcija_total('GenPed_EBV.txt', **self.setSelParam())

                # kopiraj pedigre v selection folder
                if not os.path.exists(self.AlphaSimDir + '/Selection/SelectionFolder' + str(roundNo) + '/'):
                    os.makedirs(self.AlphaSimDir + '/Selection/SelectionFolder' + str(roundNo) + '/')
                shutil.copy(self.AlphaSimDir + 'ExternalPedigree.txt',
                            self.AlphaSimDir + '/Selection/SelectionFolder' + str(roundNo) + '/')
                # TUKAJ POTEM popravis AlphaSimSpec
                # PRVIc PO BURN IN-U
                shutil.copy(self.SpecFile.genSpecFile,
                            self.AlphaSimDir)  # skopiraj generično ALphaSimSpec datoteko v AlphaSimDir
                self.SpecFile = AlphaSimSpec()  # AlphaSimSpec je class iz selection, ki omogoča nastavljanje parametrov AlphaSimSpec fila
                self.SpecFile.setPedType("ExternalPedigree.txt")
                self.SpecFile.setBurnInGen(self.StBurnInGen)
                self.SpecFile.setSelGen(self.StSelGen)
                self.SpecFile.setNoSires(0)
                self.SpecFile.setNoDams(0)
                self.SpecFile.turnOnGenFlex()
                self.SpecFile.setFlexGenToFrom((self.StBurnInGen + roundNo), (self.StBurnInGen + roundNo))
                self.SpecFile.turnOnSelFlex()
                self.SpecFile.setExtPedForGen(self.StBurnInGen + roundNo)
                self.SpecFile.setTBVComp(2)
                self.SpecFile.setNB(self.StNB)
                # pozenes ALPHASIM
                os.system('./AlphaSim1.08')

                # tukaj dodaj kategorije k PedigreeAndGeneticValues (AlphaSim File)
                PedCat = OrigPed(self.AlphaSimDir, '/home/jana/Genotipi/Genotipi_CODES')
                PedCat.addInfo() #to ti zapiše PedigreeAndGeneticValues_cat.txt v AlphaSim/SimualatedData

                #tukaj pridobi podatke za generacijske intervale
                GenInt = genInterval(self.AlphaSimDir) #tukaj preberi celoten pedigre
                if self.seltype == 'class':
                    GenInt.prepareGenInts(['vhlevljeni', 'pt']) #pri klasični so izrbrani potomci vhlevljeni (test in pripust) in plemenske telice
                if self.seltype == 'gen':
                    GenInt.prepareGenInts(['genTest', 'pt']) #pri klasični so izbrani potomci vsi genomsko testirani (pozTest in pripust) in plemenske telice
                blupNextGen = estimateBV(self.AlphaSimDir, way='milk', sel=self.seltype)
                blupNextGen.computeEBV(self.StNB) #estimate EBV with added phenotypes only of animals of certain category (here = milk)
                Acc.saveAcc()
                GenTrends.saveTrends()
                #zdaj za vsako zapiši, ker vsakič na novo prebereš
                Acc.writeAcc()
                GenTrends.writeTrends()
            print 'Process finished'

        def setText(self):
            self.message.setText('Button Clicked')




"""
if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = SelParam()
    window.show()
    window.NoNB.setText('6700')
    window.telF.setText('0.966')
    window.pt.setText('0.85')
    window.bm.setText('0.0127')
    window.cowsLactLact.setText('0.8')
    window.potomciNP.setText('0.0135')
    window.vhlevljeni.setText('0.6')
    window.mladi.setText('0.3')
    window.pb.setText('0.5')
    window.pripust1.setText('0.7')
    window.telM.setText('0.73')
    window.bik12.setText('0.12')
    window.kraveUpE.setText('4')
    window.bmOdbiraE.setText('2')
    window.bmUpE.setText('3')
    window.cakE.setText('3')
    window.pripustUpE.setText('1.4')
    window.pbUpE.setText('5')
    window.mladiDozE.setText('250')
    window.pripustDozE.setText('27')
    window.pozitivnoTestDozE.setText('220')
    window.StBurnInGenE.setText('10')
    window.NoSelGen.setText('20')
    window.SelFromGen.setText('2')
    window.SelToGen.setText('2')
    window.AlphaSimDirShow.setText('/home/jana/bin/AlphaSim1.05Linux/')
    sys.exit(app.exec_())

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = SelParam()
    window.show()
    window.NoNB.setText('100')
    window.telF.setText('0.96')
    window.cowsLactLact.setText('0.8')
    window.pt.setText('0.9')
    window.bm.setText('0.1')
    window.potomciNP.setText('0.12')
    window.vhlevljeni.setText('0.667')
    window.mladi.setText('0.5')
    window.pb.setText('0.5')
    window.genpb.setText('0.2')
    window.pripust1.setText('0.5')
    window.genpripust1.setText('0.8')
    window.telM.setText('0.68')
    window.bik12.setText('0.333')
    window.kraveUpE.setText('4')
    window.bmOdbiraE.setText('2')
    window.bmUpE.setText('3')
    window.cakE.setText('3')
    window.pripustUpE.setText('1.5')
    window.pbUpE.setText('5')
    window.mladiDozE.setText('25')
    window.pripustDozE.setText('5')
    window.pozitivnoTestDozE.setText('35')
    window.StBurnInGenE.setText('10')
    window.NoSelGen.setText('20')
    window.NoDams.setText('40')
    window.NoSires.setText('1')
    window.SelFromGen.setText('1')
    window.SelToGen.setText('5')
    window.AlphaSimDirShow.setText('/home/jana/bin/AlphaSim1.05Linux/')
    sys.exit(app.exec_())


"""
if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = SelParam()
    window.show()
    window.NoNB.setText('8640')
    window.NoSires.setText('4')
    window.NoDams.setText('3500')
    window.telF.setText('0.99')
    window.pt.setText('0.9')
    window.bm.setText('0.012')
    window.cowsLactLact.setText('0.8')
    window.potomciNP.setText('0.0105')
    window.vhlevljeni.setText('0.6')
    window.mladi.setText('0.3')
    window.pb.setText('0.5')
    window.genpb.setText('0.1')
    window.pripust1.setText('0.7')
    window.genpripust1.setText('0.9')
    window.telM.setText('0.73')
    window.bik12.setText('0.12')
    window.kraveUpE.setText('4')
    window.bmOdbiraE.setText('2')
    window.bmUpE.setText('3')
    window.cakE.setText('3')
    window.pripustUpE.setText('1.4')
    window.pbUpE.setText('5')
    window.mladiDozE.setText('250')
    window.pripustDozE.setText('27')
    window.pozitivnoTestDozE.setText('400')
    window.StBurnInGenE.setText('20')
    window.NoSelGen.setText('40')
    window.SelFromGen.setText('21')
    window.SelToGen.setText('40')
    window.potomciNP_M_genP.setText('100')
    window.k_genP.setText('10')
    window.pb_genP.setText('100')
    window.AlphaSimDirShow.setText('/home/jana/bin/AlphaSim1.05Linux/')
    sys.exit(app.exec_())
