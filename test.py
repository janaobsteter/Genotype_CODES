# -*- coding: utf-8 -*-
from __future__ import division
import sys
import os
from PyQt4 import QtGui, QtCore, uic
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import math
import selection
from selection import selekcija_ena_gen, nastavi_cat, selekcija_total
from collections import defaultdict
import shutil
import pandas as  pd
import numpy as np

# nalozi GUI za selekcijo
qtCreatorFile = '/home/jana/Genotype_CODES/SelectionParameters.ui'  # Enter file here.
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)


# SelParam je class za okno
class SelParam(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)
        self.AlphaSimDir.clicked.connect(self.choose_dir)
        self.DoMagic.clicked.connect(self.selekcija)
        # self.DoMagic.clicked.connect(self.setSelParam)
        # self.DoMagic.clicked.connect(self.setText)


    # poveži tipko za AlphaSimDir z QFileDialog (posploši za vse take tipke!)
    def choose_dir(self):
        AlphaSimDirPath = QtGui.QFileDialog.getExistingDirectory(self, 'Save Directory')
        if AlphaSimDirPath:
            self.AlphaSimDirShow.setText(AlphaSimDirPath)

    # funkcija, ki ustvari dict in vseh vnešenih parametrov okna SelParam
    # ta dict potem daš funkciji, ki dela selekcijo oz. nastavlja kategorije
    def setSelParam(self):
        setParamDict = defaultdict()
        setParamDict['stNBn'] = int(self.NoNB.text()) if not self.NoNB.text().isEmpty() else 0
        setParamDict['kraveUp'] = int(self.kraveUpE.text()) if not self.kraveUpE.text().isEmpty() else 0
        setParamDict['bmOdbira'] = int(self.bmOdbiraE.text()) if not self.bmOdbiraE.text().isEmpty() else 0
        setParamDict['bmUp'] = int(self.bmUpE.text()) if not self.bmUpE.text().isEmpty() else 0
        setParamDict['cak'] = int(self.cakE.text()) if not self.cakE.text().isEmpty() else 0
        setParamDict['pripustUp'] = (self.pripustUpE.text()) if not self.pripustUpE.text().isEmpty() else 0
        setParamDict['pbUp'] = int(self.pbUpE.text()) if not self.pbUpE.text().isEmpty() else 0

        setParamDict['pripustDoz'] = int(self.pripustDozE.text()) if not self.pripustDozE.text().isEmpty() else 0
        setParamDict['mladiDoz'] = int(self.mladiDozE.text()) if not self.mladiDozE.text().isEmpty() else 0
        setParamDict['pozitivnoTestDoz'] = int(
            self.pozitivnoTestDozE.text()) if not self.pozitivnoTestDozE.text().isEmpty() else 0
        setParamDict['BurnInYN'] = self.BurnInYNE.isChecked()
        setParamDict['StBurnInGen'] = int(self.StBurnInGenE.text()) if not self.StBurnInGenE.text().isEmpty() else 0
        setParamDict['StSelGen'] = int(self.SelToGen.text()) if not self.SelToGen.text().isEmpty() else 0
        setParamDict['NumberOfSires'] = int(self.NoSires.text()) if not self.NoSires.text().isEmpty() else 0
        setParamDict['NumberOfDams'] = int(self.NoDams.text()) if not self.NoDams.text().isEmpty() else 0
        # setParamDict['AlphaSimDir'] = choose_dir()

        setParamDict['nrFn'] = int(setParamDict['stNBn'] * 0.5)
        setParamDict['telFn'] = int(self.telF.text() * setParamDict['nrFn']) if not self.telF.text().isEmpty() else 0
        setParamDict['ptn'] = int(self.pt.text() * setParamDict['telFn']) if not self.pt.text().isEmpty() else 0
        setParamDict['bmn'] = int(
            self.bm.text() * setParamDict['ptn'] * setParamDict['kraveUp']) if not self.bm.text().isEmpty() else 0

        setParamDict['nrMn'] = int(setParamDict['stNBn'] * 0.5)
        setParamDict['potomciNPn'] = int(
            self.potomciNP.text() * setParamDict['nrMn']) if not self.potomciNP.text().isEmpty() else 0
        setParamDict['vhlevljenin'] = int(
            self.vhlevljeni.text() * setParamDict['potomciNPn']) if not self.vhlevljeni.text().isEmpty() else 0
        setParamDict['bik12n'] = int(self.bik12.text() * setParamDict['nrMn']) if not self.bik12.text().isEmpty() else 0
        setParamDict['mladin'] = int(
            self.mladi.text() * setParamDict['vhlevljenin']) if not self.mladi.text().isEmpty() else 0
        setParamDict['pbn'] = int(self.pb.text() * setParamDict['mladin']) if not self.pb.text().isEmpty() else 0
        setParamDict['pripust1n'] = int(
            self.pripust1.text() * setParamDict['vhlevljenin']) if not self.pripust1.text().isEmpty() else 0
        setParamDict['pripust2n'] = int(round(setParamDict['pripust1n'] * (setParamDict['pripustUp'] - 1)))
        pd.DataFrame.from_dict(setParamDict, orient='index').to_csv('/home/jana/SelectionParamTEST.csv')
        print setParamDict

    # funkcija selekcija
    def selekcija(self):
        self.BurnInYN = self.BurnInYNE.isChecked()  # ali izvedeš tudi BurnIn
        self.StBurnInGen = int(self.StBurnInGenE.text())  #
        self.StSelGen = int(self.NoSelGen.text())
        self.StartSelGen = int(self.SelFromGen.text())
        self.StopSelGen = int(self.SelToGen.text())
        self.NumberOfSires = int(self.NoSires.text())  # number of sires in the population
        self.NumberOfDams = int(self.NoDams.tex())
        self.AlphaSimDir = self.AlphaSimDirShow.text()
        self.AlphaSimPed = str(self.AlphaSimDir) + '/SimulatedData/PedigreeAndGeneticValues.txt'

        ##############################################################################
        # burnin        #ce NAJ se izvede tudi burn in - ce ga se nimas
        ##############################################################################
        if self.BurnInYN:
            # prestavi se v AlphaSim Dir
            os.chdir(self.AlphaSimDir)
            shutil.copy(AlphaSimSpec.genSpecFile,
                        self.AlphaSimDir)  # skopiraj generično ALphaSimSpec datoteko v AlphaSimDir
            SpecFile = AlphaSimSpec()  # AlphaSimSpec je class iz selection, ki omogoča nastavljanje parametrov AlphaSimSpec fila
            SpecFile.setPedType("Internal")  # pedigree je za burn in internal
            SpecFile.setNB(self.stNBn)  # stevilo novorojenih
            SpecFile.setBurnInGen(self.StBurnInGen)  # stevilo burnINGen
            SpecFile.setSelGen(self.NoSelGen)  # st selection gen
            SpecFile.setNoSires(self.NumberOfSires)
            SpecFile.setNoDams(self.NumberOfDams)
            SpecFile.turnOnGenFlex()
            SpecFile.setFlexGenToFrom(1, (self.StBurnInGen + 1))  # pozeni od generacije 1 do burnin+1
            SpecFile.turnOnSelFlex()
            SpecFile.setExtPedForGen(
                self.StBurnInGen + 1)  # za katero generacijo uvozi external pedigre  - ena po burn in
            SpecFile.setTBVComp(1)  # skomputiraj TBV
            # pozenes ALPHASIM
            os.system('./AlphaSim1.05')
        ##############################################################################
        # ko si ustvaril burn in - ali pa ga imaš od prej - imaš torej PEdigreeAndGeneticValues
        # za VSAK KROG selekcije
        # 1) vzami AlphaSim sproduciran PedigreeAndGeneticValues
        # 2) preračunaj TBV --> EBV, doloci zeljeno korelacijo
        # 3) nastavi kategorije (nastavi kategorije glede na EBV) ali izvedi selekcijo (določi kategorije glede na prejšenje leto)
        # obe določita starše novim živalim - torej dodata novo generacijo plus starše
        # obe ustvarita external pedigree za naslednjo generacijo

        for roundNo in range(self.StartSelGen, (self.StopSelGen + 1)):  # za vsak krog selekcije
            # precacunas EBV v Ru in zapises PEDIGRE --> naredi ti GenPed_EBV.txt v cwd
            self.computeEBV(0.8)

            # USTVARI EXTERNAL PEDIGREE
            # doloci kategorije zivalim v pedigreju - če je to prvi krog, nastavi kategorije,
            # ce pa je to eden od naslednjih krogov, pa preberi kategorije iz prejsnje generacije
            # selekcija_total zapise kategorije, sex in active za vsako generacijo
            # nastavi_cat in selekcija_total ti zapišeta ExternalPedigree.txt
            if roundNo == 1:  # če je to prvi krog - nimaš še kategorij od prej
                nastavi_cat('GenPed_EBV.txt', **self.setSelParam())

            elif roundNo > 1:
                # izvedi selekcijo, doloci kategorije zivali, dodaj novo generacijo in dodeli starse
                # pedigre se zapise v AlphaSimDir/ExternalPedigree.txt
                selekcija_total('GenPed_EBV.txt', **self.setSelParam())

            # prestavi se v AlphaSim Dir
            os.chdir(self.AlphaSimDir)
            # kopiraj pedigre v selection folder
            shutil.copy('ExternalPedigree.txt', self.AlphaSimDir + '/Selection/SelectionFolder' + str(roundNo) + '/')
            # TUKAJ POTEM popravis AlphaSimSpec
            # PRVIc PO BURN IN-U
            shutil.copy(AlphaSimSpec.genSpecFile,
                        self.AlphaSimDir)  # skopiraj generično ALphaSimSpec datoteko v AlphaSimDir
            SpecFile = AlphaSimSpec()  # AlphaSimSpec je class iz selection, ki omogoča nastavljanje parametrov AlphaSimSpec fila
            SpecFile.setPedType("ExternalPedigree.txt")
            SpecFile.setBurnInGen(self.StBurnInGen)
            SpecFile.setSelGen(self.StSelGen)
            SpecFile.setNoSires(0)
            SpecFile.setNoDams(0)
            SpecFile.turnOnGenFlex()
            SpecFile.setFlexGenToFrom((self.StBurnInGen + roundNo), (self.StBurnInGen + roundNo))
            SpecFile.turnOnSelFlex()
            SpecFile.setExtPedForGen(self.StBurnInGen + roundNo)
            SpecFile.setTBVComp(2)
            SpecFile.setNB(self.stNBn)
            # pozenes ALPHASIM
            os.system('./AlphaSim1.05')

    def setText(self):
        self.message.setText('Button Clicked')

app = QtGui.QApplication(sys.argv)
window = SelParam()

window.NoNB.setText('6700')
window.setSelParam()

