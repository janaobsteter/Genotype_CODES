# -*- coding: utf-8 -*-
from __future__ import division
import sys
import os
from PyQt4 import QtGui, QtCore, uic
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import math
import selection
from selection import *
from selection import nastavi_cat, selekcija_total
from collections import defaultdict
import shutil
import pandas as pd
import numpy as np



#nalozi GUI za selekcijo
#qtCreatorFile = '/home/jana/Genotipi/Genotipi_CODES/SelectionParameters.ui' # Enter file here.
qtCreatorFile = '/home/jana/Genotype_CODES/SelectionParameters.ui' # Enter file here.
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)



#SelParam je class za okno
class SelParam(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)
        self.AlphaSimDir.clicked.connect(self.choose_dir)
        self.DoMagic.clicked.connect(self.selekcija)
        self.DoMagic.clicked.connect(self.setSelParam)
        self.SpecFile = AlphaSimSpec()
        # AlphaSimSpec je class iz selection, ki omogoča nastavljanje parametrov AlphaSimSpec fila
        self.setParamDict = defaultdict()
        # self.DoMagic.clicked.connect(self.setSelParam)
        #self.DoMagic.clicked.connect(self.setText)

        self.gEBV_YN.stateChanged.connect(self.disableEBV)
        self.EBV_YN.stateChanged.connect(self.disablegEBV)



    #poveži tipko za AlphaSimDir z QFileDialog (posploši za vse take tipke!)
    def choose_dir(self):
        AlphaSimDirPath = QtGui.QFileDialog.getExistingDirectory(self, 'Save Directory')
        if AlphaSimDirPath:
            self.AlphaSimDirShow.setText(AlphaSimDirPath)

    def disableEBV(self):
        if self.gEBV_YN.isChecked():
            self.vhlevljeni.setEnabled(False)
            self.mladi.setEnabled(False)
            self.cakE.setEnabled(False)
            self.mladiDozE.setEnabled(False)
            self.EBV_YN.setEnabled(False)
        if not self.gEBV_YN.isChecked():
            self.vhlevljeni.setEnabled(True)
            self.mladi.setEnabled(True)
            self.cakE.setEnabled(True)
            self.mladiDozE.setEnabled(True)
            self.EBV_YN.setEnabled(True)

    def disablegEBV(self):
        if self.EBV_YN.isChecked():
            self.gEBV_YN.setEnabled(False)
            #self.genmladi.setEnabled(False)
        if not self.EBV_YN.isChecked():
            self.gEBV_YN.setEnabled(True)
            #self.genmladi.setEnabled(True)


    #funkcija, ki ustvari dict in vseh vnešenih parametrov okna SelParam
    #ta dict potem daš funkciji, ki dela selekcijo oz. nastavlja kategorije
    def setSelParam(self):
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
        self.setParamDict['ptn'] = int(float(self.pt.text()) * self.setParamDict['telFn']) if not self.pt.text().isEmpty() else 0
        self.setParamDict['bmn'] = int(
            float(self.bm.text()) * self.setParamDict['ptn'] * self.setParamDict['kraveUp']) if not self.bm.text().isEmpty() else 0


        self.setParamDict['telMn'] = int(float(self.telM.text()) * self.setParamDict['nrMn']) if not self.telF.text().isEmpty() else 0
        self.setParamDict['potomciNPn'] = int(
            float(self.potomciNP.text()) * self.setParamDict['nrMn']) if not self.potomciNP.text().isEmpty() else 0

        if self.EBV_YN.isChecked():
            self.setParamDict['vhlevljenin'] = int(
                float(self.vhlevljeni.text()) * self.setParamDict['potomciNPn']) if not self.vhlevljeni.text().isEmpty() else 0
            self.setParamDict['mladin'] = int(
                float(self.mladi.text()) * self.setParamDict['vhlevljenin']) if not self.mladi.text().isEmpty() else 0
            self.setParamDict['pbn'] = int(float(self.pb.text()) * self.setParamDict['mladin']) if not self.pb.text().isEmpty() else 0
            self.setParamDict['pripust1n'] = int(round(
                float(self.pripust1.text()) * self.setParamDict['vhlevljenin'])) if not self.pripust1.text().isEmpty() else 0

        if self.gEBV_YN.isChecked():
            self.setParamDict['pbn'] = int(float(self.pb.text()) * self.setParamDict['potomciNPn']) if not self.pb.text().isEmpty() else 0
            self.setParamDict['pripust1n'] = int(round(
                float(self.pripust1.text()) * self.setParamDict['potomciNPn'])) if not self.pripust1.text().isEmpty() else 0

        self.setParamDict['bik12n'] = int(float(self.bik12.text()) * self.setParamDict['telMn']) if not self.bik12.text().isEmpty() else 0
        self.setParamDict['pripust2n'] = int(
            round(self.setParamDict['pripust1n'] * (float(self.setParamDict['pripustUp']) - 1)))
        pd.DataFrame.from_dict(self.setParamDict, orient='index').to_csv('/home/jana/SelectionParamTEST.csv', header=False)
        return self.setParamDict

    #funkcija selekcija
    def selekcija(self):
        self.BurnInYN = self.BurnInYNE.isChecked() #ali izvedeš tudi BurnIn
        self.StNB = int(self.NoNB.text())
        self.StBurnInGen = int(self.StBurnInGenE.text())
        self.StSelGen = int(self.NoSelGen.text())
        self.StartSelGen = int(self.SelFromGen.text())
        self.StopSelGen = int(self.SelToGen.text())
        self.NumberOfSires = int(self.NoSires.text()) if not self.NoSires.text().isEmpty() else 0#number of sires in the population
        self.NumberOfDams = int(self.NoDams.text()) if not self.NoDams.text().isEmpty() else 0 # to je za burn in - NoDams in NoSires
        self.AlphaSimDir = str(self.AlphaSimDirShow.text())
        self.AlphaSimPed = str(self.AlphaSimDir).strip('/.') + '/SimulatedData/PedigreeAndGeneticValues.txt'
        self.accuracyEBV = float(self.accEBV.text())


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
            os.system('./AlphaSim1.05')
##############################################################################
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
            # precacunas EBV v Ru in zapises PEDIGRE --> naredi ti GenPed_EBV.txt v cwd
            GenPedEBV = OrigPed(self.AlphaSimDir)
            GenPedEBV.computeEBV(self.accuracyEBV)


            #USTVARI EXTERNAL PEDIGREE
            #doloci kategorije zivalim v pedigreju - če je to prvi krog, nastavi kategorije,
            #ce pa je to eden od naslednjih krogov, pa preberi kategorije iz prejsnje generacije
            #selekcija_total zapise kategorije, sex in active za vsako generacijo
            #nastavi_cat in selekcija_total ti zapišeta ExternalPedigree.txt


            if roundNo == 1: #če je to prvi krog - nimaš še kategorij od prej
                nastavi_cat('GenPed_EBV.txt', **self.setSelParam())
    
            elif roundNo > 1:
                #izvedi selekcijo, doloci kategorije zivali, dodaj novo generacijo in dodeli starse
                #pedigre se zapise v AlphaSimDir/SelectionFolder/ExternalPedigree.txt
                selekcija_total('GenPed_EBV.txt', **self.setSelParam())
    

            #kopiraj pedigre v selection folder
            if not os.path.exists(self.AlphaSimDir + '/Selection/SelectionFolder' + str(roundNo) + '/'):
                os.makedirs(self.AlphaSimDir + '/Selection/SelectionFolder' + str(roundNo) + '/')
            shutil.copy(self.AlphaSimDir + 'ExternalPedigree.txt', self.AlphaSimDir + '/Selection/SelectionFolder' + str(roundNo) + '/')
            # TUKAJ POTEM popravis AlphaSimSpec
            # PRVIc PO BURN IN-U
            shutil.copy(self.SpecFile.genSpecFile, self.AlphaSimDir) #skopiraj generično ALphaSimSpec datoteko v AlphaSimDir
            self.SpecFile = AlphaSimSpec() #AlphaSimSpec je class iz selection, ki omogoča nastavljanje parametrov AlphaSimSpec fila
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
            #pozenes ALPHASIM
            os.system('./AlphaSim1.05')

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
    window.potomciNP.setText('0.0135')
    window.vhlevljeni.setText('0.6')
    window.mladi.setText('0.3')
    window.pb.setText('0.5')
    window.pripust1.setText('0.7')
    window.telM.setText('0.73')
    window.bik12.setText('0.12')
    window.accEBV.setText('0.8')
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
"""


if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = SelParam()
    window.show()
    window.NoNB.setText('100')
    window.telF.setText('0.96')
    window.pt.setText('0.9')
    window.bm.setText('0.1')
    window.potomciNP.setText('0.12')
    window.vhlevljeni.setText('0.667')
    window.mladi.setText('0.5')
    window.pb.setText('0.5')
    window.pripust1.setText('0.5')
    window.telM.setText('0.68')
    window.bik12.setText('0.333')
    window.accEBV.setText('0.8')
    window.kraveUpE.setText('4')
    window.bmOdbiraE.setText('2')
    window.bmUpE.setText('3')
    window.cakE.setText('3')
    window.pripustUpE.setText('1.5')
    window.pbUpE.setText('5')
    window.mladiDozE.setText('250')
    window.pripustDozE.setText('27')
    window.pozitivnoTestDozE.setText('220')
    window.StBurnInGenE.setText('10')
    window.NoSelGen.setText('20')
    window.NoDams.setText('40')
    window.NoSires.setText('1')
    window.SelFromGen.setText('1')
    window.SelToGen.setText('5')
    window.AlphaSimDirShow.setText('/home/jana/bin/AlphaSim1.05Linux/')
    sys.exit(app.exec_())
