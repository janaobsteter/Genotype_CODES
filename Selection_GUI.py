import sys
from PyQt4 import QtGui, QtCore, uic
import math
import selection
from selection import selekcija_ena_gen, nastavi_cat, selekcija_total
from collections import defaultdict


qtCreatorFile = '/home/jana/Genotipi/Genotipi_CODES/SelectionParameters.ui' # Enter file here.
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)


class SelParam(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)
        self.AlphaSimDir.clicked.connect(self.choose_dir)
        #self.DoMagic.clicked.connect(self.selekcija())

    def choose_dir(self):
        AlphaSimDirPath = QtGui.QFileDialog.getExistingDirectory(self, 'Save Directory')
        if AlphaSimDirPath:
            self.AlphaSimDirShow.setText(AlphaSimDirPath)

    def setSelParam(self):
        self.setParamDict = defaultdict()
        setParamDict['stNBn'] = int(self.NoNB.text())
        setParamDict['nrFn'] = int(self.stNBn) * 0.5
        setParamDict['telFn'] = int(self.telF.text()) * self.nrFn
        setParamDict['ptn'] = int(self.pt.text()) * self.telFn
        setParamDict['bmn'] = int(self.bm.text()) * self.ptn * self.kraveUp

        setParamDict['nrMn'] = int(self.stNBn) * 0.5
        setParamDict['potomciNPn'] = int(self.potomciNP.text()) * self.nrMn
        setParamDict['vhlevljenin'] = int(self.vhlevljeni.text()) * self.potomciNPn
        setParamDict['mladin'] = int(self.mladi.text()) * self.vhlevljenin
        setParamDict['pbn'] = int(self.pb.text()) * self.mladin
        setParamDict['pripust1n'] = int(self.pripust1.text()) * self.vhlevljenin
        setParamDict['pripust2n'] = int(round(self.pripust1n * (self.pripustUp - 1)))

        setParamDict['kraveUp'] = int(self.kraveUpE.text())
        setParamDict['bmOdbira'] = int(self.bmOdbiraE.text())
        setParamDict['bmUp'] = int(self.bmUpE.text())
        setParamDict['cak'] = int(self.cakE.text())
        setParamDict['pripustUp'] = int(self.pripustUpE.text())
        setParamDict['pbUp'] = int(self.pbUpE.text())

        setParamDict['pripustDoz'] = int(self.pripustDozE.text())
        setParamDict['mladiDoz'] = int(self.mladiDozE.text())
        setParamDict['pozitivnoTestDoz'] = int(self.pozitivnoTestDozE.text())
        setParamDict['BurnInYN'] = self.BurnInYNE.isChecked()
        setParamDict['StBurnInGen'] = int(self.StBurnInGenE.text())
        setParamDict['StSelGen'] = int(self.SelToGen.text())
        setParamDict['NumberOfSires'] = int(self.NoSires.text())
        setParamDict['NumberOfDams'] = int(self.NoDams.tex())
        setParamDict['AlphaSimDir'] = self.choose_dir()
        return setParamDict

    def selekcija(self):
        self.BurnInYN = self.BurnInYNE.isChecked()
        self.StBurnInGen = int(self.StBurnInGenE.text())
        self.StSelGen = int(self.SelToGen.text())
        self.NumberOfSires = int(self.NoSires.text())
        self.NumberOfDams = int(self.NoDams.tex())
        self.AlphaSimDir = self.AlphaSimDirShow.text()
        self.AlphaSimPed = str(self.AlphaSimDir) + '/SimulatedData/PedigreeAndGeneticValues.txt'

        if self.BurnInYN:
            for roundNo in range(self.StSelGen + 1):
                if roundNo == 0:  # do burn in
                    # prestavi se v AlphaSim Dir
                    os.chdir(self.AlphaSimDir)
                    shutil.copy('/home/jana/Genotipi/Genotipi_CODES/AlphaSimSpec.txt', self.AlphaSimDir)
                    os.system('sed -i "s|PedigreeType|Internal|g" AlphaSimSpec.txt')
                    os.system('sed -i "s|EnterBurnInGenerationNumber|' + str(self.StBurnInGen) + '|g" AlphaSimSpec.txt')
                    os.system('sed -i "s|EnterSelectionGenerationNumber|' + str(self.StSelGen) + '|g" AlphaSimSpec.txt')
                    os.system('sed -i "s|EnterNumberOfSires|' + str(self.NumberOfSires) + '|g" AlphaSimSpec.txt')
                    os.system('sed -i "s|EnterNumberOfDams|' + str(self.NumberOfDams) + '|g" AlphaSimSpec.txt')
                    os.system('sed -i "s|TurnOnGenFlex|On|g" AlphaSimSpec.txt')
                    os.system('sed -i "s|StartFlexGen,StopFlexGen|1,' + str(self.StBurnInGen + 1) + '|g" AlphaSimSpec.txt')
                    os.system('sed -i "s|TurnOnSelFlex|On|g" AlphaSimSpec.txt')
                    os.system('sed -i "s|TheImportedGenerationPed|' + str(self.StBurnInGen + 1) + '|g" AlphaSimSpec.txt')
                    os.system('sed -i "s|TBVComputation|1|g" AlphaSimSpec.txt')
                    os.system('sed -i "s|EnterIndividualInPopulation|' + str(self.stNBn) + '|g" AlphaSimSpec.txt')
                    # pozenes ALPHASIM
                    os.system('./AlphaSim1.05')

                elif roundNo == 1:
                    os.chdir('/home/jana/Genotipi/Genotipi_CODES/')
                    # preracunas EBV v Ru in ZAPIses PEDIGRE
                    shutil.copy("Rcorr_PedEBV.R", "Rcorr_PedEBV_ThisGen.R")
                    os.system('sed -i "s|AlphaSimPed|' + self.AlphaSimPed + '|g" Rcorr_PedEBV_ThisGen.R')
                    call('Rscript Rcorr_PedEBV_ThisGen.R', shell=True)
                    # tukaj nastvais zacetne kategorije
                    #global ped, categories, sex, active
                    ped, categories, sex, active = nastavi_cat('GenPed_EBV.txt',*self.setSelParam())
                    # prestavi se v AlphaSim Dir
                    os.chdir(self.AlphaSimDir)
                    # kopiraj pedigre v selection folder
                    shutil.copy('ExternalPedigree.txt', self.AlphaSimDir + '/Selection/SelectionFolder' + str(roundNo) + '/')
                    # TUKAJ POTEM POPRAVIs AlphaSimSpec
                    # PRVIc PO BURN IN-U
                    shutil.copy('/home/jana/Genotipi/Genotipi_CODES/AlphaSimSpec.txt', self.AlphaSimDir)
                    os.system('sed -i "s|PedigreeType|ExternalPedigree.txt|g" AlphaSimSpec.txt')
                    os.system('sed -i "s|EnterBurnInGenerationNumber|' + str(self.StBurnInGen) + '|g" AlphaSimSpec.txt')
                    os.system('sed -i "s|EnterSelectionGenerationNumber|' + str(self.StSelGen) + '|g" AlphaSimSpec.txt')
                    os.system('sed -i "s|EnterNumberOfSires|0|g" AlphaSimSpec.txt')
                    os.system('sed -i "s|EnterNumberOfDams|0|g" AlphaSimSpec.txt')
                    os.system('sed -i "s|TurnOnGenFlex|On|g" AlphaSimSpec.txt')
                    os.system('sed -i "s|StartFlexGen,StopFlexGen|' + str(self.StBurnInGen + roundNo) + ',' + \
                              str(self.StBurnInGen + roundNo) + '|g" AlphaSimSpec.txt')
                    os.system('sed -i "s|TurnOnSelFlex|On|g" AlphaSimSpec.txt')
                    os.system('sed -i "s|TheImportedGenerationPed|' + str(self.StBurnInGen + roundNo) + \
                              '|g" AlphaSimSpec.txt')
                    os.system('sed -i "s|TBVComputation|2|g" AlphaSimSpec.txt')
                    os.system('sed -i "s|EnterIndividualInPopulation|' + str(self.stNBn) + '|g" AlphaSimSpec.txt')

                    # Pozenes ALPHASIM
                    os.system('./AlphaSim1.05')


                elif roundNo > 1:
                    os.chdir('/home/jana/Genotipi/Genotipi_CODES/')
                    # prerauncas EBV v Ru in ZAPIses PEDIGRE
                    shutil.copy("Rcorr_PedEBV.R", "Rcorr_PedEBV_ThisGen.R")
                    os.system('sed -i "s|AlphaSimPed|' + self.AlphaSimPed + '|g" Rcorr_PedEBV_ThisGen.R')
                    call('Rscript Rcorr_PedEBV_ThisGen.R', shell=True)
                    # tukaj nastavis zacetne kategorije
                    #global ped, categories, sex, active
                    ped, categories, sex, active = selekcija_ena_gen('GenPed_EBV.txt', *self.setSelParam())

                    # prestavi se v AlphaSim Dir
                    os.chdir(AlphaSimDir)
                    # kopiraj pedigre v selection folder
                    shutil.copy('ExternalPedigree.txt', self.AlphaSimDir + '/Selection/SelectionFolder' + str(roundNo) + '/')
                    # TUKAJ POTEM POPRAVIs AlphaSimSpec
                    # PRVIc PO BURN IN-U
                    shutil.copy('/home/jana/Genotipi/Genotipi_CODES/AlphaSimSpec.txt', self.AlphaSimDir)
                    os.system('sed -i "s|PedigreeType|ExternalPedigree.txt|g" AlphaSimSpec.txt')
                    os.system('sed -i "s|EnterBurnInGenerationNumber|' + str(self.StBurnInGen) + '|g" AlphaSimSpec.txt')
                    os.system('sed -i "s|EnterSelectionGenerationNumber|' + str(self.StSelGen) + '|g" AlphaSimSpec.txt')
                    os.system('sed -i "s|EnterNumberOfSires|0|g" AlphaSimSpec.txt')
                    os.system('sed -i "s|EnterNumberOfDams|0|g" AlphaSimSpec.txt')
                    os.system('sed -i "s|TurnOnGenFlex|On|g" AlphaSimSpec.txt')
                    os.system('sed -i "s|StartFlexGen,StopFlexGen|' + str(self.StBurnInGen + roundNo) + ',' + str(
                        self.StBurnInGen + roundNo) + '|g" AlphaSimSpec.txt')
                    os.system('sed -i "s|TurnOnSelFlex|On|g" AlphaSimSpec.txt')
                    os.system(
                        'sed -i "s|TheImportedGenerationPed|' + str(self.StBurnInGen + roundNo) + '|g" AlphaSimSpec.txt')
                    os.system('sed -i "s|TBVComputation|2|g" AlphaSimSpec.txt')
                    os.system('sed -i "s|EnterIndividualInPopulation|' + str(self.stNBn) + '|g" AlphaSimSpec.txt')

                    # pozenes ALPHASIM
                    os.system('./AlphaSim1.05')

            if not BurnInYN:
                for roundNo in range(1, (self.StSelGen + 1)):
                    if roundNo == 1:
                        os.chdir('/home/jana/Genotipi/Genotipi_CODES/')
                        # precacunas EBV v Ru in zapises PEDIGRE
                        shutil.copy("Rcorr_PedEBV.R", "Rcorr_PedEBV_ThisGen.R")
                        os.system('sed -i "s|AlphaSimPed|' + self.AlphaSimPed + '|g" Rcorr_PedEBV_ThisGen.R')
                        call('Rscript Rcorr_PedEBV_ThisGen.R', shell=True)
                        # tukaj nastavis zacetne kategorije
                       # global ped, categories, sex, active
                        ped, categories, sex, active = nastavi_cat('GenPed_EBV.txt', *self.setSelParam())

                        # prestavi se v AlphaSim Dir
                        os.chdir(self.AlphaSimDir)
                        # kopiraj pedigre v selection folder
                        shutil.copy('ExternalPedigree.txt',
                                    self.AlphaSimDir + '/Selection/SelectionFolder' + str(roundNo) + '/')
                        # TUKAJ POTEM popravis AlphaSimSpec
                        # PRVIc PO BURN IN-U
                        shutil.copy('/home/jana/Genotipi/Genotipi_CODES/AlphaSimSpec.txt', self.AlphaSimDir)
                        os.system('sed -i "s|PedigreeType|ExternalPedigree.txt|g" AlphaSimSpec.txt')
                        os.system('sed -i "s|EnterBurnInGenerationNumber|' + str(self.StBurnInGen) + '|g" AlphaSimSpec.txt')
                        os.system('sed -i "s|EnterSelectionGenerationNumber|' + str(self.StSelGen) + '|g" AlphaSimSpec.txt')
                        os.system('sed -i "s|EnterNumberOfSires|0|g" AlphaSimSpec.txt')
                        os.system('sed -i "s|EnterNumberOfDams|0|g" AlphaSimSpec.txt')
                        os.system('sed -i "s|TurnOnGenFlex|On|g" AlphaSimSpec.txt')
                        os.system('sed -i "s|StartFlexGen,StopFlexGen|' + str(self.StBurnInGen + roundNo) + ',' + str(
                            self.StBurnInGen + roundNo) + '|g" AlphaSimSpec.txt')
                        os.system('sed -i "s|TurnOnSelFlex|On|g" AlphaSimSpec.txt')
                        os.system(
                            'sed -i "s|TheImportedGenerationPed|' + str(self.StBurnInGen + roundNo) + '|g" AlphaSimSpec.txt')
                        os.system('sed -i "s|TBVComputation|2|g" AlphaSimSpec.txt')
                        os.system('sed -i "s|EnterIndividualInPopulation|' + str(self.stNBn) + '|g" AlphaSimSpec.txt')

                        # pozenes ALPHASIM
                        os.system('./AlphaSim1.05')

                    elif roundNo > 1:
                        os.chdir('/home/jana/Genotipi/Genotipi_CODES/')

                        shutil.copy("Rcorr_PedEBV.R", "Rcorr_PedEBV_ThisGen.R")
                        os.system('sed -i "s|AlphaSimPed|' + self.AlphaSimPed + '|g" Rcorr_PedEBV_ThisGen.R')
                        call('Rscript Rcorr_PedEBV_ThisGen.R', shell=True)
                        # global ped, categories, sex, active
                        ped, categories, sex, active = selekcija_ena_gen('GenPed_EBV.txt',*self.setSelParam())
                        # prestavi se v AlphaSim Dir
                        os.chdir(self.AlphaSimDir)
                        # kopiraj pedigre v selection folder
                        os.system('mkdir ' + self.AlphaSimDir + '/Selection/SelectionFolder' + str(roundNo))
                        shutil.copy('ExternalPedigree.txt',
                                    self.AlphaSimDir + '/Selection/SelectionFolder' + str(roundNo) + '/')

                        shutil.copy('/home/jana/Genotipi/Genotipi_CODES/AlphaSimSpec.txt', self.AlphaSimDir)
                        os.system('sed -i "s|PedigreeType|ExternalPedigree.txt|g" AlphaSimSpec.txt')
                        os.system('sed -i "s|EnterBurnInGenerationNumber|' + str(self.StBurnInGen) + '|g" AlphaSimSpec.txt')
                        os.system('sed -i "s|EnterSelectionGenerationNumber|' + str(self.StSelGen) + '|g" AlphaSimSpec.txt')
                        os.system('sed -i "s|EnterNumberOfSires|0|g" AlphaSimSpec.txt')
                        os.system('sed -i "s|EnterNumberOfDams|0|g" AlphaSimSpec.txt')
                        os.system('sed -i "s|TurnOnGenFlex|On|g" AlphaSimSpec.txt')
                        os.system('sed -i "s|StartFlexGen,StopFlexGen|' + str(self.StBurnInGen + roundNo) + ',' + str(
                            self.StBurnInGen + roundNo) + '|g" AlphaSimSpec.txt')
                        os.system('sed -i "s|TurnOnSelFlex|On|g" AlphaSimSpec.txt')
                        os.system(
                            'sed -i "s|TheImportedGenerationPed|' + str(self.StBurnInGen + roundNo) + '|g" AlphaSimSpec.txt')
                        os.system('sed -i "s|TBVComputation|2|g" AlphaSimSpec.txt')
                        os.system('sed -i "s|EnterIndividualInPopulation|' + str(self.stNBn) + '|g" AlphaSimSpec.txt')

                        # pozenes ALPHASIM
                        os.system('./AlphaSim1.05')




if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = SelParam()
    window.show()
    sys.exit(app.exec_())

    selParam = window.setSelParam()


