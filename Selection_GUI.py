import sys
from PyQt4 import QtGui, QtCore, uic
import math
import selection
from selection import selekcija_ena_gen, nastavi_cat

qtCreatorFile = '/home/jana/Genotipi/Genotipi_CODES/SelectionParameters.ui' # Enter file here.
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)


class SelParam(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)
        self.AlphaSimDir.clicked.connect(self.choose_dir)
        self.AlphaSimDir = self.choose_dir()

    def choose_dir(self):
        return QtGui.QFileDialog.getExistingDirectory(self, 'Save Directory')

    def selekcija(self):
        #nastavi parametre
        self.stNBn = int(self.NoNB.text())
        self.nrFn = int(self.stNBn) * 0.5
        self.telFn = int(self.telF.text()) * self.nrFn
        self.ptn = int(self.pt.text()) * self.telFn
        self.bmn = int(self.bm.text()) * self.ptn * self.kraveUp

        self.nrMn = int(self.stNBn) * 0.5
        self.potomciNPn = int(self.potomciNP.text()) * self.nrMn
        self.vhlevljenin = int(self.vhlevljeni.text()) * self.potomciNPn
        self.mladin = int(self.mladi.text()) * self.vhlevljenin
        self.pbn = int(self.pb.text()) * self.mladin
        self.pripust1n = int(self.pripust1.text()) * self.vhlevljenin
        self.pripust2n = int(round(self.pripust1n * (self.pripustUp - 1)))

        self.kraveUp = int(self.kraveUpE.text())
        self.bmOdbira = int(self.bmOdbiraE.text())
        self.bmUp = int(self.bmUpE.text())
        self.cak = int(self.cakE.text())
        self.pripustUp = int(self.pripustUpE.text())
        self.pbUp = int(self.pbUpE.text())

        self.pripustDoz = int(self.pripustDozE.text())
        self.mladiDoz = int(self.mladiDozE.text())
        self.pozitivnoTestDoz = int(self.pozitivnoTestDozE.text())
        self.BurnInYN = self.mergeAsk.isChecked()
        self.StBurnInGen = int(self.StBurnInGenE.text())
        self.StSelGen = int(self.SelToGen.text())
        self.NumberOfSires = int(self.NoSires.text())
        self.NumberOfDams = int(self.NoDams.tex())
        self.AlphaSimDir = self.choose_dir()
        self.AlphaSimPed = str(self.AlphaSimPed) + '/SimulatedData/PedigreeAndGeneticValues.txt'
        SpecFile = AlphaSimSpec()

        #ce se izvede tudi burn in
        if self.BurnInYN:
            for roundNo in range(self.StSelGen + 1):
                if roundNo == 0:  # do burn in
                    # prestavi se v AlphaSim Dir
                    os.chdir(self.AlphaSimDir)
                    shutil.copy(SpecFile, self.AlphaSimDir)
                    SpecFile.setPedType("Internal")
                    SpecFile.setBurnInGen(self.StBurnInGen)
                    SpecFile.setSelGen(self.StSelGen)
                    SpecFile.setNoSires(self.NumberOfSires)
                    SpecFile.setNoDams(self.NumberOfDams)
                    SpecFile.turnOnGenFlex()
                    SpecFile.setFlexGenToFrom(1, (self.StBurnInGen + 1))
                    SpecFile.turnOnSelFlex()
                    SpecFile.setExtPedForGen(self.StBurnInGen + 1)
                    SpecFile.setTBVComp(1)
                    SpecFile.setNB(self.stNBn)
                    # pozenes ALPHASIM
                    os.system('./AlphaSim1.05')

        for roundNo in range(1, (self.StSelGen + 1)):
            # precacunas EBV v Ru in zapises PEDIGRE
            ped.computeEBV(0.8)

            if roundNo > 0:
                global ped, categories, sex, active
                ped, categories, sex, active = nastavi_cat('GenPed_EBV.txt', stNB=self.stNBn, \
                                                             nrFn=self.nrFn, nrMn=self.nrMn, \
                                                             telFn=self.telFn, telMn=self.telMn,\
                                                             potomciNPn=self.potomciNPn, \
                                                             vhlevljenin=self.vhlevljenin, ptn=self.ptn, \
                                                             bmn=self.bmn, mladin=self.mladin, \
                                                             bik12n=self.bik12n, pripust1n=self.pripust1n, \
                                                             pripust2n=self.pripust2n, \
                                                             cak=self.cak, kraveUp=self.kraveUp, \
                                                             bmOdbira=self.bmOdbira, bmUp=self.bmUp, \
                                                             pripustDoz=self.pripustDoz, mladiDoz=self.mladiDoz,\
                                                             pozitivnoTestDoz=self.pozitivnoTestDoz, \
                                                             pbUp=self.pbUp)
            elif roundNo > 1:
                global ped, categories, sex, active
                ped, categories, sex, active = selekcija_ena_gen('GenPed_EBV.txt', \
                                                                 categories=categories, sex=sex, \
                                                                 active=active, stNB=self.stNBn, \
                                                                 nrFn=self.nrFn, nrMn=self.nrMn, \
                                                                 telFn=self.telFn, telMn=self.telMn, \
                                                                 potomciNPn=self.potomciNPn, \
                                                                 vhlevljenin=self.vhlevljenin, ptn=self.ptn, \
                                                                 bmn=self.bmn, mladin=self.mladin, \
                                                                 bik12n=self.bik12n, pripust1n=self.pripust1n, \
                                                                 pripust2n=self.pripust2n, \
                                                                 cak=self.cak, kraveUp=self.kraveUp, \
                                                                 bmOdbira=self.bmOdbira, bmUp=self.bmUp, \
                                                                 pripustDoz=self.pripustDoz,
                                                                 mladiDoz=self.mladiDoz, \
                                                                 pozitivnoTestDoz=self.pozitivnoTestDoz, \
                                                                 pbUp=self.pbUp)

                # prestavi se v AlphaSim Dir
                os.chdir(self.AlphaSimDir)
                # kopiraj pedigre v selection folder
                shutil.copy('ExternalPedigree.txt',
                            self.AlphaSimDir + '/Selection/SelectionFolder' + str(roundNo) + '/')
                # TUKAJ POTEM popravis AlphaSimSpec
                # PRVIc PO BURN IN-U
                shutil.copy('/home/jana/Genotipi/Genotipi_CODES/AlphaSimSpec.txt', self.AlphaSimDir)
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



if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = SelParam()
    window.show()
    sys.exit(app.exec_())


