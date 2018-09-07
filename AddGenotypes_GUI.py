# -*- coding: utf-8 -*-
from __future__ import division
import sys
import os
from PyQt4 import QtGui, QtCore, uic
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from selection import *
from collections import defaultdict
import shutil
import pandas as pd
import numpy as np
import os
import zipfile
import csv
import GenFiles
import commands
import tempfile


#nalozi GUI za selekcijo
qtCreatorFile = '/home/jana/Genotipi/Genotipi_CODES/GenoMergeBox.ui' # Enter file here.
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

class AddGenotypes(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)
        self.date = str(self.Date.text())
        self.Breed.addItems(["Rjava", "Holstein", "Lisasta"])
        self.breed = str(self.Breed.currentText())
        self.AlleleFormat.addItems(["top", "ab", "forward"])
        self.alleleFormat = str(self.AlleleFormat.currentText())
        self.ParentalVerificationYN.stateChanged.connect(self.enableNoSNPs)
        self.merge_ask = self.mergeAsk.isChecked()
        self.rmOriginalZip = self.removeZipAsk.isChecked()
        self.DoEverything.clicked.connect(self.performAddition)

        self.ZipFile.clicked.connect(self.choose_zipFile)
        self.tempDir.clicked.connect(self.choose_tempDir)
        self.ZipDir.clicked.connect(self.choose_zipDir)
        self.MergedDir.clicked.connect(self.choose_mergeDir)
        self.PeddarowDir.clicked.connect(self.choose_peddarowDir)
        self.ZanardiDir.clicked.connect(self.choose_zanardiDir)
        self.CodeDir.clicked.connect(self.choose_codeDir)



    def enableNoSNPs(self):
        if self.ParentalVerificationYN.isChecked():
            self.NoSNPs.setEnabled(True)
            self.NoSNPs.addItems(["100", "200", "800"])
            self.PVSNPs = int(self.NoSNPs.currentText())
        if not self.ParentalVerificationYN.isChecked():
            self.NoSNPs.setEnabled(False)

    def choose_zipFile(self):
        zipfile = QtGui.QFileDialog.getOpenFileName(self, 'Choose Zip File')
        if zipfile:
            self.ZipFileShow.setText(zipfile)

    def choose_tempDir(self):
        tempdir = QtGui.QFileDialog.getExistingDirectory(self, 'Choose Directory')
        if tempdir:
            self.tempDirShow.setText(tempdir)

    def choose_zipDir(self):
        zipdir = QtGui.QFileDialog.getExistingDirectory(self, 'Choose Directory')
        if zipdir:
            self.ZipDirShow.setText(zipdir)

    def choose_mergeDir(self):
        mergedir = QtGui.QFileDialog.getExistingDirectory(self, 'Choose Directory')
        if mergedir:
            self.MergeDirShow.setText(mergedir)

    def choose_peddarowDir(self):
        pdir = QtGui.QFileDialog.getExistingDirectory(self, 'Choose Directory')
        if pdir:
            self.PeddarowDirShow.setText(pdir)

    def choose_zanardiDir(self):
        zdir = QtGui.QFileDialog.getExistingDirectory(self, 'Choose Directory')
        if zdir:
            self.ZanardiDirShow.setText(zdir)

    def choose_codeDir(self):
        cdir = QtGui.QFileDialog.getExistingDirectory(self, 'Choose Directory')
        if cdir:
            self.CodeDirShow.setText(cdir)


    def remove_from_zip(zipfname, *filenames):
        tempdir = tempfile.mkdtemp()
        try:
            tempname = os.path.join(tempdir, 'new.zip')
            with zipfile.ZipFile(zipfname, 'r') as zipread:
                with zipfile.ZipFile(tempname, 'w') as zipwrite:
                    for item in zipread.infolist():
                        if item.filename not in filenames:
                            data = zipread.read(item.filename)
                            zipwrite.writestr(item, data)
            shutil.move(tempname, zipfname)
        finally:
            shutil.rmtree(tempdir)


    def performAddition(self):
        ########################################################
        # set directories and file names
        ########################################################
        self.zip = str(self.ZipFileShow.text())
        self.tempdir = str(self.tempDirShow.text())
        self.ZipDir = str(self.ZipDirShow.text())
        self.MergeDir = str(self.MergeDirShow.text())
        self.peddarow = str(self.PeddarowDirShow.text())
        self.ZanDir = str(self.ZanardiDirShow.text())
        self.codeDir = str(self.CodeDirShow.text())
        print(self.codeDir)
        print(self.ZanDir)
        print(self.peddarow)
        print(self.breed)
        print(self.alleleFormat)
        self.tempDir = self.tempDir + "/Genotipi_" + str(self.date) + "/"
        print(self.tempDir)

        # ask what action does the user want to perform
        # action = raw_input("Do you want to extract SNPs for parental verification  [Y/N] ")


        # File with a list of 800 SNPs for parentage verification
        SNP800 = self.codeDir + '/Names_800SNPs.txt'
        # file with IDs and seq for the animals
        RJ_IDSeq = self.codeDir + "/Rjave_seq_ID.csv"
        # SNP coding
        SNPSifrant = self.codeDir + "/Sifrant_SNP.csv"

        # name of the file
        zipPackage = self.zip
        #########################################################################################################
        ##########################################################################################################
        ##########################################################################################################
        # create dictionaries
        ##########################################################################################################
        ##########################################################################################################

        # create a dictionary of the number of SNPs and corresponding chip names
        chips = GenFiles.chips
        SNP800Sifrant_Dict = GenFiles.SNP800Sifrant_Dict

        GenoFile = defaultdict(set)
        SampleIDs = defaultdict(list)
        PedFiles = defaultdict(list)
        MapFiles = defaultdict(list)
        PedFilesQC = defaultdict(list)
        MapFilesQC = defaultdict(list)
        AllInfo = []
        # dictionary to hold downlaod date of the genotype package
        DateDownloaded = defaultdict(list)
        DateGenotyped = defaultdict(list)

        # list to hold the SNP800 files produced in the run
        SNP800_Peds = []
        SNP800_Maps = []

        # read in animal ID / Seq / DateOfBirth / SexCode table
        # create a dictionary
        Rj_IDSeq_Dict = defaultdict()
        with open(RJ_IDSeq, 'rb') as IDSeq:
            reader = csv.reader(IDSeq, delimiter=',')
            for line in reader:
                Rj_IDSeq_Dict[line[0]] = line[1:]

        ############################################################################################################
        #############################################################################################################
        # create a directory with the current date for temp genotype manipulation
        if not os.path.exists(self.tempDir):
            os.makedirs(self.tempDir)

        # change current working directory to the created directory
        os.chdir(self.tempDir)

        shutil.copy(self.zip, self.tempDir)

        # zipPackages = (filter(lambda x: x.endswith('.zip'), os.listdir(tempDir)))

        # now you have all the files in the same format (the six files zipped)
        # extract the FinalReport from each of them and SNPMap

        # try:
        #    zipPackages.remove(zip_file) #remove original zipfile (zipfolder) from the list
        # except:
        #    pass


        onePackage = GenFiles.genZipPackage(self.zip)
        onePackage.extractFinalReport()
        onePackage.extractSNPMap()
        onePackage.extractSampleMap()

        # check for error IDs and replace the prior identified errouneous IDs
        replaceIDs = [('SI 549458926', 'SI54945829'), ('SI549458926', 'SI54945829'), ('SI4574059', 'SI04574059'),
                      ('SI84048801', 'SI84048802'), ('SI4384195', 'SI04384195'), ('Si24289407', 'SI24289407'),
                      ('SI53595706_201851770050_R08C02', 'SI53595706'),
                      ('SI53595706_201851770081_R03C02', 'SI53595706'), ('SI15036148 (COF)', 'SI15036148'),
                      ('SI85036127 (ASUL)', 'SI85036127'), ('SI55035882 (HRABRI)', 'SI55035882'),
                      ('SI95095002 (KINGSTON', 'SI95095002'),
                      ('SI45094707 (VALDEN)', 'SI45094707'), ('SI34951462 (CAFIERO)', 'SI34951462'),
                      ('SI85026654 (VAUDEK)', 'SI85026654'), ('SI74941696 VASK', 'SI74941696')]
        errorIDs = onePackage.extractErrorNames()  # extract Sample Names if they exist - they shouldnt be in the file
        # to samo, če ti samo prav popravi!!!!!!!!!!!!!
        if errorIDs:
            print (onePackage.name, errorIDs)
            for i in errorIDs:
                os.system('sed -i  "s|' + str(i[0]) + '|' + i[
                    1] + '|g" ' + onePackage.name + "_FinalReport.txt")  # errorIDs are tuples, replace first element witht the second
                os.system('sed -i  "s|' + str(i[0]) + '|' + i[1] + '|g" ' + onePackage.name + '_Sample_Map.txt')
                ###############
        for i in replaceIDs:
            os.system('sed -i  "s|' + i[0] + '|' + i[
                1] + '|g" ' + onePackage.name + "_FinalReport.txt")  # errorIDs are tuples, replace first element witht the second
            os.system('sed -i  "s|' + i[0] + '|' + i[1] + '|g" ' + onePackage.name + '_Sample_Map.txt')

        # copy pedda.param and python script to the current directory
        shutil.copy((self.peddarow + "/peddar.param"), "peddar.param")
        shutil.copy((self.peddarow + "/pedda_row.py"), "pedda_row.py")
        # replace strings with shell command
        os.system(
            'sed -i "s|test_FinalReport.txt|' + onePackage.name + "_FinalReport.txt" + '|g" peddar.param')  # insert FinalReport name into peddar.param
        os.system(
            'sed -i "s|Dominant |Dominant_|g" ' + onePackage.name + "_FinalReport.txt")  # problem Dominant Red with a space
        os.system('sed -i "s|Dominant |Dominant_|g" ' + onePackage.name + '_SNP_Map.txt')  ##problem Dominant Red with a space
        os.system('sed -i "s/test_outputfile/"' + onePackage.name + '"/g" peddar.param')  # insert OutPut name into peddar.param
        os.system(
            'sed -i "s/test_SNPMap.txt/"' + onePackage.name + '_SNP_Map.txt' + '"/g" peddar.param')  # insert SNPMap name into peddar.param
        os.system(
            'sed -i "s/AlleleFormat/"' + self.alleleFormat + '"/g" peddar.param')  # insert desired AlleleFormat name into peddar.param
        os.system('sed -i "s/TEST/"' + self.breed + '"/g" peddar.param')
        os.system("python pedda_row.py")  # transform into ped and map file

        # ABFORMAT
        # shutil.copy((self.peddarow + "/peddar.param"), "peddar.param")
        # shutil.copy((self.peddarow + "/pedda_row.py"), "pedda_row.py")
        # # replace strings with shell command
        # os.system(
        #     'sed -i "s|test_FinalReport.txt|' + onePackage.name + "_FinalReport.txt" + '|g" peddar.param')  # insert FinalReport name into peddar.param
        # os.system(
        #     'sed -i "s|Dominant |Dominant_|g" ' + onePackage.name + "_FinalReport.txt")  # problem Dominant Red with a space
        # os.system('sed -i "s|Dominant |Dominant_|g" ' + onePackage.name + '_SNP_Map.txt')  ##problem Dominant Red with a space
        # os.system(
        #     'sed -i "s/test_outputfile/"' + onePackage.name + "_AB" + '"/g" peddar.param')  # insert OutPut name into peddar.param
        # os.system(
        #     'sed -i "s/test_SNPMap.txt/"' + onePackage.name + '_SNP_Map.txt' + '"/g" peddar.param')  # insert SNPMap name into peddar.param
        # os.system('sed -i "s/AlleleFormat/"' + "ab" + '"/g" peddar.param')  # insert desired AlleleFormat name into peddar.param
        # os.system('sed -i "s/TEST/"' + pasma + '"/g" peddar.param')
        # os.system("python pedda_row.py")  # transform into ped and map file

        # create a new zip file with corrected error names
        # shutil.move(onePackage.name+'_Sample_Map.txt', 'Sample_Map.txt') #rename extracted SampleMap
        with zipfile.ZipFile(onePackage.name + '_FinalReport.zip', 'w', zipfile.ZIP_DEFLATED) as myzip:
            myzip.write(onePackage.name + '_FinalReport.txt')  # create new FinalReport zip
        with zipfile.ZipFile('Sample_Map.zip', 'w', zipfile.ZIP_DEFLATED) as myzip:
            myzip.write(onePackage.name + '_Sample_Map.txt')  # create new Sample_Map.zip

        self.remove_from_zip(onePackage.zipname, onePackage.finalreportname)
        self.remove_from_zip(onePackage.zipname, onePackage.samplemapname)

        shutil.move(onePackage.name + '_FinalReport.zip', onePackage.finalreportname)
        shutil.move('Sample_Map.zip', onePackage.samplemapname)

        with zipfile.ZipFile(onePackage.zipname, 'a', zipfile.ZIP_DEFLATED) as z:
            z.write(onePackage.finalreportname)
        with zipfile.ZipFile(onePackage.zipname, 'a', zipfile.ZIP_DEFLATED) as z:
            z.write(onePackage.samplemapname)

        # make pedfile a GenFiles pedFile object
        pedfile = GenFiles.pedFile(onePackage.name + '.ped')
        mapfile = GenFiles.mapFile(onePackage.name + '.map')

        # Perform QC!
        os.system("bash " + self.codeDir + "/1_QC_FileArgs.sh " + pedfile.name + " " + pedfile.chip)
        PedFilesQC[pedfile.chip].append(self.tempDir + pedfile.name + "_" + pedfile.chip + "_CleanIndsMarkers.ped")
        MapFilesQC[pedfile.chip].append(self.tempDir + pedfile.name + "_" + pedfile.chip + "_CleanIndsMarkers.map")

        # add file to the dictionary of chip files
        PedFiles[pedfile.chip].append(self.tempDir + pedfile.pedname)
        MapFiles[pedfile.chip].append(self.tempDir + mapfile.mapname)
        GenoFile[pedfile.chip].add(pedfile.name)
        DateDownloaded[date] += (pedfile.name)
        DateGenotyped[onePackage.genodate] += [(x, pedfile.chip) for x in (pedfile.samples)]
        AllInfo += [(x, pedfile.chip, pedfile.name, onePackage.genodate) for x in (pedfile.samples)]
        for i in pedfile.samples:
            if i in Rj_IDSeq_Dict:
                SampleIDs[i] = [i, Rj_IDSeq_Dict.get(i)[0], onePackage.genodate, pedfile.chip, self.date]
            else:
                print "Sample ID " + i + " in " + pedfile.name + " not found!!!"

        ################################################################################################
        ###############################################################################################
        # END OF THE LOOP
        # merge produced SNP800 files
        # merge ped files if merge_ask = Y
        # create table for govedo
        #############################################################################################
        ###############################################################################################
        print "The number of genotyped animals is {}.".format(len(SampleIDs))
        print "The number of genotype packages (different date of genotyping) is {}.".format(len(DateGenotyped))
        print "The number of different genotyping chips is {0}: {1}.".format(len(PedFiles), PedFiles.keys())

        # Perform QC!!!

        # #create a table of individuals for govedo
        # #columns are seq, chip, date genotyped
        GenotypedInd = pd.DataFrame.from_dict(SampleIDs, orient='index', dtype=None)
        GenotypedInd.columns = ['ID', 'ZIV_ID_SEQ', 'GenoDate', 'Chip', 'DownloadDate']
        imiss = pd.read_table(self.tempDir + pedfile.name + "_" + pedfile.chip + ".imiss", sep="\s+")[["IID", "F_MISS"]]
        imiss.columns = ['ID', "F_MISS"]
        Tabela = pd.merge(GenotypedInd, imiss, on="ID")
        Tabela.to_csv(path_or_buf=self.tempDir + str(onePackage.genodate) + 'GovedoInd.csv', sep=",", index=False)
        print("Created table for Govedo.")
        #



        # merge is outside the loop
        # merge all the chips needed updating

        for i in PedFiles:
            if not os.path.exists(self.MergeDir + str(i)):
                os.makedirs(self.MergeDir + str(i))
            for pedfile, mapfile in zip(PedFiles[i], MapFiles[i]):
                shutil.copy(pedfile, self.MergeDir + str(i))
                shutil.copy(mapfile, self.MergeDir + str(i))
            os.chdir(self.MergeDir + str(i))
            shutil.copy(self.codeDir + "/PARAMFILE.txt", self.MergeDir + i)
            pedToMerge = ",".join(PedFiles[i]).strip("'")
            mapToMerge = ",".join(MapFiles[i]).strip("'")
            if not os.path.isfile(self.MergeDir + i + '/PLINK_MERGED.ped'):
                mergeChipCommand = "plink --file {0} --cow --merge-list {1} --recode --out PLINK_MERGED".format(
                    (PedFiles[i][0].strip(".ped")), 'MergeChip.txt')
                with open('MergeChip.txt', 'w') as csvfile:
                    writer = csv.writer(csvfile, delimiter=" ")
                    [writer.writerow(r) for r in
                     zip(PedFiles[i][1:], MapFiles[i][1:])]  # leave the first one out - that goes in the plink command line
            if os.path.isfile(self.MergeDir + i + '/PLINK_MERGED.ped'):
                mergeChipCommand = "plink --file PLINK_MERGED --cow --merge-list {0} --recode --out PLINK_MERGED".format(
                    'MergeChip.txt')
                with open('MergeChip.txt', 'w') as csvfile:
                    writer = csv.writer(csvfile, delimiter=" ")
                    [writer.writerow(r) for r in zip(PedFiles[i], MapFiles[i])]

            status, output = commands.getstatusoutput(mergeChipCommand)  # merge with plink

            if status == 0:
                print "Successfully merged " + str(i) + " " + self.MergeDir + " " + i
            else:
                print "Merging went wrong, error: " + str(status)

        for chip in PedFiles:
            PedFiles[chip] = [i.replace("ZipGenoFiles", "ZipGenoFiles/") for i in PedFiles[chip]]

        for chip in MapFiles:
            MapFiles[chip] = [i.replace("ZipGenoFiles", "ZipGenoFiles/") for i in MapFiles[chip]]




        # MERGE FOR QC-ed data!!!!

        # for i in PedFiles:
        #    if not os.path.exists(PLINKDIR+str(i)):
        #        os.makedirs(PLINKDIR+str(i))
        #    for pedfile, mapfile in zip (PedFilesQC[i], MapFilesQC[i]):
        #        shutil.copy(pedfile, PLINKDIR+str(i))
        #        shutil.copy(mapfile, PLINKDIR+str(i))
        #    os.chdir(PLINKDIR+str(i))
        #    shutil.copy("/home/jana/Genotipi/Genotipi_CODES/PARAMFILE.txt", PLINKDIR+i)
        #    pedToMerge = ",".join(PedFilesQC[i]).strip("'")
        #    mapToMerge = ",".join(MapFilesQC[i]).strip("'")
        #    if not os.path.isfile(PLINKDIR+i+'/PLINK_MERGED_' + i + '_CleanIndsMarkers.ped'):
        #        mergeChipCommand = "plink --file {0} --cow --merge-list {1} --recode --out {2}".format((PedFilesQC[i][0].strip(".ped")), 'MergeChip.txt', "PLINK_MERGED_" + i + "_CleanIndsMarkers")
        #        with open('MergeChip.txt', 'w') as csvfile:
        #            writer = csv.writer(csvfile, delimiter=" ")
        #            [writer.writerow(r) for r in zip(PedFilesQC[i][1:], MapFilesQC[i][1:])] #leave the first one out - that goes in the plink command line
        #    if os.path.isfile(PLINKDIR+i+'/PLINK_MERGED_' + i + '_CleanIndsMarkers.ped'):
        #        mergeChipCommand = 'plink --file PLINK_MERGED_{0}_CleanIndsMarkers --cow --merge-list {1} --recode --out PLINK_MERGED_{0}_CleanIndsMarkers'.format(i, 'MergeChip.txt')
        #        with open('MergeChip.txt', 'w') as csvfile:
        #            writer = csv.writer(csvfile, delimiter=" ")
        #            [writer.writerow(r) for r in zip(PedFilesQC[i], MapFilesQC[i])]
        #
        #    status, output = commands.getstatusoutput(mergeChipCommand) #merge with plink
        #
        #    if status == 0:
        #        print "Successfully merged " + str(i) + " " + PLINKDIR + " " + i + "_CleanIndsMarkers"
        #    else:
        #        print "Merging went wrong, error: " + str(status)



if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = AddGenotypes()
    window.show()
    sys.exit(app.exec_())