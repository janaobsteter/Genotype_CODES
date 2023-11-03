# -*- coding: utf-8 -*-
from __future__ import division
import os
from selection10 import *
import resource
import pandas as pd


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
                                    self.AlphaSimDir + '/renumf90_noPE.par')  # set levels of random aniaml effect, add var and res var
        # the paramfile is now set
        if self.sel == 'gen':
            GenFiles = snpFiles(self.AlphaSimDir)
            GenFiles.createBlupf90SNPFile()

        os.system("./renumf90 < renumParam_noPE")  # run renumf90

	resource.setrlimit(resource.RLIMIT_STACK, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
        os.system('./blupf90 renf90.par')



homeDir = os.getcwd()
WorkingDir = "/home/v1jobste/JanaO/"

for reference in ["1K", "2K", "3K", "4K", "5K", "6K", "7K", "8K", "9K"]:
    for rep in [0, 1, 2, 3, 4]:
        AlphaSimDir = homeDir + "/FillInBurnIn" + str(rep) + "_permEnv/"
#        AlphaSimDir = homeDir + "/FillInBurnIn" + str(rep) + "/"

        os.chdir(AlphaSimDir)
        print str(os.getcwd())

        #move correspodning IndForGeno_xK.txt to IndForGeno.txt
	os.system("rm IndForGeno.txt")
        os.system("rm GenoFile.txt")
        os.system("cp IndForGeno_" + reference + ".txt  IndForGeno.txt")

	#use module to estimate BV
        blupNextGen = estimateBV(AlphaSimDir, WorkingDir + "/CodeDir", way='milk', sel='gen')
	blupNextGen.computeEBV()
#        blupNextGen.computeEBV_permEnv_herd(setVar=True, varPE=varPE, varE=varEest, varH=varHY,
#                                            repeats=repeats)


        #run MatchAfterRenum --> label as renumbered_SOlutions_xK
        os.system("bash Match_AFTERRenum.sh")
        os.system("mv renumbered_Solutions renumbered_Solutions_" + reference + "_noPE")




