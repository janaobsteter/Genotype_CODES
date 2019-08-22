# -*- coding: utf-8 -*-
from __future__ import division
import os
from selection10 import *
import resource
import pandas as pd

variances = sys.argv[1].split(",")
varPE = float(variances[0])
varH = float(variances[1])
varHY = float(variances[2])
varHTD = float(variances[3])
varE = float(variances[4])


class estimateBV:
    def __init__(self, AlphaSimDir, codeDir, way, sel):
        self.AlphaSimDir = AlphaSimDir
        self.way = way
        self.sel = sel
        self.codeDir = codeDir


    def computeEBV_permEnv_herd(self, setVar=False, varPE=0.0, varE=0.0, herd=True, varH = 0.0, repeats=1):
        """
        This function prepares blupf90 phenotypic .dat and pedigree .ped file according to the specification for a model with permanent Environemtn
        It prepared different files whether we are in fill in or whether in selection gen
        It also prepares different blupf90 spec file whether we are running conventional or genomic prediction
        :param setVar: are you setting the variances manually (permanent env and residual)
        :param varE: what is the proportion of residual variance towards additive genetic variance (Ve / Va)
        :param varPE:what is the proportion of variance for permanent environment towards additive genetic variance (Vpe / Va)
        :param repeats: how many times is the phenotypes measures on a animal in one selection cycle
        :return: prepares .dat, .ped and parameter file and runs blupf90
        """
        # uredi blupparam file
        # get variance components from AlphaSim Output Files
        OutputFiles = AlphaSim_OutputFile(self.AlphaSimDir)
        genvar = OutputFiles.getAddVar()  # dobi additivno varianco
        resvar = OutputFiles.getResVar()  # dobi varianco za ostanek
        permEvar = 0
        herdvar = 0
        if setVar:
            resvar = genvar * varE
            permEvar = genvar * varPE
            herdvar = genvar * varH

        # pripravi fajle za blupf90
        blupFiles = blupf90(self.AlphaSimDir, self.codeDir, way=self.way, permEnv=True, varPE=permEvar, herd=True, varH=herdvar)

        # skopiraj paramFile za renumf90
        if self.sel == 'gen':
            shutil.copy(blupFiles.blupgenParamFile_permEnv_herd,
                        blupFiles.AlphaSimDir + 'renumf90.par')  # skopiraj template blupparam file
        elif self.sel == 'class':
            shutil.copy(blupFiles.blupgenParamFile_Clas_permEnv_herd,
                        blupFiles.AlphaSimDir + 'renumf90.par')  # skopiraj template blupparam file

        blupFiles.prepareParamFiles_permEnv_herd(genvar, permEvar, resvar, herdvar,
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



homeDir = os.getcwd()


for reference in ["1K", "2K", "3K", "4K", "5K"]:
    for rep in [0, 1, 2]:
        AlphaSimDir = homeDir + "/FillInBurnIn" + str(rep) + "_permEnv/"
        #move correspodning IndForGeno_xK.txt to IndForGeno.txt

        GenFiles = snpFiles(AlphaSimDir)
        GenFiles.createBlupf90SNPFile()

        blupNextGen = estimateBV(AlphaSimDir, WorkingDir + "/CodeDir", way='milk', sel=seltype)
        varEest = varE + varH + varHTD
        blupNextGen.computeEBV_permEnv_herd(setVar=True, varPE=varPE, varE=varEest, varH=varHY,
                                            repeats=repeats)

        #run MatchAfterRenum --> label as renumbered_SOlutions_xK
        os.system("bash Match_AFTERRenum.sh")
        os.system("mv renumbered_Solutions renumbered_Solutions_" + reference)




