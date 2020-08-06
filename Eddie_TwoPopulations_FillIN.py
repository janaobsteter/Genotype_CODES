# -*- coding: utf-8 -*-
from __future__ import division
import sys
import os
import math
from selection10 import *
from selection10 import nastavi_cat_TGV, selekcija_total_TGV
from collections import defaultdict
import shutil
import pandas as pd
import numpy as np
import resource
import ast
from random import randint


WorkingDir = "/home/v1jobste/jobsteter/"


######################################################################################
######################################################################################
#scenarios = ['Class_LargePop'] #, 'GenSLO', 'BmGen', 'OtherCowsGen', 'Gen']
REP = sys.argv[1]



for rep in [REP]:
    #####################################################################################################
    #####################################################################################################
    #FILL IN - 2x newborns
    #####################################################################################################
    if not os.path.isdir("FillIn_TwoPop_" + str(rep)):
        os.makedirs("FillIn_TwoPop_" + str(rep))
    os.chdir("FillIn_TwoPop_" + str(rep)) #prestavi se v FillInBurnin za ta replikat
    os.system('cp -r ' + WorkingDir + '/Essentials/* .') # skopiraj vse iz Esentials
    os.system('cp -r ' + WorkingDir + '/CodeDir/* .') # skopiraj vse iz CodeDir
    seed =  randint(-100000000, -1)
    os.system("echo " + str(seed) + " > Seed.txt")

    #first make a FILLIN
    #nastavi AlphaSimSpec
    print(os.getcwd())

    SpecFile = AlphaSimSpec(os.getcwd(), WorkingDir + "/CodeDir", type="Multitrait")  # AlphaSimSpec je class iz selection, ki omogoča nastavljanje parametrov AlphaSimSpec fila
    SpecFile.setPedType("Internal")  # pedigree je za burn in internal
    SpecFile.setNB(17280)  # stevilo novorojenih
    SpecFile.setBurnInGen(20)  # stevilo burnINGen
    SpecFile.setSelGen(40)  # st selection gen
    SpecFile.setNoSires(20)
    SpecFile.setNoDams(8640)
    SpecFile.turnOnGenFlex()
    SpecFile.setFlexGenToFrom(1, 21)  # pozeni od generacije 1 do burnin+1
    SpecFile.turnOnSelFlex()
    SpecFile.setExtPedForGen(21)  # za katero generacijo uvozi external pedigre  - ena po burn in
    SpecFile.setTBVComp(1)  # skomputiraj TBV
    # pozenes ALPHASIM
    os.system('./AlphaSim1.08')






