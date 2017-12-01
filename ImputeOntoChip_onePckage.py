# -*- coding: UTF-8 -*-
#this is a script to impute all genotypes onto IDBv03


import os
import sys
import zipfile
import shutil
from collections import defaultdict
import csv
import GenFiles_Gen
from GenFiles_Gen import *
import commands
import tempfile
import pandas as pd
from itertools import chain

class ZanardiParamFile():
    def __init__(self, Path, WorkingDir):
        self.WorkingDir = WorkingDir
        self.paramFile_orig = Path + 'PARAMFILE.txt'  
        shutil.copy(Path +'PARAMFILE.txt', WorkingDir + 'PARAMFILE.txt')
        self.paramFile_working = WorkingDir + 'PARAMFILE.txt'  
        
    def setPedMap(self, name):
        os.system('sed -i "s|PathToPed|' + name + '.ped|g" ' + self.paramFile_working)
        os.system('sed -i "s|PathToMap|' + name + '.map|g" ' + self.paramFile_working)
                        
    def setOutName(self, outName):
        os.system('sed -i "s|OutputName|' + outName + '|g" ' + self.paramFile_working)
    
    def setPedigree(self, pedigree):
        os.system('sed -i "s|PedigreePath|' + pedigree + '|g" ' + self.paramFile_working)

    def setRef(self, refList):
        os.system('''sed -i 's|FImputeOptions|ref="''' + refList + '''";|g' ''' + self.paramFile_working)

reload(sys.modules['GenFiles_Gen'])
#set code directory
CodeDir = '/home/jana/Genotipi/Genotipi_CODES/'
#zanardi direcory
ZanDir = '/home/jana/Genotipi/Genotipi_CODES/Zanardi/'
#set working directory
WorkingDir = '/home/jana/Genotipi/Genotipi_WORK/ImputationOntoIDBv03/'
#set directory with the chip PLINK files
GenDir = '/run/user/1000/gvfs/smb-share:server=kis-h2.si,share=kisdfs/ZIV/vol1/ZIV/VSI/JanaO/Genotipi/Genotipi_latest/Rjava/'
#set reference chip to impute on
RefChip = 'IDBv03'
#this are the reference SNPs - the ones you have to extract from all other chips
RefMap = mapFile(GenDir + RefChip + '/PLINK_MERGED_'+RefChip+'_CleanIndsMarkers.map')
RefPed = pedFile(GenDir + RefChip + '/PLINK_MERGED_'+RefChip+'_CleanIndsMarkers.ped')
RefSNPs = list(RefMap.snps)

#a list of chips to impute
chipImpute = ['GGPv02', 'GGPv03', 'GGPv04', 'HD', 'HDv02', '50Kv01', '50Kv02']

#here extract reference SNPs from each chip PLINK file, write them to chip_RefChip.ped
#also extract all individuals that need to be imputed - ImputeInds
ImputeInds = []
os.chdir(WorkingDir)
for chip in chipImpute:
    chipPed = pedFile(GenDir + chip + '/PLINK_MERGED_'+chip+'_CleanIndsMarkers.ped')
    ImputeInds.append(chipPed.samples) 
    chipPed.extractSNPList_Binary(RefSNPs, chip+'_'+RefChip)

#merge all the extracted chip PLINK files together (all contain only RefChip SNPs)
mergeList(RefMap.name, [chip+'_'+RefChip for chip in chipImpute])

#now make Cluster SNPs files, make cluster number 1 for all of them
for clustN in range(10):
    pd.DataFrame({1: RefSNPs[clustN::10], 2: [1] * len(RefSNPs[clustN::10])}).to_csv(WorkingDir + 'ClusterSNPs_' + str(clustN) + '.txt', header=None, index=None, sep=" ")

#make one Cluster Individuals file with all the indivuduals except for the ones on IDBv03
pd.DataFrame({'FID': ['Rjava'] * len(set(list(chain.from_iterable(ImputeInds)))), 'Inds': list(set(list(chain.from_iterable(ImputeInds)))), 'Cluster': [1] * len(set(list(chain.from_iterable(ImputeInds))))}).to_csv(WorkingDir + 'ImputeInds.txt', header=None, index=None, sep=" ")
#make a file with IDs of the reference chip
pd.DataFrame({'Inds': RefPed.samples}).to_csv(WorkingDir + 'RefInds.txt', header=None, index=None, sep=" ")



#compute allelic concordance!
Concordance = []
#from here if masking and imputation to obtain accuracies of imputation
for masking in range(10):
    os.system('plink --file MERGEDforImputation --cow --chr 1-29 --zero-cluster ' + WorkingDir + 'ClusterSNPs_' + str(masking) + '.txt --within ' + WorkingDir + 'ImputeInds.txt --make-bed --out Masked' + str(masking))
    os.system('plink --bfile  Masked' + str(masking) + ' --cow --recode --out Masked' + str(masking))
    ParamFile = ZanardiParamFile(CodeDir, WorkingDir)
    ParamFile.setPedMap(WorkingDir + 'Masked' + str(masking))
    ParamFile.setOutName('Imputed' + str(masking))
    ParamFile.setPedigree(CodeDir + 'RJPedigre_ZANARDI')
    ParamFile.setRef(WorkingDir + 'RefInds.txt')
    try:
        os.system('python Zanardi.py --fimpute --save')
    except:
        os.system('cut -f1 OUTPUT/Error_SNP_position.txt -d";"| tail -n +2 > tmp && mv tmp ErrorSNPs.txt')
        os.system('plink --file MERGEDforImputation --cow --exclude ErrorSNPs.txt --recode --out MERGEDforImputation')
        os.system('python Zanardi.py --fimpute --save')

    #Now extract only reference individuals from the imputed file
    os.system('plink --file ' + WorkingDir + '/OUTPUT/FIMPUTE_Imputed' + str(
        masking) + ' --cow --keep ImputedIndsConc.txt --recode --out ImputedMasking' + str(masking))
    os.system(
        'plink --file ' + WorkingDir + '/MERGEDforImputation --cow --keep ImputedIndsConc.txt --recode --out OrigMasking' + str(
            masking))

    # Now compute the allelic concordance
    chkConc = AllelicConcordance('OrigMasking' + str(masking) + '.ped', 'ImputedMasking' + str(masking) + '.ped')
    os.system('cut -f1 ClusterSNPs_' + str(masking) + '.txt -d " " > ClusterList' + str(masking) + '.txt')
    chkConc.concSNPs('ClusterList' + str(masking) + '.txt', 'Concordance' + str(masking))
    Concordance.append(chkConc.extractConc('Concordance' + str(masking) + '.txt'))

Concordance.append(mean(Concordance))
pd.DataFrame({"Conc": Concordance}).to_csv(WorkingDir + "FinalAllelicConcordance.txt", index=None)

#get genotypic concordance - bash skript
os.system ("bash ~/Genotipi/Genotipi_CODES/OneMaskingGenCor.sh")
#output files are names GenCorr_samples and GenCorr_snps
os.system("cat GenCorr_samples* > GenCorr_SAMPLES.txt")
print("Genetic Samples Correlation is " + str(mean(pd.read_table(WorkingDir + "GenCorr_SAMPLES.txt", sep=" ", header=None)[1])))
os.system("cat GenCorr_snps* > GenCorr_SNPS.txt")
print("Genetic SNP Correlation is " + str(mean(pd.read_table(WorkingDir + "GenCorr_snps.txt", sep=" ", header=None)[1])))


#potem pa pripravi še blupf90 files
blupFile = blupf90(WorkingDir + 'OUTPUT/FIMPUTE_ImputedFULL.ped')
blupFile.createBlupGeno("BlupIMPUTED")