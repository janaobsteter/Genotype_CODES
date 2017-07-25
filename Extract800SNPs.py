import os
import sys
import zipfile
import shutil
from collections import defaultdict
import csv
import GenFiles
import commands
import tempfile
import pandas as pd
import os
import sys


os.chdir("/home/jana/Mreza/JanaO/Genotipi/Genotipi_latest/Rjava/GGPv03")
PLINK_name = "PLINK_MERGED_GGPv03_CleanIndsMarkers"
#PEDDAROW directory
peddarow="/home/jana/Genotipi/Genotipi_CODES/SNPchimpRepo/source_codes/PEDDA_ROW"
#File with a list of 800 SNPs for parentage verification
SNP800="/home/jana/Genotipi/ParentalVerification_SNPSNP/Names_800SNPs.txt"
#file with IDs and seq for the animals
RJ_IDSeq="/home/jana/Genotipi/Genotipi_CODES/Rjave_seq_ID.csv"
#SNP coding
SNPSifrant="/home/jana/Genotipi/ParentalVerification_SNPSNP/Sifrant_SNP.csv"

#########################################################################################################
##########################################################################################################
##########################################################################################################
#create dictionaries
##########################################################################################################
##########################################################################################################

#create a dictionary of the number of SNPs and corresponding chip names
chips = GenFiles.chips
SNP800Sifrant_Dict = GenFiles.SNP800Sifrant_Dict


#list to hold the SNP800 files produced in the run
SNP800_Peds=[]
SNP800_Maps=[]


#read in animal ID / Seq / DateOfBirth / SexCode table
#create a dictionary
Rj_IDSeq_Dict = defaultdict()
with open(RJ_IDSeq, 'rb') as IDSeq:
    reader = csv.reader(IDSeq, delimiter=',')
    for line in reader:
        Rj_IDSeq_Dict[line[0]] = line[1:]
        





#make pedfile a GenFiles pedFile object
pedfile=GenFiles.pedFile(PLINK_name +  '.ped')
mapfile=GenFiles.mapFile(PLINK_name +  '.map')

SampleIDs = defaultdict(list)

#add file to the dictionary of chip files
for i in pedfile.samples:
    if i in Rj_IDSeq_Dict:
        SampleIDs[i] = [Rj_IDSeq_Dict.get(i)[0]]
    else: 
        print "Sample ID " + i + " not found!!!"

#extract the 800 SNPs
pedfile.extractParentalSNPs(800)

   
#create a table of SNPs for parental verification
parentalTable = pd.read_table('ParentalSNP800_' + str(pedfile.sernum) + '.ped', sep=" ", header=None)
parentalTable = parentalTable.drop(parentalTable.columns[[0,2,3,4,5]], axis=1)

#ind800SNP = pd.DataFrame({'ID' = Rj_IDSeq_Dict['SI04640608']*1600, 
ped800 = GenFiles.pedFile('ParentalSNP800_' + str(pedfile.sernum) + '.ped')
map800 = GenFiles.mapFile('ParentalSNP800_' + str(pedfile.sernum) + '.map')

AlleleCodes = []
for i in map800.snps:
    AlleleCodes.append((i, min(SNP800Sifrant_Dict[i])))
    AlleleCodes.append((i, max(SNP800Sifrant_Dict[i])))
    
    
IndSNPs = []
for sampleID in ped800.samples:
    for (value, (snpname, code)) in zip(ped800.individualSNPsList(sampleID), AlleleCodes):
        try:
            IndSNPs.append((Rj_IDSeq_Dict.get(sampleID)[0], value, code))
        except:
            print sampleID


IndSNPsDF = pd.DataFrame(IndSNPs)
IndSNPsDF.columns =["ZIV_ID_SEQ", "VREDNOST", "SIFRA"]
IndSNPsDF["VREDNOST"].replace(to_replace='0', value='', inplace=True)
IndSNPsDF.to_csv('GovedoParSNPs.csv', sep=",", index=False )





 
#check whether the number of individuals in all genoPackages corresponds to the number of individuals with extracted 800SNPs for parentage testing
if len(SampleIDs) == len(ped800.samples):
    print "The number of individuals corresponds to the number of individuals in SNP800 ped file"
else: 
    print "The number of individuals DOES NOT correspond to the number of individuals in SNP800 ped file, check for ERRORS!"


print "Succesfully created SNPfile for Govedo Database: GovedoParSNPs.csv"





