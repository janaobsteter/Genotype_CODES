#This is a script to add newly genotyped individuals and downloaded GeneSeek zip file (Final Reports)
#to the existing database of the latest genotypes

#Pipeline:
#1)define chip, dictionary, dictionary to hold chip: genotype package: animal ids, dictionary to hold genotype package : download date
#2)create temp directory within breed_TEMP/DownloadDate
#3)create directory if not existin, unzip file
#4)for each genotype package: unzip Final report and SNP_Map, change spurious strings within the files and adjust peddar.param file
#run peddar_row to transform FinalReports to PLINK and MAP formats
#write names to dictionaries


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
########################################################
#set directories and file names
########################################################
date=23102017
pasma="Rjava"
AlleleFormat="top"
zip_file=""
merge_ask='Y'
#Ask the user for the current date (date of download) and breed
#date = raw_input("Enter the date (today): ")
#pasma = raw_input("Enter the breed [Rjava/Holstein/Lisasta]: ") 
#AlleleFormat=raw_input("Enter the desired allele coding [top / forward / ab]: ")
#zip_file = raw_input("Enter the name of the downloaded zip file: ")
#merge_ask=raw_input("Do you want to merge newly downloaded genotypes to the Latest Genotypes files (by chip)? [Y/N] ")

#ask what action does the user want to perform
#action = raw_input("Do you want to extract SNPs for parental verification  [Y/N] ")
action = 'N'
if action == 'Y':
    PVSNPs = input("How many SNPs would you like to use for parental verification? ")
#ask whether you want to remove original zip
rmOriginalZip=raw_input('Remove original zip? [Y/N] ')
#create directory path to hold current temp genotype files within Genotipi_DATA and breed directory
tempDir = "/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/Top/ZipGenoFiles/"
#PEDDAROW directory
peddarow="/home/jana/Genotipi/Genotipi_CODES/SNPchimpRepo/source_codes/PEDDA_ROW/"
#Zip latest
#Zip_lat="/run/user/1000/gvfs/smb-share:server=kis-h2.si,share=kisdfs/ZIV/vol1/ZIV/VSI/jana/"+pasma+"/ZipGenoFiles/"
Zip_lat="/home/jana/Genotipi/TEST/"
PLINKDIR = '/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/Top/Top_PLINK/'
#Genotipi_latest directory
Gen_lat = "/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/"+pasma+"/"
#path to Zanardi
ZanDir="/home/jana/Genotipi/Genotipi_CODES/Zanardi/"

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

GenoFile = defaultdict(set)
SampleIDs = defaultdict(list)
PedFiles = defaultdict(list)
MapFiles = defaultdict(list)
AllInfo = []
#dictionary to hold downlaod date of the genotype package   
DateDownloaded = defaultdict(list)
DateGenotyped = defaultdict(list)

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
        

############################################################################################################
#############################################################################################################
#create a directory with the current date for temp genotype manipulation
if not os.path.exists(tempDir):
    os.makedirs(tempDir)


#change current working directory to the created directory
os.chdir(tempDir)



    
    
zipPackages = (filter(lambda x: x.endswith('.zip'), os.listdir(tempDir)))

#now you have all the files in the same format (the six files zipped)
#extract the FinalReport from each of them and SNPMap
  
try:
    zipPackages.remove(zip_file) #remove original zipfile (zipfolder) from the list
except:
    pass

for zipPackage in zipPackages:
    onePackage=GenFiles.genZipPackage(zipPackage)
    onePackage.extractFinalReport()
    onePackage.extractSNPMap()
    onePackage.extractSampleMap()
    
    #check for error IDs and replace the prior identified errouneous IDs
    replaceIDs = [('SI4574059','SI04574059'),('SI84048801','SI84048802'),('SI4384195','SI04384195'),('Si24289407','SI24289407')]
    errorIDs = onePackage.extractErrorNames() #extract Sample Names if they exist - they shouldnt be in the file
    if errorIDs:
        print (onePackage.name, errorIDs)
        for i in errorIDs:
            os.system('sed -i  "s|' +str(i[0])+ '|' + i[1] + '|g" ' + onePackage.name+"_FinalReport.txt") #errorIDs are tuples, replace first element witht the second
            os.system('sed -i  "s|' +str(i[0])+ '|' + i[1] + '|g" '+onePackage.name+'_Sample_Map.txt') 
    for i in replaceIDs:
        os.system('sed -i  "s|' +i[0]+ '|' + i[1] + '|g" ' + onePackage.finalreportname.strip('zip') + 'txt') #errorIDs are tuples, replace first element witht the second
        os.system('sed -i  "s|' +i[0]+ '|' + i[1] + '|g" '+onePackage.name + '_Sample_Map.txt')
                 
                    
    #copy pedda.param and python script to the current directory
    shutil.copy((peddarow+"/peddar.param"), "peddar.param")
    shutil.copy((peddarow+"/pedda_row.py"), "pedda_row.py")
    #replace strings with shell command
    os.system('sed -i "s/test_FinalReport.txt/'+ onePackage.finalreportname.strip('zip') + 'txt/g" peddar.param') #insert FinalReport name into peddar.param
    os.system('sed -i "s|Dominant |Dominant_|g" ' + onePackage.finalreportname.strip('zip') + 'txt') #problem Dominant Red with a space
    os.system('sed -i "s|Dominant |Dominant_|g" ' + onePackage.name+'_SNP_Map.txt') ##problem Dominant Red with a space
    os.system('sed -i "s/test_outputfile/"'+onePackage.name+'"/g" peddar.param') #insert OutPut name into peddar.param
    os.system('sed -i "s/test_SNPMap.txt/"'+onePackage.name+'_SNP_Map.txt'+'"/g" peddar.param') #insert SNPMap name into peddar.param
    os.system('sed -i "s/AlleleFormat/"'+AlleleFormat+'"/g" peddar.param') #insert desired AlleleFormat name into peddar.param
    os.system('sed -i "s/TEST/"'+pasma+'"/g" peddar.param')
    os.system("python pedda_row.py") #transform into ped and map file
    
    
    #make pedfile a GenFiles pedFile object
    pedfile=GenFiles.pedFile(onePackage.name +  '.ped')
    mapfile=GenFiles.mapFile(onePackage.name +  '.map')
    
    #add file to the dictionary of chip files
    PedFiles[pedfile.chip].append(tempDir+pedfile.pedname)
    MapFiles[pedfile.chip].append(tempDir+mapfile.mapname)
    GenoFile[pedfile.chip].add(pedfile.name)
    DateDownloaded[date] += (pedfile.name)
    DateGenotyped[onePackage.genodate] += [(x, pedfile.chip) for x in (pedfile.samples)]
    AllInfo += [(x, pedfile.chip, pedfile.name, onePackage.genodate) for x in (pedfile.samples)]
    for i in pedfile.samples:
        if i in Rj_IDSeq_Dict:
            SampleIDs[i] = [Rj_IDSeq_Dict.get(i)[0], onePackage.genodate, pedfile.chip, date]
        else: 
            print "Sample ID " + i + " not found!!!"
   


    

################################################################################################
###############################################################################################
#END OF THE LOOP
#merge produced SNP800 files
#merge ped files if merge_ask = Y
#create table for govedo
#############################################################################################
###############################################################################################
print "The number of genotyped animals is {}.".format(len(SampleIDs))
print "The number of genotype packages (different date of genotyping) is {}.".format(len(DateGenotyped))
print "The number of different genotyping chips is {0}: {1}.".format(len(PedFiles), PedFiles.keys())


#create a table of individuals for govedo   
#columns are seq, chip, date genotyped  
#GenotypedInd = pd.DataFrame.from_dict(SampleIDs, orient='index', dtype=None)
#GenotypedInd.columns = ['ZIV_ID_SEQ','GenoDate','Chip','DownloadDate']
#GenotypedInd.to_csv(path_or_buf = tempDir+str(date)+'GovedoInd.csv', sep=",", index=False  )



#merge is outside the loop
#merge all the chips needed updating

for i in PedFiles:
    if not os.path.exists(PLINKDIR+str(i)):
        os.makedirs(PLINKDIR+str(i))
    for pedfile, mapfile in zip (PedFiles[i], MapFiles[i]):
        shutil.copy(pedfile, PLINKDIR+str(i))
        shutil.copy(mapfile, PLINKDIR+str(i))
    os.chdir(PLINKDIR+str(i))
    shutil.copy("/home/jana/Genotipi/Genotipi_CODES/PARAMFILE.txt", PLINKDIR+i)
    pedToMerge = ",".join(PedFiles[i]).strip("'")
    mapToMerge = ",".join(MapFiles[i]).strip("'")
    if not os.path.isfile(PLINKDIR+i+'/PLINK_MERGED.ped'):
        mergeChipCommand = "plink --file {0} --cow --merge-list {1} --recode --out PLINK_MERGED".format((PedFiles[i][0].strip(".ped")), 'MergeChip.txt')
        with open('MergeChip.txt', 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter=" ")
            [writer.writerow(r) for r in zip(PedFiles[i][1:], MapFiles[i][1:])] #leave the first one out - that goes in the plink command line
    if os.path.isfile(PLINKDIR+i+'/PLINK_MERGED.ped'):
        mergeChipCommand = "plink --file PLINK_MERGED --cow --merge-list {0} --recode --out PLINK_MERGED".format('MergeChip.txt')
        with open('MergeChip.txt', 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter=" ")
            [writer.writerow(r) for r in zip(PedFiles[i], MapFiles[i])] 
        
    status, output = commands.getstatusoutput(mergeChipCommand) #merge with plink
    
    if status == 0:
        print "Successfully merged " + str(i)
    else:
        print "Merging went wrong, error: " + str(status)

for chip in PedFiles:
    PedFiles[chip] = [i.replace("ZipGenoFiles", "ZipGenoFiles/") for i in PedFiles[chip]]
    
    
for chip in MapFiles:
    MapFiles[chip] = [i.replace("ZipGenoFiles", "ZipGenoFiles/") for i in MapFiles[chip]]
    
    