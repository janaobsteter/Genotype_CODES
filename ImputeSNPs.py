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


#######################################################
#create dictionaries
#######################################################
#create a dictionary of the number of SNPs and corresponding chip names
chips=GenFiles.chips

GenoFile = defaultdict()
PedFiles = defaultdict()
MapFiles = defaultdict()

for each in set(chips.values()):
    GenoFile[each] = {}
    PedFiles[each] = []
    MapFiles[each] = []

#dictionary to hold download date of the genotype package   
DateDownloaded = defaultdict()

#list to hold the SNP800 files produced in the run
SNP800_Peds=[]
SNP800_Maps=[]

########################################################
#set directories and file names
########################################################
date=13012017
pasma="Rjava"
AlleleFormat="ab"
zip_file="Matija_Rigler.zip"
#Ask the user for the current date (date of download) and breed
#date = input("Enter the date (today): ")
#pasma = raw_input("Enter the breed: ") 
#AlleleFormat=raw_input("Enter the desired allele coding [top / forward / ab]: ")
#zip_file = raw_input("Enter the name of the downloaded zip file: ")
merge_ask=raw_input("Do you want to merge newly downloaded genotypes to the Latest Genotypes files (by chip)? [Y/N] ")
#create directory path to hold current temp genotype files within Genotipi_DATA and breed directory
tempDir = "/home/janao/Genotipi/Genotipi_DATA/"+pasma+"_TEMP/Genotipi_"+str(date)
#PEDDAROW directory
peddarow="/home/janao/Genotipi/Genotipi_CODES/SNPchimpRepo/source_codes/PEDDA_ROW"
#Genotipi_latest directory
Gen_lat = "/home/janao/Genotipi/Genotipi_DATA/Genotipi_latest/"+pasma
#File with a list of 800 SNPs for parentage verification
SNP800="/home/janao/Genotipi/ParentalVerification_SNPSNP/Names_800SNPs.txt"
#path to Zanardi
ZanDir="/home/janao/Genotipi/Genotipi_CODES/Zanardi"
#file with IDs and seq for the animals
RJ_IDSeq="/home/janao/Genotipi/Genotipi_CODES/Rjave_seq_ID.csv"
#SNP coding
SNPSifrant="/home/janao/Genotipi/ParentalVerification_SNPSNP/Sifrant_SNP.csv"

#########################################################
##########################################################
#create a directory with the current date for temp genotype manipulation
if not os.path.exists(tempDir):
    os.makedirs(tempDir)


#extract downloaded files into created temp directory
os.chdir("/home/janao/Downloads")
zf = zipfile.ZipFile(zip_file, 'r')
print zf.namelist()
zf.extractall(tempDir)
zf.close()

#change current working directory to the created directory
os.chdir(tempDir)
#list all the directories within - all the downloaded "genotype packages"
genPacks = filter(lambda x: x.startswith('Matija_Rigler'), os.listdir(os.getcwd()))


################################################################
#loop through all subdirectories
#################################################################
#extract FinalReport and SNP_Map  into the directory of the genotype package
for genPack in genPacks:
    os.chdir(tempDir+"/"+genPack)
    FinalFile = filter(lambda x: x.endswith('FinalReport.zip'), os.listdir(os.getcwd()))
    FinalReport = zipfile.ZipFile(FinalFile[0], 'r')
    FinalReport.extractall(os.getcwd())
    FinalTxt=filter(lambda x: x.endswith('FinalReport.txt'), os.listdir(os.getcwd()))[0]
    SERNUM=genPack[-9:].strip("_") # take last nine, since some date are in 22apr2016 form
    SNP_Map = zipfile.ZipFile("SNP_Map.zip", 'r')
    SNP_Map.extractall(os.getcwd())
    
    #copy pedda.param and python script to the current directory
    shutil.copy((peddarow+"/peddar.param"), "peddar.param")
    shutil.copy((peddarow+"/pedda_row.py"), "pedda_row.py")
    #replace strings with shell command
    os.system('sed -i "s%SI %SI%g" ' + FinalTxt) #remove space in sample IDs
    os.system('sed -i "s%  % %g" ' + FinalTxt) #remove double spacing
    os.system('sed -i "s/test_FinalReport.txt/"'+FinalTxt+'"/g" peddar.param') #insert FinalReport name into peddar.param
    os.system('sed -i "s/test_outputfile/"'+genPack+'"/g" peddar.param') #insert OutPut name into peddar.param
    os.system('sed -i "s/AlleleFormat/"'+AlleleFormat+'"/g" peddar.param') #insert desired AlleleFormat name into peddar.param
    os.system('sed -i "s/TEST/"'+pasma+'"/g" peddar.param')
    os.system("python pedda_row.py") #transform into ped and map file
    
    #read in the ped file, check whether there are names in the third column
    pedname=filter(lambda x: x.endswith(".ped"), os.listdir(os.getcwd()))[0] #get the name of the ped file
    pedfile=open(pedname).read().strip("\n").split("\n") #read the file, strip the last new line and split by new lines
    names = []
    ids=[]
    for line in pedfile:
        ids.append(line.split(" ")[1])
        if line.split(" ")[2] != "0": #split by space and read the third column
            print "Names in the ped file!"
            names.append(line.split(" ")[2]) #create a list of animal names in the ped file
    if names: #if there are any names found, remove them / else don't do anything
        for name in names:
            os.system('sed -i "s/'+name+' //g" '+pedname) #remove all the names from the ped file
        print "Names removed from the ped file"
    
    mapname=filter(lambda x: x.endswith(".map"), os.listdir(os.getcwd()))[0] #get the name of the ped file
    mapfile=open(mapname).read().strip("\n").split("\n") #read the file, strip the last new line and split by new lines
    chip=chips[len(mapfile)] # get the chip name according to the number of SNPs
    GenoFile[chip][genPack] = ids #add animals ids under chip - genotyping package name
    DateDownloaded[genPack] = date #date of the genotyping package download
    PedFiles[chip].append(pedname)
    MapFiles[chip].append(mapname)
    
    #extract the 800 SNPs
    command = "plink --file {0} --cow --extract {1} --recode --out {2}".format(genPack, SNP800, genPack+"_SNP800")
    os.system(command)
    SNP800_Peds.append(os.getcwd()+"/"+genPack+"_SNP800.ped")
    SNP800_Maps.append(os.getcwd()+"/"+genPack+"_SNP800.map")      
                  
    #copy ped and map file to the Genotipi_latest
    shutil.copy(pedname, Gen_lat+"/"+chip)
    shutil.copy(mapname, Gen_lat+"/"+chip)     

    
        
################################################################################################
###############################################################################################
#END OF THE LOOP
#merge produced SNP800 files
#merge ped files if merge_ask = Y
#create table for govedo
#############################################################################################
###############################################################################################
    
#CHANGE THE DIRECTORY back to the tempDir
os.chdir(tempDir)

#merge SNP800 files
#first check whether there are more than one SNP800 files
if len(SNP800_Peds) == len(SNP800_Maps) & len(SNP800_Peds) != 1:
    shutil.copy("/home/janao/Genotipi/Genotipi_CODES/PARAMFILE.txt", Gen_lat+"/"+chip) #copy the PARAMFILE 
    # write a --merge-list file for plink with all the SNP800 file on the download date
    with open('SNP800_List.txt', 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter=" ")
        [writer.writerow(r) for r in zip(SNP800_Peds[1:], SNP800_Maps[1:])] #leave the first one out - that goes in the plink command line
    #merge with plink
    mergecommand = "plink --file {0} --cow --merge-list {1} --recode --out {2}".format((SNP800_Peds[0].strip(".ped")),"SNP800_List.txt", str(date)+"_SNP800") 
    os.system(mergecommand)
elif len(SNP800_Peds)==1 & len(SNP800_Maps) == 1: #if there is only one genoPackage, then just copy the SNP800 file and rename it to date_SNP800
    print "Only one SNP800 file"
    shutil.move(SNP800_Peds[0], tempDir+"/"+str(date)+"_SNP800.ped")
    shutil.move(SNP800_Maps[0], tempDir+"/"+str(date)+"_SNP800.map")
elif len(SNP800_Peds) != len(SNP800_Maps): #if the number of SNP800 ped and map files is not the same, report an error
    print "The length of Ped and Map SNP800 List are not equal, check for inconsistences and errors"



#employing outer R script to prepare Govedo table for SNPs
shutil.copy("/home/janao/Genotipi/Genotipi_CODES/Add_alleles_Govedo_PARENTAL_SNP800_python.R", os.getcwd())
os.system('sed -i "s|SNP800genoFile|"'+str(date)+'_SNP800.ped"|g" Add_alleles_Govedo_PARENTAL_SNP800_python.R')
os.system('sed -i "s|SNP800mapFile|"'+str(date)+'_SNP800.map"|g" Add_alleles_Govedo_PARENTAL_SNP800_python.R')
os.system('sed -i "s|SifrantFile|"'+SNPSifrant+'"|g" Add_alleles_Govedo_PARENTAL_SNP800_python.R')
os.system('sed -i "s|Govedo800SNPFile|"'+str(date)+'_GovedoSNP800.csv"|g" Add_alleles_Govedo_PARENTAL_SNP800_python.R')
os.system('sed -i "s|IDSeqFile|"'+RJ_IDSeq+'"|g" Add_alleles_Govedo_PARENTAL_SNP800_python.R')
os.system("Rscript Add_alleles_Govedo_PARENTAL_SNP800_python.R")



#read in animal ID / Seq / DateOfBirth / SexCode table
#create a dictionary
Rj_IDSeq_Dict = defaultdict()
with open(RJ_IDSeq, 'rb') as IDSeq:
    reader = csv.reader(IDSeq, delimiter=',')
    for line in reader:
        Rj_IDSeq_Dict[line[0]] = line[1:]
        
        
#create a table of individuals for govedo     
numInd = 0 #count ind
IndList = []
for chip in GenoFile:
    for genpackage in GenoFile[chip]:
        numInd += len(GenoFile[chip][genpackage])
        for ind in GenoFile[chip][genpackage]:
            IndList.append((Rj_IDSeq_Dict[ind][0], chip, genpackage[-9:].strip("_"), date))#create a list (seq, chip, dategenotyped, datedownloaded)
            
#write the table with individuals to the tempDit      
with open(tempDir+"/" +str(date)+'_IndGovedo.csv', 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter=" ")
        [writer.writerow(r) for r in IndList]

 
#check whether the number of individuals in all genoPackages corresponds to the number of individuals with extracted 800SNPs for parentage testing
if numInd == len():
    print "The number of individuals corresponds to the number of individuals in SNP800 ped file"
else: 
    print "The number of individuals DOES NOT correspond to the number of individuals in SNP800 ped file, check for ERRORS!"






""" TRYING TO CREATE GOVEDO TABLE - CONCLUDED IT IS MUCH EASIER WITH R
#read in the SNP800 file
SNP800Ped=open(str(date)+"_SNP800.ped").read().strip("\n").split("\n")    
SNP800Map=open(str(date)+"_SNP800.map").read().strip("\n").split("\n")
SNPAlleles = []
for line in SNP800Map:
    SNPAlleles.append(line.split("\t")[1])
    SNPAlleles.append(line.split("\t")[1])






#read in SNPSifrant and create a dictionary: {SNP: set(code1, code2)}
SNPSifrant_Dict = defaultdict(set)
with open(SNPSifrant, 'rb') as SNP_Sifrant:
    reader = csv.reader(SNP_Sifrant, delimiter=',')
    for line in reader:
        SNPSifrant_Dict[line[0]].add(line[1])

for line in SNP800Ped:
    alleles=line.split(" ")[6:]
tableSNP800 = zip(
"""
#merge is outside the loop
#merge
if merge_ask == 'Y':
    shutil.copy("/home/janao/Genotipi/Genotipi_CODES/PARAMFILE.txt", Gen_lat+"/"+chip)
    pedToMerge = ",".join(SNP800_Peds).strip("'")
    os.system('sed -i "s/PathToped/"'+pedToMerge+'"/g" PARAMFILE.txt')
#AT THE END!
#prepare a table with individuals for Govedo

    
    
    

        
#PARAMFILE for merging
#Table -genotyped ind for govedo
#SNP800 extraction and table
#at the end print a report - number of genotypes animals, number with sequences ...