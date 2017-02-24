import shutil
import GenFiles
import os
import zipfile
import tempfile
import pandas as pd
from collections import defaultdict
import csv


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


Zip_lat='/run/user/1000/gvfs/smb-share:server=kis-h2.si,share=kisdfs/ZIV/vol1/ZIV/VSI/JanaO/Rjava/ZipFiles/'
tempDir='/home/janao/Genotipi/SampleMaps/'
zipPackages = (filter(lambda x: x.endswith('.zip'), os.listdir(Zip_lat)))
os.chdir(tempDir)

#correct sample names and rezip files
for zipPackage in zipPackages:
    shutil.copy(Zip_lat+zipPackage, tempDir) #copy zipPackage into temp dir
    onePackage=GenFiles.genZipPackage(zipPackage)
    errorIDs = onePackage.extractErrorNames() #extract Sample Names if they exist - they shouldnt be in the file
    if errorIDs:
        shutil.move(onePackage.name+'_Sample_Map.txt', 'Sample_Map.txt') #rename extracted SampleMap
        onePackage.extractFinalReport() #extract Finalreport to replace the spurious names
        for i in errorIDs:
            os.system('sed -i  "s|' +str(i[0])+ '|' + i[1] + '|g" ' + onePackage.name+"_FinalReport.txt") #errorIDs are tuples, replace first element witht the second
            os.system('sed -i  "s|' +str(i[0])+ '|' + i[1] + '|g" '+'Sample_Map.txt') 
        #remove old Sample_Map and FinalReport from the zip archive and put the new one into the archive
        #zip_deflated to compress the zip (reduce in size)
        with zipfile.ZipFile(onePackage.name+'_FinalReport.zip', 'w', zipfile.ZIP_DEFLATED) as myzip:
            myzip.write(onePackage.name+'_FinalReport.txt') #create new FinalReport zip  
        with zipfile.ZipFile('Sample_Map.zip', 'w', zipfile.ZIP_DEFLATED) as myzip:
            myzip.write('Sample_Map.txt')      #create new Sample_Map.zip

        remove_from_zip(onePackage.zipname, onePackage.name+'_FinalReport.zip')
        remove_from_zip(onePackage.zipname, 'Sample_Map.zip')
        
        with zipfile.ZipFile(onePackage.zipname, 'a', zipfile.ZIP_DEFLATED) as z:
            z.write(onePackage.name+'_FinalReport.zip')
        with zipfile.ZipFile(onePackage.zipname, 'a', zipfile.ZIP_DEFLATED) as z:
            z.write('Sample_Map.zip')     
        os.remove(onePackage.name+'_FinalReport.zip') #remove temp extracted FInalReports (with wrong names)
        os.remove(onePackage.name+'_FinalReport.txt')
        os.remove('Sample_Map.txt')
        os.remove('Sample_Map.zip')
    else: #if no errorIDs are found, remove the zip package thath you've just copied
        os.remove(zipPackage) 
        os.remove(onePackage.name+'_Sample_Map.txt')
 
######################################################################       
######################################################################       
#this is to check whether all the names in the files got corrected right     
tempDir='/home/janao/Genotipi/SampleMaps/'
os.chdir(tempDir)
AllGenInd = []
zipPackages = (filter(lambda x: x.endswith('.zip'), os.listdir(tempDir)))

for zipPackage in zipPackages:
    #shutil.copy(Zip_lat+zipPackage, tempDir)
    onePackage=GenFiles.genZipPackage(zipPackage)
    onePackage.extractSampleMap()
    table=pd.read_table(onePackage.name + '_Sample_Map.txt')
    AllGenInd += [(onePackage.name,x) for x in table['ID']]
    os.remove(onePackage.name + '_Sample_Map.txt')


##################################################################################3
##################################################################################
#check how many IDs not found in pedigree (sequences)
RJ_IDSeq="/home/janao/Genotipi/Genotipi_CODES/Rjave_seq_ID.csv"
Rj_IDSeq_Dict = defaultdict()
with open(RJ_IDSeq, 'rb') as IDSeq:
    reader = csv.reader(IDSeq, delimiter=',')
    for line in reader:
        Rj_IDSeq_Dict[line[0]] = line[1:]
         
lala = AllGenInd


    
    
AllGenIDs = [x[1] for x in AllGenInd]
AllGenIDs = pd.DataFrame(AllGenIDs)
AllGenIDs.to_csv(tempDir+'/AllGenInd.csv', sep=",")
All = pd.read_csv('AllGenInd.csv')
AllGenID = [] #vsi iz novih RJ Finalreports
errorIDs = []
for ind in [x.upper() for x in All['0']]:
    try:
        AllGenID.append((Rj_IDSeq_Dict.get(ind)[0]))
    except:
        errorIDs.append(ind)
 
 ############################################################################
 ############################################################################     
#the errouneous IDs inserted by GeneSeek        
replaceIDs = [('SI4574059','SI04574059'),('SI84048801','SI84048802'),('SI4384195','SI04384195'),('Si24289407','SI24289407')]
spPackages =[]
for i in AllGenInd:
    for errorID in errorIDs:
        if errorID in i:
            spPackages.append(i)
spPackages=list(set(spPackages))
zipErrorPackages=[i[0]+'.zip' for i in spPackages]         
zipErrorPackages=['Matija_Rigler_BOVGP4V01-2_20160926-2.zip'] #not all capitals
#replace only in this

for zipPackage in zipErrorPackages:
    shutil.copy(Zip_lat+zipPackage, tempDir) #copy zipPackage into temp dir
    onePackage=GenFiles.genZipPackage(zipPackage)
    errorIDs = onePackage.extractErrorNames() #extract Sample Names if they exist - they shouldnt be in the file
    if errorIDs:
        shutil.move(onePackage.name+'_Sample_Map.txt', 'Sample_Map.txt') #rename extracted SampleMap
        onePackage.extractFinalReport() #extract Finalreport to replace the spurious names
        for i in errorIDs:
            os.system('sed -i  "s|' +str(i[0])+ '|' + i[1] + '|g" ' + onePackage.name+"_FinalReport.txt") #errorIDs are tuples, replace first element witht the second
            os.system('sed -i  "s|' +str(i[0])+ '|' + i[1] + '|g" '+'Sample_Map.txt') 
        #replace spurious IDs
        for i in replaceIDs:
            os.system('sed -i  "s|' +i[0]+ '|' + i[1] + '|g" ' + onePackage.name+"_FinalReport.txt") #errorIDs are tuples, replace first element witht the second
            os.system('sed -i  "s|' +i[0]+ '|' + i[1] + '|g" '+'Sample_Map.txt') 
        
        #remove old Sample_Map and FinalReport from the zip archive and put the new one into the archive
        #zip_deflated to compress the zip (reduce in size)
        with zipfile.ZipFile(onePackage.name+'_FinalReport.zip', 'w', zipfile.ZIP_DEFLATED) as myzip:
            myzip.write(onePackage.name+'_FinalReport.txt') #create new FinalReport zip  
        with zipfile.ZipFile('Sample_Map.zip', 'w', zipfile.ZIP_DEFLATED) as myzip:
            myzip.write('Sample_Map.txt')      #create new Sample_Map.zip

        remove_from_zip(onePackage.zipname, onePackage.name+'_FinalReport.zip')
        remove_from_zip(onePackage.zipname, 'Sample_Map.zip')
        
        with zipfile.ZipFile(onePackage.zipname, 'a', zipfile.ZIP_DEFLATED) as z:
            z.write(onePackage.name+'_FinalReport.zip')
        with zipfile.ZipFile(onePackage.zipname, 'a', zipfile.ZIP_DEFLATED) as z:
            z.write('Sample_Map.zip')     
        os.remove(onePackage.name+'_FinalReport.zip') #remove temp extracted FInalReports (with wrong names)
        os.remove(onePackage.name+'_FinalReport.txt')
        os.remove('Sample_Map.txt')
        os.remove('Sample_Map.zip')
    else: #if no errorIDs are found, remove the zip package thath you've just copied
        os.remove(zipPackage) 
        os.remove(onePackage.name+'_Sample_Map.txt')
#after this check again for the errorIDs
                                                                            
SQLGen = pd.read_csv('/home/janao/Genotipi/AllGen_27012017.csv')
SQlGenInd = [str(x) for x in SQLGen['ZIV_ID_SEQ']] #vsi iz sqla

RealAllInd = list(set(SQlGenInd) | set(AllGenID))
RealAllInd1 = pd.DataFrame({'ZIV_ID_SEQ':RealAllInd})
RealAllInd1.to_csv('AllIndSet_27012017.csv', sep=" ", index=None)