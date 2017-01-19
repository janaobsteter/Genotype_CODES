import os
import zipfile
import commands
"""
genDir="/home/janao/bsw/"
genFiles=os.listdir(genDir)
SampleDir="/home/janao/SampleJurij"


for genFile in genFiles:
    os.chdir(genDir+"/"+genFile)
    MatijaZip = filter(lambda x: x.endswith('Sample_Map.zip'), os.listdir(os.getcwd()))
    if MatijaZip:
        MatijaZip = zipfile.ZipFile(MatijaZip[0], 'r')
        MatijaZip.extractall(SampleDir)
        os.rename(SampleDir+"/Sample_Map.txt", SampleDir+"/Sample_Map"+genFile+".txt")
    if not MatijaZip:
        MatijaZipM = filter(lambda x: x.endswith('.zip'), os.listdir(os.getcwd()))
        MatijaZipM = zipfile.ZipFile(MatijaZipM[0], 'r')
        MatijaZipM.extractall(os.getcwd())
        subDirs=filter(os.path.isdir, os.listdir(os.getcwd()))
        for subdir in subDirs:
            os.chdir(genDir+"/"+genFile)
            MatijaZip = filter(lambda x: x.endswith('Sample_Map.zip'), os.listdir(os.getcwd()))
            if MatijaZip:
                    MatijaZip = zipfile.ZipFile(MatijaZip[0], 'r')
                    MatijaZip.extractall(SampleDir)
                    os.rename(SampleDir+"/Sample_Map.txt", SampleDir+"/Sample_Map"+subdir+".txt") 
                    
 """                     
                        
genDir='/home/janao/GENZIP/bsw/'

os.chdir(genDir)
subDirs=filter(os.path.isdir, os.listdir(genDir))
GenZipDir='/home/janao/GENZIP/ZipFiles/'

import shutil
for subDir in subDirs:
    print subDir
    os.chdir(genDir+subDir)
    if len(filter(lambda x: x.endswith(".zip"), os.listdir(os.getcwd()))) >= 6: #check how many files there are in the directory
        name=filter(lambda x: x.endswith('FinalReport.zip'), os.listdir(os.getcwd()))[0].strip("_FinalReport.zip") #find the finalreport to create a name
        if not os.path.isfile(name + ".zip"): #check whether the files are already ziped and zip only if they arent'
            status, output = commands.getstatusoutput("zip -r " + name + " *")  #check also the status
            if status != 0: #is it not 0, then pass
                print "Unable to zip files in " + subDir 
                pass
        if not os.path.isfile(GenZipDir + name + ".zip"): #check whether you have already moved the file
            shutil.move(genDir+subDir+"/"+name+".zip", GenZipDir) #else move
    elif len(filter(lambda x: x.endswith(".zip"), os.listdir(os.getcwd()))) == 1:
        MatijaZip = filter(lambda x: x.endswith('.zip'), os.listdir(os.getcwd()))[0]
        status, output = commands.getstatusoutput("unzip " + MatijaZip)
        if status != 0:
            print "Unable to unzip " + subDir +"/"+ MatijaZip
            pass
        #remove blank spaces from file names
        status, output = commands.getstatusoutput("bash /home/janao/GENZIP/rename.sh")
        if status != 0:
            print "Unable to perform 'bash /home/janao/GENZIP/rename.sh' command"
            exit()
        withDirs=filter(os.path.isdir, os.listdir(os.getcwd()))
        for withDir in withDirs:
            os.chdir(genDir + subDir + "/" + withDir)
            if len(filter(lambda x: x.endswith('zip'), os.listdir(os.getcwd()))) == 6:
                name=filter(lambda x: x.endswith('FinalReport.zip'), os.listdir(os.getcwd()))[0].strip("_FinalReport.zip")
                if not os.path.isfile(name + ".zip"):
                    status, output = commands.getstatusoutput("zip -r " + name + " *")
                    if status != 0:
                        print "Unable to zip files in " + subDir + withDir 
                        pass
                if not os.path.isfile(GenZipDir + name + ".zip"):
                    shutil.move(genDir+subDir+"/"+withDir+"/"+name+".zip", GenZipDir)
    else: 
        print "Some error occured, suspicious number of files in the directory "+ subDir
            
        