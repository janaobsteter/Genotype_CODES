#script to perform initial extraction of all up-to-date FinalReport package from Jurij bsw folder
#it zips all the files for one genotype package and moves it to a new ZipFile directory

import os
import zipfile
import commands
                  
                        
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
    elif len(filter(lambda x: x.endswith(".zip"), os.listdir(os.getcwd()))) == 1: #if there is only one zip file --> contains many packages, unzip it
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
        withDirs=filter(os.path.isdir, os.listdir(os.getcwd())) #now repeat the zip and move for each directory created with unzip
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
    else: #if a strange number of files in a directory
        print "Some error occured, suspicious number of files in the directory "+ subDir
            
        