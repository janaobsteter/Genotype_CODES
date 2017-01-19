import os
ZipDir = "/home/janao/GENZIP/ZipFiles/"
SampleDir = '/home/janao/GENZIP/SampleMaps/'
Zips = os.listdir(ZipDir)
os.chdir(ZipDir)


for zipFile in Zips:
    name = zipFile.strip(".zip")+"_SampleMap"
    os.system("unzip -j "+zipFile+"  Sample_Map.zip -d .")
    os.system("mv Sample_Map.zip " + SampleDir + name + ".zip")
    os.system("unzip " + SampleDir+name)
    os.system("mv Sample_Map.txt " + name + ".txt")
    
    
os.system("mv *SampleMap.txt " + SampleDir)
    