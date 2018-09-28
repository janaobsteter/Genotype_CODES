import os
import GenFiles
import pandas as pd
from collections import defaultdict

workdir = "/home/jana/Genotipi/Genotipi_DATA/Rjava_TEMP/"
os.chdir("/home/jana/Genotipi/Genotipi_DATA/Rjava_TEMP/")

os.system("ls -d Genotipi*/ > Dirs.txt")

dirs = list(pd.read_table("Dirs.txt", header=None).loc[:,0])
print(dirs)


chips = {19720: "GGPv02",
26145: "GGPv03",
26151: "GGPv03",
30105: "GGPv04",
30106: "GGPv04",
76883: "HD" ,
138892: "HDv02",
139376: "HDv02",
54001:"50Kv01" ,
54609: "50Kv02",
51274: "IDBv03",
52445: "IDBv03"
         }

print(chips.values())

pedPars = []
for dir in dirs:
    peds = [x.strip(".ped") for x in os.listdir(workdir + "/" + dir + "/") if x.endswith(".ped") and "Clean" not in x]
    for ped in peds:
        try:
            chip = GenFiles.mapFile(workdir + "/" + dir + "/" + ped + ".map").chip
            if chip in chips.values() and os.path.isfile(workdir + "/" + dir + "/" + ped + "_" + chip + ".ped"):
                print(workdir + "/" + dir + "/" + ped + "_" + chip + ".ped")
                pedPars.append((workdir + "/" + dir + "/" + ped,
                                workdir + "/" + dir + "/" + ped + "_" + chip))
        except:
            pass


print(pedPars)

compareDF = pd.DataFrame(columns=["File", "Concordance", "Format"])
comparedir = "/home/jana/Genotipi/Genotipi_DATA/Compare/"
os.chdir(comparedir)
for ped in pedPars:
    (ped1, ped2) = ped
    os.system("cut " + ped1 + '.ped -f1500 -d" " |  sort | uniq > Alleles.txt')
    alleles = open("Alleles.txt").read().strip().split("\n")
    if 'G' in alleles:
        format = "Top"
    if 'B' in alleles:
        format = "AB"
    if alleles == ['A']:
        print(ped + ": only A in alleles.")
        pass
    print("Format: " + format)
    try:
        os.system("plink --file " + ped1 + " --merge " +
                  ped2  + ".ped " + ped2 + ".map --merge-mode 7 --cow --recode --out DIFF > DIFFtmp.txt")
        a = open("DIFFtmp.txt").read().strip().split("\n")
        c = [x for x in a if "concordance rate" in x][0].split(" ")[-1].strip(".")
        compareDF = compareDF.append(pd.DataFrame({"File": [ped1.split("/")[-1]], "Concordance": [c], "Format": [format]}))
    except:
        pass

compareDF.to_csv(comparedir + "ComparePlink1PLink2DF.csv")