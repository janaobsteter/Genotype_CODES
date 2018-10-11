import os
import GenFiles
import pandas as pd
from collections import defaultdict
import subprocess
import re
from itertools import chain

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

pedChip = defaultdict()
for dir in dirs:
    peds = [x.strip(".ped") for x in os.listdir(workdir + "/" + dir + "/") if x.endswith(".ped") and "Clean" not in x]
    for ped in peds:
        try:
            if GenFiles.mapFile(workdir + "/" + dir + "/" + ped + ".map").chip in chips.values():
                pedChip[workdir + "/" + dir + "/" + ped] = GenFiles.mapFile(workdir + "/" + dir + "/" + ped + ".map").chip
        except:
            pass


print(pedChip)

mergedir = "/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/Top/"

compareDF = pd.DataFrame(columns=["File", "Concordance", "Format", "NumIndiv"])
comparedir = "/home/jana/Genotipi/Genotipi_DATA/Compare/"
os.chdir(comparedir)
for ped in pedChip.keys():
    print(ped)
    os.system("cut " + ped + '.ped -f1,2 -d" " > Inds.txt')
    numindiv = os.popen("less Inds.txt | wc -l").read().strip()
    os.system("cut " + ped + '.ped -f1500 -d" " |  sort | uniq > Alleles.txt')
    alleles = open("Alleles.txt").read().strip().split("\n")
    print(alleles)
    if 'G' in alleles:
        format = "Top"
    if 'B' in alleles:
        format = "AB"
    if alleles == ['A']:
        print(ped + ": only A in alleles.")
        pass
    print("Format: " + format)
    mergedir = "/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/" + format + "/"
    try:
        os.system("plink --file " + ped + " --merge " +
                  mergedir + pedChip[ped] + "/" +  "PLINK_MERGED.ped " +
                  mergedir + pedChip[ped] + "/" +  "PLINK_MERGED.map --cow --recode --out DIFF0")
        print("plink --file " + ped + " --merge " +
                  mergedir + pedChip[ped] + "/" +  "PLINK_MERGED.ped " +
                  mergedir + pedChip[ped] + "/" +  "PLINK_MERGED.map --cow --recode --out " + ped.split("/")[-1] + "DIFF0")

        out = open(comparedir + "DIFF0.log").read().strip().split("\n")
        print(out)
        c = [x for x in chain.from_iterable([x.split(" ") for x in out if "Variant" in x]) if
             bool(re.search(r'\d', x))]
        print(c)
        print(os.getcwd())
        pd.DataFrame({"Variant": c}).to_csv(comparedir + "SpuriousSNPs.txt", header=None, index=None)
        os.system("""sed -i "s/'//g" """ + comparedir + "SpuriousSNPs.txt")
        os.system("grep -Fwf " + comparedir + "SpuriousSNPs.txt " + ped + ".map > " + comparedir + "RemoveSNPs.txt")

        os.system("plink --file " + ped + " --merge " +
                  mergedir + pedChip[ped] + "/" + "PLINK_MERGED.ped " +
                  mergedir + pedChip[
                      ped] + "/" + "PLINK_MERGED.map --exclude RemoveSNPs.txt --merge-mode 7 --cow "
                                   "--keep Inds.txt --recode --out DIFF > DIFFtmp.txt")

        a = open("DIFFtmp.txt").read().strip().split("\n")
        c = [x for x in a if "concordance rate" in x][0].split(" ")[-1].strip(".")
        compareDF = compareDF.append(pd.DataFrame({"File": [ped.split("/")[-1]], "Concordance": [c], "Format": [format],
                                                   "NumIndiv": numindiv}))
    except:
        pass

compareDF.to_csv(comparedir + "CompareDF.csv")