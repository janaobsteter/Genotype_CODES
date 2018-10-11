import os
import GenFiles
import pandas as pd
from collections import defaultdict
from itertools import chain
import re

workdir = "/home/jana/Genotipi/Genotipi_DATA/Compare/CompareAB/"
os.chdir(workdir)


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

dirs = ["GGPv04", "HDv02"]



comparedir = "/home/jana/Genotipi/Genotipi_DATA/Compare/"
os.chdir(comparedir)
pedPars = pd.DataFrame()
for dir in dirs:
    compareDF = pd.DataFrame(columns=["File", "Concordance", "Reference"])
    os.chdir(workdir + "/" + dir + "/")
    peds = [x.strip(".ped") for x in os.listdir(os.getcwd()) if x.endswith(".ped") and x.startswith("Matija")]
    for ped in peds:
        print(ped)
        os.system('cut -f1,2 -d" " ' + ped + '.ped > Inds.txt')
        for merged in ["PLINK_MERGEDOne", "PLINK_MERGED"]:
            os.system("plink --file " + ped + " --merge " +
                      merged + ".ped " + merged + ".map --cow --recode --out DIFF0 > DIFF0.txt")
            #wrong SNPs
            out = open(merged + ".log").read().strip().split("\n")
            c = [x for x in chain.from_iterable([x.split(" ") for x in out if "Variant" in x]) if
                 bool(re.search(r'\d', x))]
            pd.DataFrame({"Variant": c}).to_csv("SpuriousSNPs.txt", header=None, index=None)
            os.system("""sed -i "s/'//g" SpuriousSNPs.txt""")
            os.system("grep -Fwf SpuriousSNPs.txt " + ped + ".map > RemoveSNPs.txt")
            os.system("plink --file " + ped + " --cow --exclude RemoveSNPs.txt --recode --out " + ped)
            os.system("plink --file " + merged + " --cow --exclude RemoveSNPs.txt --recode --out " + merged)
            os.system("plink --file " + ped + " --merge " +
                      merged + ".ped " + merged + ".map --merge-mode 7 --geno 0 --keep Inds.txt "
                                                    "--cow --recode --out DIFF" + ped + " > DIFFtmp.txt")
            #get concordance
            a = open("DIFFtmp.txt").read().strip().split("\n")
            c = [x for x in a if "concordance rate" in x][0].split(" ")[-1].strip(".")

            compareDF = compareDF.append(pd.DataFrame({"File": [ped.split("/")[-1]], "Concordance": [c], "Reference": [merged]}))


    compareDF.to_csv(workdir + "CompareMerged" + dir + ".csv")