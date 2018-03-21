# -*- coding: utf-8 -*-
import os
import pandas as pd
from collections import defaultdict

os.chdir("/home/jana/Documents/PhD/CompBio/")
herds = pd.read_table("/home/jana/Documents/PhD/CompBio/TestingGBLUP/PedCows_HERDS.txt", sep=" ")
IndGeno = pd.read_table("/home/jana/Documents/PhD/CompBio/IndForGeno_5gen.txt", header=None)

RefAsum = defaultdict()


number = 1

for herd1 in range(1, 3):
        for herd2 in range(herd1,3):

            ref = sorted(list(herds.Indiv[herds.cluster.isin([herd1, herd2])])) #tukaj odberi živali v obeh čredah

            pd.DataFrame({"ID": ref}).to_csv("/home/jana/Documents/PhD/CompBio/IndMatrix.txt", index=None, header=None)


            os.system("grep -Fwf IndMatrix.txt PedigreeNrm.txt  > RefMatrix")
            a = pd.read_table("/home/jana/Documents/PhD/CompBio/RefMatrix", sep="\s+", header=None)
            a.columns =  ["Indiv"] + list(IndGeno.loc[:,0])


            refA = a.loc[:, ref]


            sumRef = sum(refA).sum()

            RefAsum[number] = [herd1, herd2, sumRef]
            number = number + 1


RefDF = pd.DataFrame.from_dict(RefAsum, orient="index")
RefADF = RefDF.drop_duplicates()
RefADF.to_csv("Ref)


NapAsum = defaultdict()
PbAsum = defaultdict()
number = 1         
for herd in range(1,101):
            ref = sorted(list(herds.Indiv[herds.cluster == herd])) #tukaj odberi živali v obeh čredah

            pd.DataFrame({"ID": ref}).to_csv("/home/jana/Documents/PhD/CompBio/IndHerd.txt", index=None, header=None)


            os.system("grep -Fwf IndHerd.txt PedigreeNrm.txt  > HerdMatrix")
            a = pd.read_table("/home/jana/Documents/PhD/CompBio/HerdMatrix", sep="\s+", header=None)
            a.columns =  ["Indiv"] + list(IndGeno.loc[:,0])

            refnapA = a.loc[:, list(nr)] # sorodstvo z napovedno populacijo
            refpbA = a.loc[:, list(pb)] # orodstvo s plemenskimi biki
            sumRefNap = sum(refnapA).sum()
            sumRefPb = sum(refpbA).sum()

            NapAsum[number] = [herd, sumRefNap]
            PbAsum[number] = [herd, sumRefPb]
            number = number + 1