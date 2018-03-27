# -*- coding: utf-8 -*-
import os
import pandas as pd
from collections import defaultdict

os.chdir("/home/jana/Documents/PhD/CompBio/")
herds = pd.read_table("/home/jana/Documents/PhD/CompBio/TestingGBLUP/PedCows_HERDS.txt", sep=" ")
IndGeno = pd.read_table("/home/jana/Documents/PhD/CompBio/IndForGeno_5gen.txt", header=None)

RefAmean = defaultdict()


number = 1

for herd1 in range(1, 101):
        for herd2 in range(herd1,101):

            ref = sorted(list(herds.Indiv[herds.cluster.isin([herd1, herd2])])) #tukaj odberi živali v obeh čredah

            pd.DataFrame({"ID": ref}).to_csv("/home/jana/Documents/PhD/CompBio/IndMatrix.txt", index=None, header=None)


            os.system("grep -Fwf IndMatrix.txt PedigreeNrm.txt  > RefMatrix")
            a = pd.read_table("/home/jana/Documents/PhD/CompBio/RefMatrix", sep="\s+", header=None)
            a.columns =  ["Indiv"] + list(IndGeno.loc[:,0])


            refA = a.loc[:, ref]


            meanRef = mean(refA).mean()

            RefAmean[number] = [herd1, herd2, meanRef]
            number = number + 1


RefDF = pd.DataFrame.from_dict(RefAmean, orient="index")
RefADF = RefDF.drop_duplicates()
RefADF.columns = ["Herd1", "Herd2", "A"]
RefADF.to_csv("RefADF_mean.csv", index=None)


ped = pd.read_table("/home/jana/Documents/PhD/CompBio/PedigreeAndGeneticValues_cat.txt", sep=" ")
nr = ped.Indiv[ped.cat.isin(['potomciNP'])]
pb = ped.Indiv[ped.cat == 'pb']

NapAmean = defaultdict()
PbAmean = defaultdict()
number = 1         
for herd in range(1,101):
            ref = sorted(list(herds.Indiv[herds.cluster == herd])) #tukaj odberi živali v obeh čredah

            pd.DataFrame({"ID": ref}).to_csv("/home/jana/Documents/PhD/CompBio/IndHerd.txt", index=None, header=None)


            os.system("grep -Fwf IndHerd.txt PedigreeNrm.txt  > HerdMatrix")
            a = pd.read_table("/home/jana/Documents/PhD/CompBio/HerdMatrix", sep="\s+", header=None)
            a.columns =  ["Indiv"] + list(IndGeno.loc[:,0])

            refnapA = a.loc[:, list(nr)] # sorodstvo z napovedno populacijo
            refpbA = a.loc[:, list(pb)] # orodstvo s plemenskimi biki
            meanRefNap = mean(refnapA).mean()
            meanRefPb = mean(refpbA).mean()

            NapAmean[number] = [herd, meanRefNap]
            PbAmean[number] = [herd, meanRefPb]
            number = number + 1
            
            
            
NapADF = pd.DataFrame.from_dict(NapAmean, orient="index")
NapADF.columns = ["Herd","A"]
NapADF.to_csv("NapADF_mean.csv", index=None)


PbADF = pd.DataFrame.from_dict(PbAmean, orient="index")
PbADF.columns = ["Herd","A"]
PbADF.to_csv("PbADF_mean.csv", index=None)