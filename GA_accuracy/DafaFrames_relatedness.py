# -*- coding: utf-8 -*-
import os
import pandas as pd
from collections import defaultdict
import numpy as np

herds = pd.read_table("PedCows_HERDS.txt", sep=" ")
IndGeno = pd.read_table("IndForGeno.txt", header=None)

RefAmean = defaultdict()


number = 1

for herd1 in range(1, 101):
        for herd2 in range(herd1,101):

            ref = sorted(list(herds.Indiv[herds.cluster.isin([herd1, herd2])])) #tukaj odberi živali v obeh čredah

            pd.DataFrame({"ID": ref}).to_csv("IndMatrix.txt", index=None, header=None)


            os.system("grep -Fwf IndMatrix.txt PedigreeNrm.txt.txt  > RefMatrix")
            a = pd.read_table("RefMatrix", sep="\s+", header=None)
            a.columns =  ["Indiv"] + list(IndGeno.loc[:,0])


            refA = a.loc[:, ref]


            meanRef = np.mean(refA).mean()

            RefAmean[number] = [herd1, herd2, meanRef]
            number = number + 1


RefDF = pd.DataFrame.from_dict(RefAmean, orient="index")
RefADF = RefDF.drop_duplicates()
RefADF.columns = ["Herd1", "Herd2", "A"]
RefADF.to_csv("RefADF_mean.csv", index=None)


ped = pd.read_table("PedigreeAndGeneticValues_cat.txt", sep=" ")
nr = ped.Indiv[ped.cat.isin(['potomciNP'])]
pb = ped.Indiv[ped.cat == 'pb']

NapAmean = defaultdict()
PbAmean = defaultdict()
number = 1         
for herd in range(1,101):
            ref = sorted(list(herds.Indiv[herds.cluster == herd])) #tukaj odberi živali v obeh čredah

            pd.DataFrame({"ID": ref}).to_csv("IndHerd.txt", index=None, header=None)


            os.system("grep -Fwf IndHerd.txt Hmatrix.txt  > HerdMatrix")
            a = pd.read_table("HerdMatrix", sep="\s+", header=None)
            a.columns =  ["Indiv"] + list(IndGeno.loc[:,0])

            refnapA = a.loc[:, list(nr)] # sorodstvo z napovedno populacijo
            refpbA = a.loc[:, list(pb)] # orodstvo s plemenskimi biki
            meanRefNap = np.mean(refnapA).mean()
            meanRefPb = np.mean(refpbA).mean()

            NapAmean[number] = [herd, meanRefNap]
            PbAmean[number] = [herd, meanRefPb]
            number = number + 1
            
            
            
NapADF = pd.DataFrame.from_dict(NapAmean, orient="index")
NapADF.columns = ["Herd","A"]
NapADF.to_csv("NapADF_mean.csv", index=None)


PbADF = pd.DataFrame.from_dict(PbAmean, orient="index")
PbADF.columns = ["Herd","A"]
PbADF.to_csv("PbADF_mean.csv", index=None)