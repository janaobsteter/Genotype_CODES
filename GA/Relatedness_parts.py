# -*- coding: utf-8 -*-
from __future__ import division
import pandas as pd
import numpy as np
from collections import defaultdict
import sys
import os
from itertools import combinations_with_replacement


start = int(sys.argv[1])


#Calculate relatedness according to herds
herds = pd.read_table("PedCows_HERDS.txt", sep=" ")
IndGeno = pd.read_table("INDPED.txt", header=None)


#Tukaj izračunaj sorodstvo med živalmi v obema čredama
RefAmean = defaultdict()

number = 1

comb = [x for x in combinations_with_replacement(range(1, 101), 2)]

for pair in comb[start:(start+505)]:
    herd1, herd2 = pair
    print(str(herd1) + "_" + str(herd2))
    ref = sorted(list(herds.Indiv[herds.cluster.isin([herd1, herd2])]))  # tukaj odberi živali v obeh čredah

    pd.DataFrame({"ID": ref}).to_csv("IndMatrix" + str(herd1) + str(herd2) + ".txt", index=None, header=None)

    os.system("grep -Fwf IndMatrix" + str(herd1) + str(herd2) + ".txt PedigreeNrm.txt > RefMatrix" + str(herd1) + str(herd2))
    a = pd.read_table("RefMatrix" + str(herd1) + str(herd2), sep="\s+", header=None)
    a.columns = ["Indiv"] + list(IndGeno.loc[:, 0])

    refA = a.loc[:, ref]
    meanRef = np.mean(refA).mean()
    RefAmean[number] = [herd1, herd2, meanRef]


    os.system("rm RefMatrix" + str(herd1) + str(herd2))
    os.system("rm IndMatrix" + str(herd1) + str(herd2) + ".txt")
    number = number + 1

RefDF = pd.DataFrame.from_dict(RefAmean, orient="index")
RefADF = RefDF.drop_duplicates()
RefADF.columns = ["Herd1", "Herd2", "A"]
RefADF.to_csv("RefADF_mean" + str(start) + ".csv", index=None)
