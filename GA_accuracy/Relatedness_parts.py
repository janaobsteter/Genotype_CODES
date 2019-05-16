from __future__ import division

from itertools import combinations_with_replacement

import pandas as pd
import numpy as np
from collections import defaultdict
import sys
import os

start = sys.argv[1]


#Calculate relatedness according to herds
herds = pd.read_table("PedCows_HERDS.txt", sep=" ")
IndGeno = pd.read_table("INDPED.txt", header=None)


#Tukaj izračunaj sorodstvo med živalmi v obema čredama
RefAmean = defaultdict()

number = 1

comb = [x for x in combinations_with_replacement(range(1, 101), 2)]

for pair in comb[start:(start+505)]:
    herd1, herd2 = pair
    ref = sorted(list(herds.Indiv[herds.cluster.isin([herd1, herd2])]))  # tukaj odberi živali v obeh čredah

    pd.DataFrame({"ID": ref}).to_csv("IndMatrix.txt", index=None, header=None)

    os.system("grep -Fwf IndMatrix.txt PedigreeNrm.txt > RefMatrix")
    a = pd.read_table("RefMatrix", sep="\s+", header=None)
    a.columns = ["Indiv"] + list(IndGeno.loc[:, 0])

    refA = a.loc[:, ref]

    meanRef = np.mean(refA).mean()

    RefAmean[number] = [herd1, herd2, meanRef]
    number = number + 1

RefDF = pd.DataFrame.from_dict(RefAmean, orient="index")
RefADF = RefDF.drop_duplicates()
RefADF.columns = ["Herd1", "Herd2", "A"]
RefADF.to_csv("RefADF_mean" + str(start) + ".csv", index=None)