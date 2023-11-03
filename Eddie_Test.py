# -*- coding: utf-8 -*-
from __future__ import division
import sys

rep = sys.argv[1]
scenario = sys.argv[2]
strategy = sys.argv[3]
refSize = sys.argv[4]
repeats = int(sys.argv[5])
variances = sys.argv[6].split(",")
varPE = float(variances[0])
varH = float(variances[1])
varHY = float(variances[2])
varHTD = float(variances[3])
varE = float(variances[4])
name = sys.argv[7]



print("Rep in " + str(rep))
print("Repeats is " + str(repeats))
print("RefSize is " + str(refSize))
print(variances)

if refSize  != "0":
    print('RefSize is not 0!')
elif refSize == "0":
    print("refSize is 0")

