from selection10 import genInterval
import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt
import os




genInt = genInterval(os.getcwd() + '/')
genInt.plotGenInt()

print 'Created plot GenInt_PLOT.png'
