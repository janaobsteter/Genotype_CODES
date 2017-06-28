import random
import pandas as pd
import os
from itertools import chain

os.chdir('/home/jana/bin/AlphaSim1.05Linux/REALFillIn20BurnIn20/')
pd.DataFrame({0: sorted(list(set(chain.from_iterable([random.sample((list(pd.read_table('Categories_gen' + str(i) + 'DF.csv', sep=",")['k'].
dropna().astype(int))), 2000) for i in range(34, 41)]))))}).to_csv('//home/jana/bin/AlphaSim1.05Linux/IndForGeno.txt', index=None, header=None, sep='\n')