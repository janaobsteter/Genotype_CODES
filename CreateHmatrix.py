import numpy as np
import pandas as pd

a = pd.read_table("PedigreeNrm.txt", sep="\s+", header=None)
Aanimals = list(a.loc[:,0])
a = a.drop(0, axis=1)
a.columns = list(Aanimals)
a.index = Aanimals


g = pd.read_table("GenotypeNrm.txt", header=None)
Ganimals = list(g.loc[:, 0])
g = g.drop(0, axis=1)
g.index = Ganimals
g.columns = Ganimals


sorted = sorted(Aanimals)
a = a[sorted]
a = a.reindex([sorted])

a11Animals = list(set(Aanimals) - set(Ganimals))
order = a11Animals + Ganimals
a = a.ix[order, order]

a11 = a.ix[a11Animals , a11Animals]
a12 = a.ix[a11Animals , Ganimals]
a21 = a.ix[Ganimals, a11Animals]
a22 = a.ix[Ganimals, Ganimals]

a22inv = pd.DataFrame(np.linalg.pinv(a22.values), a22.columns, a22.index)

gb = 0.95*g + 0.05*a22
d = gb- a22

upLeft = np.matmul(np.matmul(np.matmul(np.matmul(a12, a22inv), d) , a22inv), a21)
downLeft = np.matmul(np.matmul(d, a22inv), a21)
upRight = np.matmul(np.matmul(a12, a22inv), d)
downRight = d


Up = np.hstack((upLeft, upRight))
Down = np.hstack((downLeft, downRight))
Hadd = np.vstack((Up, Down))

H = a + Hadd
H = H.ix[sorted, sorted]

H.to_csv("Hmatrix.txt", header=None, sep=" ")



