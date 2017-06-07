ped = pd.read_table(AlphaSimDir + '/Blupf90.ped', header=None, sep=" ")
gens = list(chain.from_iterable([[i] * 100 for i in range(1970, 1999)]))

ped.loc[:, 'gen'] = gens
ped[[0,1,2,'gen']].to_csv(AlphaSimDir + '/PedYOB.txt', header=None, index=None, sep=" ")


dat = pd.read_table(AlphaSimDir + '/blupDat', header=None, sep=" ")
dat.loc[:, 'gen'] = gens
dat[[1,2,3,'gen']].to_csv(AlphaSimDir + '/DatYOB.txt', header=None, index=None, sep=" ")
