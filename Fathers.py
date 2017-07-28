pedT = pd.read_csv(AlphaSimDir +"/SimulatedData/PedigreeAndGeneticValues_cat.txt", sep=" ")
ped = pedT.loc[pedT.Generation.isin(range(58))]
pb = ped.loc[ped.cat=='gpb']
pb.loc[:, "NbOf"] = [len(ped.loc[(ped.Father == i)]) for i in list(pb.Indiv)] # & (ped.sex=='F')

plt.plot(pb.Generation, pb.NbOf)
plt.xlabel('Leto rojstva bika')
plt.ylabel('Stevilo potomcev')


counts =  pd.DataFrame({'Father': list(ped.loc[ped.Father, "Father"].value_counts().index), 'Offspring': list(ped.loc[ped.cat==cat, "Father"].value_counts())})
counts =  pd.DataFrame({'Father': list(ped.Father.value_counts().index), 'Offspring': list(ped.Father.value_counts())})

cakOcetjeBest = list(chain.from_iterable([Ped.catCurrent_indiv_age_CriteriaEBV('cak', (2 + x), 4) for x in range(1, cak + 1)]))


#KOLIKO potomcev imajo pb na tranziciji
ped = pedT.loc[pedT.Generation.isin(range(43,44))]
FatT = pd.DataFrame({"Indiv": list(set(ped.Father))})
FatT = pd.merge(FatT, pedT, on="Indiv")
FatT.loc[:, "NbOf"] = [len(ped.loc[(ped.Father == i)]) for i in list(FatT.Indiv)] # & (ped.sex=='F')
FatT[FatT.cat=='gpb']