ped = pd.read_csv(AlphaSimDir +"/SimulatedData/PedigreeAndGeneticValues_cat.txt", sep=" ")
pb = ped.loc[ped.cat=='mladi']
pb.loc[:, "NbOf"] = [len(ped.loc[(ped.Father == i) & (ped.sex=='F')]) for i in list(pb.Indiv)]

plt.plot(pb.Generation, pb.NbOf)
plt.xlabel('Leto rojstva bika')
plt.ylabel('Stevilo potomcev')


counts =  pd.DataFrame({'Father': list(ped.loc[ped.Father, "Father"].value_counts().index), 'Offspring': list(ped.loc[ped.cat==cat, "Father"].value_counts())})
counts =  pd.DataFrame({'Father': list(ped.Father.value_counts().index), 'Offspring': list(ped.Father.value_counts())})

cakOcetjeBest = list(chain.from_iterable([Ped.catCurrent_indiv_age_CriteriaEBV('cak', (2 + x), 4) for x in range(1, cak + 1)]))
