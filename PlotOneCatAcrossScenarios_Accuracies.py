cats = ['pb']

os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Class/')
accuraciesEBV = pd.read_table('AccuraciesEBVPA.csv', sep=",")
accuraciesEBV = accuraciesEBV[accuraciesEBV.Cycle.isin(range(40, 61))]
for cat in cats:
    trendCat = accuraciesEBV.loc[accuraciesEBV.Cat == cat]
    plt.plot(trendCat.Cycle, trendCat.corEBV, label = 'genClass')
plt.xlabel('Selected Generation')
plt.ylabel('Mean Generation TBV')
legend(loc='upper left')
plt.xticks(accuraciesEBV.Cycle)
plt.title(cats)

cats='genTest'
os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenSplosnaPop/')
accuraciesEBV = pd.read_table('AccuraciesEBVPA.csv', sep=",")
accuraciesEBV = accuraciesEBV[accuraciesEBV.Cycle.isin(range(40, 61))]
for cat in cats:
    trendCat = accuraciesEBV.loc[accuraciesEBV.Cat == cat]
    plt.plot(trendCat.Cycle, trendCat.corEBV, label = 'corTGV_EBV')
    plt.plot(trendCat.Cycle, trendCat.corPA, label = 'corTGV_PA')
plt.xlabel('Selected Generation')
plt.ylabel('Mean Generation TBV')
legend(loc='upper left')
plt.xticks(accuraciesEBV.Cycle)
plt.title('genTest')



os.chdir('/home/jana/bin/AlphaSim1.05Linux//')
accuraciesEBV = pd.read_table('AccuraciesEBVPA.csv', sep=",")
accuraciesEBV = accuraciesEBV[accuraciesEBV.Cycle.isin(range(40, 61))]
for cat in cats:
    trendCat = accuraciesEBV.loc[accuraciesEBV.Cat == cat]
    plt.plot(trendCat.Cycle, trendCat.corEBV, label = 'genNoRef')
plt.xlabel('Selected Generation')
plt.ylabel('Mean Generation TBV')
legend(loc='upper left')


os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenReference/')
accuraciesEBV = pd.read_table('AccuraciesEBVPA.csv', sep=",")
accuraciesEBV = accuraciesEBV[accuraciesEBV.Cycle.isin(range(40, 61))]
for cat in cats:
    trendCat = accuraciesEBV.loc[accuraciesEBV.Cat == cat]
    plt.plot(trendCat.Cycle, trendCat.corEBV, label = 'genRef')
plt.xlabel('Selected Generation')
plt.ylabel('Mean Generation TBV')
legend(loc='upper left')

os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenFather/')
accuraciesEBV = pd.read_table('AccuraciesEBVPA.csv', sep=",")
accuraciesEBV = accuraciesEBV[accuraciesEBV.Cycle.isin(range(40, 61))]
for cat in cats:
    trendCat = accuraciesEBV.loc[accuraciesEBV.Cat == cat]
    plt.plot(trendCat.Cycle, trendCat.corEBV, label = 'genFatherNoRef')
plt.xlabel('Selected Generation')
plt.ylabel('Mean Generation TBV')
legend(loc='upper left')

os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenFatherReference/')
accuraciesEBV = pd.read_table('AccuraciesEBVPA.csv', sep=",")
accuraciesEBV = accuraciesEBV[accuraciesEBV.Cycle.isin(range(40, 61))]
for cat in cats:
    trendCat = accuraciesEBV.loc[accuraciesEBV.Cat == cat]
    plt.plot(trendCat.Cycle, trendCat.corEBV, label = 'genFatherRef')
plt.xlabel('Selected Generation')
plt.ylabel('Mean Generation TBV')
legend(loc='upper left')






