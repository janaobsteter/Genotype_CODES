cats = ['gpb']

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

os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Gen/')
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






