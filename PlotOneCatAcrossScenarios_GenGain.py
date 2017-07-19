#PLOT GENETIC GAIN
f, axarr = plt.subplots(nrows=3, ncols=2, sharex=True,sharey=True)
f.suptitle('GenGain')
cats = ['pb', 'mladi', 'pBM', 'gpb', 'bm', 'k']
x = 0
y = 0
r = 0
xs = [0,1,2, 0,1, 2]
ys = [0,0,0,1,1,1]
os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Class/')
accuraciesEBV = pd.read_table('AccuraciesEBVPA.csv', sep=",")
accuraciesEBV = accuraciesEBV[accuraciesEBV.Cycle.isin(range(40, 61))]
for cat in cats:
    trendCat = accuraciesEBV.loc[accuraciesEBV.Cat == cat]
    axarr[xs[r],ys[r]].plot(trendCat.Cycle, trendCat.meanTBV, label = 'genClass')
    axarr[xs[r],ys[r]].set_title(cat)
    r = r+1


"""plt.xlabel('Selected Generation')
plt.ylabel('Mean Generation TBV')
legend(loc='upper left')
plt.xticks(accuraciesEBV.Cycle)
plt.title(cats)
"""
x = 0
y = 0
r=0
os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenSLO/')
accuraciesEBV = pd.read_table('AccuraciesEBVPA.csv', sep=",")
accuraciesEBV = accuraciesEBV[accuraciesEBV.Cycle.isin(range(40, 61))]
for cat in cats:
    trendCat = accuraciesEBV.loc[accuraciesEBV.Cat == cat]
    axarr[xs[r],ys[r]].plot(trendCat.Cycle, trendCat.meanTBV, label = 'genFather_UpdatedRef')
    r = r + 1
    
    
    
    
#PLOT ACCURACIES
f, axarr = plt.subplots(nrows=3, ncols=2, sharex=True,sharey=True)
f.suptitle('GenGain')
cats = ['pb', 'mladi', 'pBM', 'gpb', 'bm', 'k']
x = 0
y = 0
r = 0
xs = [0,1,2, 0,1, 2]
ys = [0,0,0,1,1,1]
os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Class/')
accuraciesEBV = pd.read_table('AccuraciesEBVPA.csv', sep=",")
accuraciesEBV = accuraciesEBV[accuraciesEBV.Cycle.isin(range(40, 61))]
for cat in cats:
    trendCat = accuraciesEBV.loc[accuraciesEBV.Cat == cat]
    axarr[xs[r],ys[r]].plot(trendCat.Cycle, trendCat.corEBV, label = 'genClass')
    axarr[xs[r],ys[r]].set_title(cat)
    r = r+1


"""plt.xlabel('Selected Generation')
plt.ylabel('Mean Generation TBV')
legend(loc='upper left')
plt.xticks(accuraciesEBV.Cycle)
plt.title(cats)
"""
x = 0
y = 0
r=0
os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenSLO/')
accuraciesEBV = pd.read_table('AccuraciesEBVPA.csv', sep=",")
accuraciesEBV = accuraciesEBV[accuraciesEBV.Cycle.isin(range(40, 61))]
for cat in cats:
    trendCat = accuraciesEBV.loc[accuraciesEBV.Cat == cat]
    axarr[xs[r],ys[r]].plot(trendCat.Cycle, trendCat.corEBV, label = 'genFather_UpdatedRef')
    r = r + 1


plt.xlabel('Selected Generation')
plt.ylabel('Mean Generation TBV')
legend(loc='upper left')
plt.xticks(accuraciesEBV.Cycle)
plt.title(cats)



os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenBM_UpdatedRef')

#os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Gen/')
accuraciesEBV = pd.read_table('AccuraciesEBVPA.csv', sep=",")
accuraciesEBV = accuraciesEBV[accuraciesEBV.Cycle.isin(range(40, 61))]
for cat in cats:
    trendCat = accuraciesEBV.loc[accuraciesEBV.Cat == cat]
    plt.plot(trendCat.Cycle, trendCat.meanTBV, label = 'genNoRef')
plt.xlabel('Selected Generation')
plt.ylabel('Mean Generation TBV')
legend(loc='upper left')


os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenReference/')
accuraciesEBV = pd.read_table('AccuraciesEBVPA.csv', sep=",")
accuraciesEBV = accuraciesEBV[accuraciesEBV.Cycle.isin(range(40, 61))]
for cat in cats:
    trendCat = accuraciesEBV.loc[accuraciesEBV.Cat == cat]
    plt.plot(trendCat.Cycle, trendCat.meanTBV, label = 'genRef')
plt.xlabel('Selected Generation')
plt.ylabel('Mean Generation TBV')
legend(loc='upper left')

os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenFather/')
accuraciesEBV = pd.read_table('AccuraciesEBVPA.csv', sep=",")
accuraciesEBV = accuraciesEBV[accuraciesEBV.Cycle.isin(range(40, 61))]
for cat in cats:
    trendCat = accuraciesEBV.loc[accuraciesEBV.Cat == cat]
    plt.plot(trendCat.Cycle, trendCat.meanTBV, label = 'genFatherNoRef')
plt.xlabel('Selected Generation')
plt.ylabel('Mean Generation TBV')
legend(loc='upper left')

os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenFatherReference/')
accuraciesEBV = pd.read_table('AccuraciesEBVPA.csv', sep=",")
accuraciesEBV = accuraciesEBV[accuraciesEBV.Cycle.isin(range(40, 61))]
for cat in cats:
    trendCat = accuraciesEBV.loc[accuraciesEBV.Cat == cat]
    plt.plot(trendCat.Cycle, trendCat.meanTBV, label = 'genFatherRef')
plt.xlabel('Selected Generation')
plt.ylabel('Mean Generation TBV')
legend(loc='upper left')






