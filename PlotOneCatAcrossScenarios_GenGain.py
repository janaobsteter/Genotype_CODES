# -*- coding: utf-8 -*-
#PLOT GENETIC GAIN
f, axarr = plt.subplots(nrows=3, ncols=2, sharex=True,sharey=True)
f.suptitle('GenGain')
cats = ['pb', 'mladi', 'pBM', 'gpb', 'bm', 'k']
catnames = ["Progeno test. plemenski biki", "Mladi biki", "Bikovske matere", "Gen Test plemenski biki", "Bikovske matere - zadnje leto", "Krave"]
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
    axarr[xs[r],ys[r]].plot(trendCat.Cycle, trendCat.meanTBV, label = 'Conventional')
    axarr[xs[r],ys[r]].set_title(catnames[cats.index(cat)])
    r = r+1



x = 0
y = 0
r=0
os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenSLO/')
accuraciesEBV = pd.read_table('AccuraciesEBVPA.csv', sep=",")
accuraciesEBV = accuraciesEBV[accuraciesEBV.Cycle.isin(range(40, 61))]
for cat in cats:
    trendCat = accuraciesEBV.loc[accuraciesEBV.Cat == cat]
    axarr[xs[r],ys[r]].plot(trendCat.Cycle, trendCat.meanTBV, label = 'Conventional_SLO')
    r = r + 1
    
x = 0
y = 0
r=0
os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenSplosnaPop/')
accuraciesEBV = pd.read_table('AccuraciesEBVPA.csv', sep=",")
accuraciesEBV = accuraciesEBV[accuraciesEBV.Cycle.isin(range(40, 61))]
for cat in cats:
    trendCat = accuraciesEBV.loc[accuraciesEBV.Cat == cat]
    axarr[xs[r],ys[r]].plot(trendCat.Cycle, trendCat.meanTBV, label = 'GenBulls On Other Cows')
    r = r + 1
    
x = 0
y = 0
r=0
os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenSLO_BmGen/')
accuraciesEBV = pd.read_table('AccuraciesEBVPA.csv', sep=",")
accuraciesEBV = accuraciesEBV[accuraciesEBV.Cycle.isin(range(40, 61))]
for cat in cats:
    trendCat = accuraciesEBV.loc[accuraciesEBV.Cat == cat]
    axarr[xs[r],ys[r]].plot(trendCat.Cycle, trendCat.meanTBV, label = 'GenBulls on Bull Dams')
    r = r + 1
      
    
x = 0
y = 0
r=0
os.chdir('/home/jana/bin/AlphaSim1.05Linux/')
accuraciesEBV = pd.read_table('AccuraciesEBVPA.csv', sep=",")
accuraciesEBV = accuraciesEBV[accuraciesEBV.Cycle.isin(range(40, 57))]
for cat in cats:
    trendCat = accuraciesEBV.loc[accuraciesEBV.Cat == cat]
    axarr[xs[r],ys[r]].plot(trendCat.Cycle, trendCat.meanTBV, label = 'GenBulls on All Cows')
    r = r + 1
    
plt.xlabel('Year of birth')
plt.ylabel('Mean Generation TBV')
legend(loc='upper left')
plt.xticks(accuraciesEBV.Cycle)
#plt.title(cats)

###############################################    
cats = ['pb', 'mladi', 'pBM', 'gpb', 'bm', 'k']
catnames = ["Progeno test. plemenski biki", "Mladi biki", "Bikovske matere", "Gen Test plemenski biki", "Bikovske matere - zadnje leto", "Krave"]

#cats=['pb']
#PLOT ACCURACIES
f, axarr = plt.subplots(nrows=3, ncols=2, sharex=True,sharey=True)
f.suptitle('Accuracies of EBVs')
#cat = ['pb', 'mladi', 'genTest', 'gpb', 'potomciNP', 'k']
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
    axarr[xs[r],ys[r]].plot(trendCat.YOB, trendCat.corEBV, label = 'Conventional')
    axarr[xs[r],ys[r]].set_title(catnames[cats.index(cat)])
    r = r+1




x = 0
y = 0
r=0
os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenSLO/')
accuraciesEBV = pd.read_table('AccuraciesEBVPA.csv', sep=",")
accuraciesEBV = accuraciesEBV[accuraciesEBV.Cycle.isin(range(40, 61))]
for cat in cats:
    trendCat = accuraciesEBV.loc[accuraciesEBV.Cat == cat]
    axarr[xs[r],ys[r]].plot(trendCat.YOB, trendCat.corEBV, label = 'Conventional_SLO')
    r = r + 1
    
    


x = 0
y = 0
r=0
os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenSplosnaPop/')
accuraciesEBV = pd.read_table('AccuraciesEBVPA.csv', sep=",")
accuraciesEBV = accuraciesEBV[accuraciesEBV.Cycle.isin(range(40, 61))]
for cat in cats:
    trendCat = accuraciesEBV.loc[accuraciesEBV.Cat == cat]
    axarr[xs[r],ys[r]].plot(trendCat.YOB, trendCat.corEBV, label = 'GenBulls on Other Cows')
    r = r + 1

    
x = 0
y = 0
r=0
os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_GenSLO_BmGen/')
accuraciesEBV = pd.read_table('AccuraciesEBVPA.csv', sep=",")
accuraciesEBV = accuraciesEBV[accuraciesEBV.Cycle.isin(range(40, 61))]
for cat in cats:
    trendCat = accuraciesEBV.loc[accuraciesEBV.Cat == cat]
    axarr[xs[r],ys[r]].plot(trendCat.YOB, trendCat.corEBV, label = 'GenBulls on Bull Dams')
    r = r + 1

x = 0
y = 0
r=0
os.chdir('/home/jana/bin/AlphaSim1.05Linux//')
accuraciesEBV = pd.read_table('AccuraciesEBVPA.csv', sep=",")
accuraciesEBV = accuraciesEBV[accuraciesEBV.Cycle.isin(range(40, 57))]
for cat in cats:
    trendCat = accuraciesEBV.loc[accuraciesEBV.Cat == cat]
    axarr[xs[r],ys[r]].plot(trendCat.YOB, trendCat.corEBV, label = 'GenBulls on All Cows')
    r = r + 1

plt.xlabel('Year of Birth')
plt.ylabel('Mean Generation TBV')
legend(loc='upper left')
plt.xticks(accuraciesEBV.YOB)
#plt.title(cats)




####################################################################
#več kategorij - en scenarij - en plot!
os.chdir('/home/jana/bin/AlphaSim1.05Linux/')

#os.chdir('/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_Gen/')
accuraciesEBV = pd.read_table('AccuraciesEBVPA.csv', sep=",")
accuraciesEBV = accuraciesEBV[accuraciesEBV.Cycle.isin(range(40, 58))]
cats=[ 'gpb']
for cat in cats:
    trendCat = accuraciesEBV.loc[accuraciesEBV.Cat == cat]
    plt.plot(trendCat.Cycle, trendCat.corEBV, label = "Leto odbire")
accuraciesEBV = pd.read_table('AccuraciesEBVPA_min.csv', sep=",")
accuraciesEBV = accuraciesEBV[accuraciesEBV.Cycle.isin(range(40, 58))]

for cat in cats:
    trendCat = accuraciesEBV.loc[accuraciesEBV.Cat == cat]
    plt.plot(trendCat.Cycle, trendCat.corEBV, label = "Zadnje leto uporabe")
plt.xlabel('YOB')
plt.ylabel('Accuracy of EBV (cor)')
plt.title("EBV Accuracy - genomically tested Bulls")
legend(loc='lower left')




