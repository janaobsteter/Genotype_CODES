# -*- coding: utf-8 -*-
import ast
WorkingDir = "/home/jana/"
os.chdir(WorkingDir)
scenario = "GenSLO"
par = pd.read_csv(WorkingDir + "/Essentials/SelectionParam_" + scenario + ".csv", header=None, names=["Keys", "Vals"])
par.to_dict()
selPar = defaultdict()
for key, val in zip(par.Keys, par.Vals):
    if key not in ['BurnInYN', 'EBV', 'gEBV', 'PA', 'AlphaSimDir', 'genotyped', 'EliteDamsPTBulls',
                   'EliteDamsPABulls', 'UpdateGenRef', 'sexToUpdate', 'EliteDamsGenBulls', 'gpb_pb',
                   'genTest_mladi', 'genTest_gpb']:
        try:
            selPar[key] = int(val)
        except:
            selPar[key] = float(val)
    if key in ['BurnInYN', 'EBV', 'gEBV', 'PA', 'AlphaSimDir', 'EliteDamsPTBulls',
               'EliteDamsPABulls', 'UpdateGenRef', 'sexToUpdate', 'EliteDamsGenBulls', 'gpb_pb',
               'genTest_mladi', 'genTest_gpb']:
        if val in ['False', 'True']:
            selPar[key] = bool(val == 'True')
        else:
            selPar[key] = val
    if key == 'genotyped':
        selPar[key] = ast.literal_eval(val)


BurnInYN = "False" #ali izvedeš tudi BurnIn
SelYN = "True" #ali izvedeš tudi BurnIn
StNB = 8640
StBurnInGen = 20
StFillInBurnIn = 40
StSelGen = 40
StartSelGen = 21
StopSelGen = 40
NumberOfSires = 12
NumberOfDams = 3500
selPar['AlphaSimDir'] = os.getcwd() + '/'
AlphaSimDir = os.getcwd() + '/'
AlphaSimPed = selPar['AlphaSimDir'] + '/SimulatedData/PedigreeAndGeneticValues.txt'
if selPar['EBV']:
    seltype = 'class'
if selPar['gEBV']:
    seltype = 'gen'

selPar['pbn'] = 1

selPar['genpbn'] = 1
selPar['pozitivnoTestDoz'] = 2000

selPar['pbUp'] = 1

IndCat = pd.DataFrame()
ped0, c0, s0, a0 = selekcija_total("/home/jana/GenPed_EBV.txt", **selPar)

inds = ped0.catCurrent_indiv('')
Cats = []
for ind in inds:
    for (key, value) in categories.iteritems():
        if ind in value:
            print(key)
            
for ind in inds:
    for (key, value) in categories.iteritems():
        if ind in value:
            print (key)



IndCat['Indiv'] = ped.ped.Indiv
IndCat['catBurnIN'] = ped.ped.cat


for krog in range(krogov): #ponavljaj kolikor krogov selekcije hočeš
    ped, c, s, a = selekcija_total("/home/jana/PedTotal.txt", **selPar)
    IndCat[str('cat' + str(krog))] = ped.ped.cat
    
    
    eclf1 = VotingClassifier(estimators=[('lr1', CellWall), ('lr2', Vacuoles)],voting='hard'))