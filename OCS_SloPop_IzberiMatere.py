# -*- coding: utf-8 -*-
import ast
import random
from itertools import chain

WorkingDir = "/home/jana/"
os.chdir(WorkingDir)
scenario = "Gen"
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


os.chdir("/home/jana/bin/AlphaMateLinux/OCSSloPop/CowSample/")
#GenPed_EBV.txt in PedigreeAndGeneticValues_cat.txt sta ista pedigreja- samo z rauličnimi informacijami
pedCat = pd.read_table("PedigreeAndGeneticValues_cat.txt", sep=" ")
#tukaj odberi metere
ped, (bm, motherOther) = odbira_mater("/home/jana/bin/AlphaMateLinux/OCSSloPop/GenPed_EBV.txt", **selPar)
#tukaj odberi random sample krav
motherSample = random.sample(motherOther, int(len(motherOther)*0.2)) + random.sample(bm, int(len(bm)*0.2))

nMating = int(len(motherOther)*0.2 + int(len(bm)*0.2))

#tukaj spusti čez skript, da ti odbere plemenske bike --> TO JE ZA PRVO LETO PO PREHODU, KO GREJO MLADI IN ČAKAJOČI V GPB
#mora pridit nujno za mamami, ker doopolni pedigre!
#ped, (genTest, (classOce, genOce)) = odberi_testOce(ped, **selPar)
#tukaj pridobi potomceNP pred odbiro
#potomciNP = pedCat.Indiv[(pedCat.cat == "potomciNP") & (pedCat.sex == "M")] #če izbereš potomceNP, moške, tukaj --> ti gredo vsi v genomsko testiranje!
#združi kandidate za optimizacijo
#potOce = list(set(list(potomciNP) + list(genOce)))

#PRVO LETO ODBEREŠ SAMO GENOMSKO TESTIRANE - plus gestirane (kandidati)
#pri testiranih daš tastarejše stran
#potOce = list(genTest) + list(genOce)
#select bulls + candidates + random sample of cows




#PRVO LETO daš v optimizacijo (pred odbiro!!!) vhljevljene, mlade in čakajoče
#naslednja leta dajes "osnovo" in genTest (pred odbiro!!!)
#nato spremeniš kategorije v gpb in semeniš z izbranimi (= kot pa prvo VSE - brez odbire - prestaviš v gpb in potem optimiziraš gpb)
potOce = pedCat.Indiv[pedCat.cat.isin(["vhlevljeni", "mladi", "cakajoci"])]

pd.DataFrame({"ID": list(motherSample) + list(potOce)}).to_csv("/home/jana/bin/AlphaMateLinux/OCSSloPop/CowSample/IndOpt.txt", index=None, header=None)



Cont = pd.read_table("/home/jana/bin/AlphaMateLinux/OCSSloPop/CowSample/ContributionsModeOptTarget1.txt", sep="\s+")
selM = Cont.loc[Cont.Gender==1]
selM.loc[:, "nMatingPop"] = (selM.loc[:,"nMating"] / 0.2).astype(int)
Ocetje = list(chain.from_iterable([[int(father)] * dose for (father, dose) in zip(selM.Id, selM.nMatingPop)]))
shuffle(Ocetje)

#dodaj novo generacijo
ped.add_new_gen_naive(selPar['stNBn'], selPar['potomciNPn'] * 2)
ped.compute_age()
# dodaj matere
ped.doloci_matere(kwargs.get('stNBn'), kwargs.get('potomciNPn'), kwargs.get('ptn'), kwargs.get('kraveUp'))
# preveri - mora biti nič!!! - oz. če mater še ni dovolj, potem še ne!
ped.mother_nr_blank()
# dodaj očete
ped.ped.loc[ped.ped.cat.isin(['nr', 'potomciNP']), 'Father'] = Ocetje
ped.ped.Father = ped.ped.Father.astype(int)
ped.write_ped("/home/jana/bin/AlphaMateLinux/OCSSloPop/CowSample//ExternalPedigree.txt")
ped.write_pedTotal("/home/jana/bin/AlphaMateLinux/OCSSloPop/CowSample//ExternalPedigreeTotal.txt")
