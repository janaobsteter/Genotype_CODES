import pandas as pd
import os
import sys

startRep = int(sys.argv[1])
numRep = int(sys.argv[2])


par = pd.read_csv("ParameterFile_Simulation.csv")


homeDir = os.getcwd()
SelParDir = homeDir + "/SelPar/10K_Pheno/"

for index, row in par.iterrows():

    name, ref, pg, nocontrol, potomciNPpy, telMpy, femalepy, totalpy, malesel, update, diffGenoSize, yearGeno = list(row)
    print(name)
    print(yearGeno)

    if ref:
	pass
        #change the variables in the SelectionParam file --> write it to SelectionParam_Name to current directory
        os.system("cp " + SelParDir + "/SelectionParam_Gen_MaleGS.csv " + SelParDir + "/SelectionParam_Gen_MaleGS_" + name + ".csv")

        #change NumberCowsGenotyped and NumberMalesGenotyped (male newborns/calves)
        os.system("sed -i 's/NumberPotomciNPGenotyped/" + str(potomciNPpy) + "/g' " + SelParDir + "/SelectionParam_Gen_MaleGS_" + name + ".csv")
        os.system("sed -i 's/NumberTelMGenotyped/" + str(telMpy) + "/g' " + SelParDir + "/SelectionParam_Gen_MaleGS_" + name + ".csv")
        os.system("sed -i 's/NumberCowsGenotyped/" + str(femalepy) + "/g' " + SelParDir + "/SelectionParam_Gen_MaleGS_" + name + ".csv")
        #change UpdateOrNot, MaleSelection
        os.system("sed -i 's/UpdateOrNot/" + str(update) + "/g' " + SelParDir + "/SelectionParam_Gen_MaleGS_" + name + ".csv")
        os.system("sed -i 's/MaleSelection/" + str(malesel) + "/g' " + SelParDir + "/SelectionParam_Gen_MaleGS_" + name + ".csv")
        os.system("sed -i 's/yearStartGeno/" + str(yearGeno) + "/g' " + SelParDir + "/SelectionParam_Gen_MaleGS_" + name + ".csv")



    elif not ref:
        # change the variables in the SelectionParam file --> write it to SelectionParam_Name to current directory
        os.system("cp " + SelParDir + "/SelectionParam_Gen_MaleGS_noRef.csv " + SelParDir + "/SelectionParam_Gen_MaleGS_" + name + ".csv")
        # change NumberCowsGenotyped and NumberMalesGenotyped (male newborns/calves)
        os.system("sed -i 's/NumberCowsGenotypedInit/" + str(totalpy) + "/g' " + SelParDir + "/SelectionParam_Gen_MaleGS_" + name + ".csv")
        os.system("sed -i 's/NumberPotomciNPGenotyped/" + str(potomciNPpy) + "/g' " + SelParDir + "/SelectionParam_Gen_MaleGS_" + name + ".csv")
        os.system("sed -i 's/NumberTelMGenotyped/" + str(telMpy) + "/g' " + SelParDir + "/SelectionParam_Gen_MaleGS_" + name + ".csv")
        os.system("sed -i 's/NumberCowsGenotyped/" + str(femalepy) + "/g' " + SelParDir + "/SelectionParam_Gen_MaleGS_" + name + ".csv")
        # change UpdateOrNot, MaleSelection
        os.system("sed -i 's/UpdateOrNot/" + str(update) + "/g' " + SelParDir + "/SelectionParam_Gen_MaleGS_" + name + ".csv")
        os.system("sed -i 's/MaleSelection/" + str(malesel) + "/g' " + SelParDir + "/SelectionParam_Gen_MaleGS_" + name + ".csv")
        os.system("sed -i 's/SizeStartGeno/" + str(diffGenoSize) + "/g' " + SelParDir + "/SelectionParam_Gen_MaleGS_" + name + ".csv")
        os.system("sed -i 's/yearStartGeno/" + str(yearGeno) + "/g' " + SelParDir + "/SelectionParam_Gen_MaleGS_" + name + ".csv")



    #change the Eddie.py script to read the SelectionParam_Name.csv
    #change the variables in the qsub script
    #for each rep
    for rep in range(startRep, numRep):
        os.system("cp " + homeDir + "/SimulationEddie_pheno_generic.sh " + homeDir + "/SimulationEddie_pheno_" + str(rep) + "_" + name + ".sh")
        os.system("sed -i 's/_REP_/" + str(rep) + "/g' " + homeDir + "/SimulationEddie_pheno_" + str(rep) + "_" + name + ".sh")
        os.system("sed -i 's/_NAME_/" + str(name) + "/g' " + homeDir + "/SimulationEddie_pheno_" + str(rep) + "_" + name + ".sh")
        os.system("sed -i 's/_REPEATS_/" + str(nocontrol) + "/g' " + homeDir + "/SimulationEddie_pheno_" + str(rep) + "_" + name + ".sh")
        os.system("sed -i 's/_REFERENCE_/" + str("10K" if ref else 0) + "/g' " + homeDir + "/SimulationEddie_pheno_" + str(rep) + "_" + name + ".sh")



        #change the Eddie.py script to make the directory of Gen_Name
        #add Name to qsub file

        #qsub the file
#        os.system("qsub " + homeDir + "/SimulationEddie_pheno_" + str(rep) + "_" + name + ".sh")

