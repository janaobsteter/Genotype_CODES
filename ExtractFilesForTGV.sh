#mkdir MissingFiles
scenario="ClassGen"
for rep in 5
do
    for import in 0_0 10_10 25_25 50_50 100_100 0_100
    do
        for trait in 2 3
        do
            mkdir MissingFiles/${scenario}${rep}_${import}${trait}
            tar -zxvf ${scenario}${rep}${import}${trait}.tar.gz home/v1jobste/jobsteter/10K/SU55_import/${scenario}${rep}_${import}1${trait}/PopulationSplit.txt
            mv home/v1jobste/jobsteter/10K/SU55_import/${scenario}${rep}_${import}1${trait}/PopulationSplit.txt MissingFiles/${scenario}${rep}_${import}${trait}/PopulationSplit.txt
#            tar -zxvf ClassGen${rep}${import}${trait}.tar.gz SimulatedData/PedigreeAndGeneticValues.txt
#            mv PedigreeAndGeneticValues.txt MissingFiles/ClassGen${rep}_${import}${trait}/PedigreeAndGeneticValues.txt
#            tar -xvf ClassGen${rep}${import}${trait}.tar.gz GenicVariance_import.csv
#            mv GenicVariance_import.csv MissingFiles/ClassGen${rep}_${import}${trait}/GenicVariance_import.csv
        done
    done
done
	
