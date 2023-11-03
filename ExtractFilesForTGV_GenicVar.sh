#mkdir MissingFiles
scenario="GenGen"
for rep in 7
do
    for import in 0_100
    do
        for trait in 2 3
        do
#            mkdir MissingFiles/ClassGen${rep}_${import}${trait}
#            tar -zxvf ClassGen${rep}${import}${trait}.tar.gz home/v1jobste/jobsteter/10K/SU55_import/ClassGen${rep}_${import}1${trait}/PopulationSplit.txt
#            mv home/v1jobste/jobsteter/10K/SU55_import/ClassGen${rep}_${import}1${trait}/PopulationSplit.txt MissingFiles/ClassGen${rep}_${import}${trait}/PopulationSplit.txt
#            tar -zxvf ClassGen${rep}${import}${trait}.tar.gz SimulatedData/PedigreeAndGeneticValues.txt
#            mv PedigreeAndGeneticValues.txt MissingFiles/ClassGen${rep}_${import}${trait}/PedigreeAndGeneticValues.txt
            tar -xvf ${scenario}${rep}${import}${trait}.tar.gz home/v1jobste/jobsteter/10K/SU55_import/${scenario}${rep}_${import}1${trait}/GenicVariance_import.csv
            mv home/v1jobste/jobsteter/10K/SU55_import/${scenario}${rep}_${import}1${trait}/GenicVariance_import.csv MissingFiles/${scenario}${rep}_${import}${trait}/GenicVariance_import.csv
        done
    done
done
	
