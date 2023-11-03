HOMEDIR=$PWD
for scenario in "Class" "GenSLO" "OtherCowsGen" "BmGen" "Gen"
do
	for rep in $(seq 0 19)
	do 
		cd ${scenario}${rep}/SimulatedData/
		tar -xf SimulatedData.tar.gz ./TotalGenicAndGeneticVariancesPerGeneration.txt
		cd $HOMEDIR
	done
done
