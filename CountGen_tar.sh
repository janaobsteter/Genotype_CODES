for scenario in "Class" "GenSLO" "OtherCowsGen" "BmGen" "Gen"
do
	for rep in $(seq 0 19)
	do
		gen=$(tar -x -f ${scenario}${rep}/SimulatedData/SimulatedData.tar.gz -O ./PedigreeAndGeneticValues_cat.txt | tail -1 | cut  -f1 -d" ")
		echo $scenario$rep
		echo $gen
	done
done
