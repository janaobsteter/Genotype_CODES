for scenario in "Class" "GenSLO" "OtherCowsGen" "BmGen" "Gen"
do
	for rep in $(seq 0 19)
	do
		gen=$(echo $(tail ${strategy}/${scenario}${rep}/SimulatedData/PedigreeAndGeneticValues.txt | cut -f9 -d" " | sort | uniq))
		echo ${strategy}_$scenario$rep
		echo $gen
	done
done
