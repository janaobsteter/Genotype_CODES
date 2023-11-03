for strategy in "SU52" "SU53" "SU54"
do
	for scenario in "Class" "Gen"  ##"Class" "GenSLO" "OtherCowsGen" "BmGen" "Gen"
	do
		for rep in $(seq 0 19)
		do
			gen=$(echo $(tail 10K/${strategy}/${scenario}${rep}/SimulatedData/PedigreeAndGeneticValues.txt | cut -f9 -d" " | sort | uniq))
			echo ${strategy}_$scenario$rep
			echo $gen
		done
	done
done
