for strategy in "ClassGen" "GenGen"
do
	for scenario in "100_100" "50_50" "25_25" "10_10" "0_100" "0_0"
	do
		for rep in $(seq 11 12)
		do
			for trait in 2 3
			do
				gen=$(echo $(tail 10K/SU55_import/${strategy}${rep}_${scenario}1${trait}/SimulatedData/PedigreeAndGeneticValues.txt | cut -f9 -d" " | sort | uniq))
				echo ${strategy}_${scenario}_${rep}_${trait}
				echo $gen
			done
		done
	done
done
