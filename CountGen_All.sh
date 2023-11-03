for scenario in "15" "30" "45" "60" "75"
do
	for rep in $(seq 10 19)
	do
		gen=$(echo $(tail Gen${rep}_${scenario}OCS/SimulatedData/PedigreeAndGeneticValues.txt | cut -f9 -d" " | sort | uniq))
		gen1=$(echo $(tail Gen${rep}_${scenario}OCS/GenPed_EBV.txt | cut -f1 -d"," | sort | uniq))
		echo ${rep}_${scenario} 
		echo $gen
		if [ $gen1 -ne 60 ]
			then
			echo ${rep}_${scenario}
			echo $gen1
		fi
	done
done
