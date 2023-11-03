for str in "SU55" 
do
	for scenario in "Class" "GenSLO" "OtherCowsGen" "BmGen" "Gen"
	do
		for rep in $(seq 0 19)
		do
			echo $str$scenario$rep
			ls ${str}/${scenario}${rep}/SimulatedData/QTNGenotype_Last20Gen.txt
		done
	done
done
