for scenario in "Class" "GenSLO" "OtherCowsGen" "BmGen" "Gen"
do
	for rep in $(seq 0 19)
	do
		rm $scenario$rep/SimulatedData/AllIndividualsSnpChips/Chip1Phase.txt
	done
done
