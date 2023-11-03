for scenario in "Class" "GenSLO" "BmGen" "OtherCowsGen" "Gen"
do
	for rep in $(seq 0 19)
	do
		echo $scenario$rep
		ls $scenario$rep/*Marker*
	done
done
