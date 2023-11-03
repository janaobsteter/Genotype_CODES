for str in "SU55" "SU51" "SU15"
do
	for scenario in "Class" "GenSLO" "OtherCowsGen" "BmGen" "Gen"
	do
		for rep in $(seq 0 19)
		do
			echo $str$scenario$rep
			ls ${str}/${scenario}${rep}/*Marker*
		done
	done
done
