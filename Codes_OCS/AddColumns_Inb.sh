for rep in $(seq 0 10)
do
	for scenario in "Class" "GenSLO" "OtherCowsGen" "BmGen" "Gen"
	do
		awk -v var1=${rep}_${scenario} '{print $0, var1}' Inbreeding_${scenario}${rep}.txt > Inbreeding_${scenario}${rep}A.txt
	done
done
