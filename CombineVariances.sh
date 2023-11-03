for scenario in "Class" 
do
	for rep in $(seq 0 9)
	do
		awk -v temp1=$rep -v temp3=$scenario 'BEGIN{OFS=IFS="\t"} {print $0, temp1,temp3}' ${scenario}${rep}/SimulatedData/TotalGenicAndGeneticVariancesPerGeneration.txt > tmpVariance.txt
		cat VarianceFile.txt tmpVariance.txt > tmp1 && mv tmp1 VarianceFile.txt
	done
done
