for str in "SU15"
do
	for scenario in "Class" "GenSLO" "OtherCowsGen" "BmGen" Gen"
	do
		for rep in $(seq 0 19)
		do
			cd /home/v1jobste/JanaO/$str/${scenario}${rep}/SimulatedData/AllIndividualsSnpChips/
			cp /home/v1jobste/JanaO/Essentials/Inbreeding_Ind* .
			grep -Fwf Inbreeding_Individuals.txt Chip1Genotype.txt > Chip1Genotype_Last20Gen.txt
			rm Chip2*
			rm Chip1Genotype.txt
			zip Chip1.zip Chip1Genotype_Last20Gen.txt
		done
	done
done
