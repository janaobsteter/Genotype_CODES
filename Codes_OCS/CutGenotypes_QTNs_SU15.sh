for str in "SU15"
do
	for scenario in "Class" "GenSLO" "OtherCowsGen" "BmGen" "Gen"
	do
		for rep in $(seq 0 19)
		do
			cd /home/v1jobste/JanaO/$str/${scenario}${rep}/SimulatedData/
			cp /home/v1jobste/JanaO/Essentials/Inbreeding_Ind* .
			grep -Fwf Inbreeding_Individuals.txt UnrestrictedQtnIndivGenotypes.txt > Last20Gen_QTNs.txt
		done
	done
done
