for str in "SU55" "SU15" "SU51"
do
	for scenario in "OtherCowsGen"
	do
		for rep in $(seq 0 19)
		do
			cd /home/v1jobste/JanaO/$str/${scenario}${rep}/SimulatedData/
			cp /home/v1jobste/JanaO/Essentials/Inbreeding_Ind* .
			grep -Fwf Inbreeding_Individuals.txt UnrestrictedQtnIndivGenotypes.txt > Last20Gen_QTNs.txt
		done
	done
done
