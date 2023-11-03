for str in "SU55"
do
	for scenario in "Class" "GenSLO" "OtherCowsGen" "BmGen" Gen"
	do
		for rep in $(seq 0 19)
		do
			cd /home/v1jobste/JanaO/$str/${scenario}${rep}/SimulatedData/
		        files=$(find . -maxdepth 1 -type f | grep -v UnrestrictedQtnIndivGenotypes.txt)	
			tar -czvf SimulatedData.tar.gz $files
			rm $files			
		done
	done
done
