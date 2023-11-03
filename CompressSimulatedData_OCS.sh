for str in "OCS"
do
	for scenario in 15 30 45 50 55 60 75
	do
		for rep in $(seq 0 18)
		do
			cd /home/v1jobste/JanaO/Gen${rep}_${scenario}OCS/SimulatedData/
			FILE=SimulatedData.tar.gz
			if [ ! -f "$FILE" ]; then
				files=$(find . -maxdepth 1 -type f | grep -v UnrestrictedQtnIndivGenotypes.txt)	
				tar -czvf SimulatedData.tar.gz $files
				rm $files			
			fi
		done
	done
done
