home="/home/v1jobste/JanaO/"
for rep in $(seq 14 14)
	do
	for degree in 15 30  45 60 75 ## 75 ##50 55 ##15 30 45 60 75
	do
		mkdir Gen${rep}_${degree}OCS
		cd Gen${rep}_${degree}OCS
		cp ${home}/*new1* .
		cp ${home}/selection10.py .
		sed "s/DEGREE/$degree/g" SimulationEddie_0GenOCS_SU55_new1.sh > SimulationEddie_${degree}GenOCS.sh
		sed "s/_REPEAT_/$rep/g" SimulationEddie_${degree}GenOCS.sh > SimulationEddie_${degree}_${rep}GenOCS.sh
		qsub -pe sharedmem 8 -l h_vmem=4G SimulationEddie_${degree}_${rep}GenOCS.sh
		cd ${home}
	done
done

