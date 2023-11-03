for rep in $(seq 13 13)
	do
	for degree in  15  45 60 75 ## 75 ##50 55 ##15 30 45 60 75
	do
		sed "s/DEGREE/$degree/g" SimulationEddie_0GenOCS_SU55_new.sh > SimulationEddie_${degree}GenOCS.sh
		sed "s/_REPEAT_/$rep/g" SimulationEddie_${degree}GenOCS.sh > SimulationEddie_${degree}_${rep}GenOCS.sh
		qsub -pe sharedmem 8 -l h_vmem=4G SimulationEddie_${degree}_${rep}GenOCS.sh
	done
done

