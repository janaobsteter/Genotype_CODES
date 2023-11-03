for rep in $(seq 11 15)
	do
	for degree in  15 ##55 ##15 45 60 75 ##30
	do
		sed "s/DEGREE/$degree/g" SimulationEddie_0GenOCS_SU55_drugaPol.sh > SimulationEddie_${degree}GenOCS.sh
		sed "s/_REPEAT_/$rep/g" SimulationEddie_${degree}GenOCS.sh > SimulationEddie_${degree}_${rep}GenOCS_drugaPol.sh
		qsub -pe sharedmem 8 -l h_vmem=4G SimulationEddie_${degree}_${rep}GenOCS_drugaPol.sh
	done
done

